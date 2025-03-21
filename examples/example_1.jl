using MembraneRD

#here define two examples of measures
	
using ProgressMeter, JLD, Random

struct Measure
	time::Float64
    phi_av::Float64
    binder_cumulant::Float64
    cytoEA::Float64
    cytoEB::Float64
    rho::Float64
end

#make measures here
function Measurer(M::Model; times, name, Nsave)
    isempty(times) && return
	measures = Measure[]
    Nc = length(M)
	next = 1
    function take_measure(t, s)
        #calculate phi at each cell (normalized, goes from -1, to 1)
        ϕ(i) = (s.nB[i] - s.nA[i])/(s.nA[i] + s.nB[i] + 1e-10)
        while next <= lastindex(times) && times[next] <= t
            totA, totB = sum(s.nA), sum(s.nB)
            #measure phi
            ϕav = (totB - totA) / (totA + totB)
            #in some denominators I added + 10^-10 in order to avoid NaNs
            m2 = sum((ϕ(i)-ϕav)^2 for i in 1:Nc) / Nc
            m4 = sum((ϕ(i)-ϕav)^4 for i in 1:Nc) / Nc
            #measure Binder cumulant
            bc = 1 - (m4 / (m2 ^ 2 + 1e-10)) / 3			
            #the following is to check for convergence of rho to 1
            rho = M.rho_0 * ((M.kAd / M.kAa) + totA) / ((M.kBd / M.kBa) + totB)
            push!(measures, Measure(t,ϕav, bc, s.cytoEA[], s.cytoEB[], rho))
			if next % Nsave ==0
				println("T = $(times[next]) and <ϕ>/c = $(ϕav)")
				println("Binder cumulant: $bc")
				save("$name/config_T=$(round(times[next],digits=3)).jld",compress=true, "state", s)
				save("$name/measures.jld",compress=true, "measures", measures)
			end
            next += 1
        end
    end
end

function ProgressMeasurer(T)
    p = Progress(100)
    oldval = 0
    function stat(t, _)
        newval = floor(Int, 100*t/T)
        if newval > oldval
            update!(p, newval)
            oldval = newval
        end
    end
end

function build_model_state(L; rng = Random.default_rng())
    g,posx,posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    memb0 = zeros(Int,N,4)	
    unit_length = 1.0

    theta = 0.2126944621086619#Stot/V


    Ac = unit_length^2#area associated to each cell, set unit lengthscale
    Stot = N*Ac #area of each cell is assumed to be 1, setting unit lengthscale
    V = Stot/theta
    d_timescale = 0.01#this sets unit timescale
    dA = d_timescale
    dB = d_timescale
    dEA = d_timescale
    dEB = d_timescale
    #rates in the theory (do not correspond excatly to those used in the simulations), some slight dimensional changes are needed
    kAa_th = 1.0
    kAd_th = 1.0*kAa_th
    kAc_th = 1.0
    kBa_th = 1.0
    kBd_th = 1.0*kBa_th
    kBc_th = 1.3*kAc_th
    KMM_th = 1.0
    #rates to implement
    kAc = kAc_th
    kBc = kBc_th
    kAa = kAa_th/V
    kAd = kAd_th
    kBa = kBa_th/V
    kBd = kBd_th
    KMM = KMM_th*Ac

    M = Model(g, posx, posy, dA, dB, dEA, dEB, kAc, kBc, kAa, kAd, kBa, kBd, KMM, 0.0)

    N = length(M)
    totmol = N * 10
    totA = 0.5 * totmol
    totB = 0.5 * totmol
    theta = 0.2126944621086619#Stot/V

    nA, nB, nEA, nEB = fill(0,N), fill(0,N), fill(0,N), fill(0,N)

	for i in 1:floor(Int,totA)
		nA[rand(rng, 1:N)]+=1
	end
	for i in 1:floor(Int,totB)
		nB[rand(rng, 1:N)]+=1
	end
    #nA[rand(1:N, floor(Int, totA))] .+= 1
    #nB[rand(1:N, floor(Int, totB))] .+= 1	
    EA_tot_n = floor(Int, 0.1*N)
    EB_tot_n = floor(Int, 0.1*N)
    EA_mem = floor(Int, EA_tot_n*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+((theta/(kAd_th/kAa_th))*(totA/Stot))))
    EB_mem = floor(Int, EB_tot_n*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+((theta/(kBd_th/kBa_th))*(totB/Stot))))
    #nEA[rand(1:N, EA_mem)] .+= 1
    #nEB[rand(1:N, EB_mem)] .+= 1
    for i in 1:EA_mem
		nEA[rand(rng, 1:N)]+=1
	end
	for i in 1:EB_mem
		nEB[rand(rng, 1:N)]+=1
	end
    cytoEA = EA_tot_n - EA_mem
    cytoEB = EB_tot_n - EB_mem
    s = State(nA, nB, nEA, nEB, Ref(cytoEA), Ref(cytoEB))
    M,s
end

T = 2000.0
Tmeas=10.0
Nsave=10
fld_name="test"
L = 100


#for riproducibility
seed=22
ran_ng = Random.Xoshiro(seed)
M,s = build_model_state(L,rng=ran_ng)
p = ProgressMeasurer(T)
m = Measurer(M; times=Tmeas:Tmeas:T,name=fld_name, Nsave)
stats = (s,t)->(p(s,t); m(s,t))


@time run_RD!(s, M, T; stats = m, rng=ran_ng) 