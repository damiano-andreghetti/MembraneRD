using RD

function build_model_state(L)
    g,posx,posy = RD.gen_hex_lattice(L)
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
		nA[rand(1:N)]+=1
	end
	for i in 1:floor(Int,totB)
		nB[rand(1:N)]+=1
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
		nEA[rand(1:N)]+=1
	end
	for i in 1:EB_mem
		nEB[rand(1:N)]+=1
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

M,s = build_model_state(L)
p = ProgressMeasurer(T)
m = Measurer(M; times=0:Tmeas:T,name=fld_name, Nsave)
stats = (s,t)->(p(s,t); m(s,t))

@time run_RD!(s, M, T; stats = m) 