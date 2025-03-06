using ExponentialQueues, JLD, Random

struct Model
    neig::Matrix{Int}
    centers::Matrix{Float64}
    dA::Float64
    dB::Float64
    dEA::Float64
    dEB::Float64
    kAc::Float64
    kBc::Float64
    kAa::Float64
    kAd::Float64
    kBa::Float64
    kBd::Float64
    KMM::Float64
    rho_0::Float64
end

struct State{T}
	nA::Vector{T}
	nB::Vector{T}
	nEA::Vector{T}
	nEB::Vector{T}
	cytoEA::Base.RefValue{T}
	cytoEB::Base.RefValue{T}
end



Base.length(s::State) = length(s.nEA)
Base.length(M::Model) = size(M.neig, 1)

#make measures here
function Measurer(M::Model; times)
	measures = zeros(length(times),5)
    Nc = size(M.neig,1)
	phi_field = zeros(Nc)
	phi_av = 0.0
	binder_cumulant = 0.0
	next = 1
    function measure(t, s)
        while t > times[next]  
            totA, totB = sum(s.nA), sum(s.nB)
            #calculate phi at each cell (normalized, goes from -1, to 1)
            for i in eachindex(phi_field)
                #in some denominators I added + 10^-10 in order to avoid NaNs
                phi_field[i] = (s.nB[i] - s.nA[i])/(s.nA[i] + s.nB[i] + 1e-10)
            end
            phi_av = (totB - totA) / (totA + totB)
            println("T = $t and <ϕ>/c = $phi_av")
            #measure phi
            measures[next,1] = phi_av
            #measure Binder cumulant
            binder_cumulant = 1-
                ((sum(x->(x-phi_av)^4, phi_field)/Nc) /
                ((sum(x->(x-phi_av)^2, phi_field)/Nc)^2 + 10^-10))/3
            measures[next,2] = binder_cumulant
            #measure enzymes in cytosolic reservoir
            measures[next,3] = s.cytoEA
            measures[next,4] = s.cytoEB
            println("Binder cumulant: ", measures[next,2])			
            #the following is to check for convergence of rho to 1
            rho = M.rho_0 * ((M.kAd / M.kAa) + totA) / ((M.kBd / M.kBa) + totB)
            measures[next,5] = rho
            next += 1
        end
    end
end

function run_RD!(s::State, M::Model, T; 
        stats = (x...)->nothing, 
        rng = Random.default_rng())

    n = length(M)
    QA,QB,QEA,QEB,QcatA,QcatB = (ExponentialQueue(n) for _ in 1:6)
    QattEA = QA * 1.0
    QattEB = QB * 1.0

    function update(i)
        QA[i] = s.nA[i]
        QB[i] = s.nB[i]
        QEA[i] = s.nEA[i]
        QEB[i] = s.nEB[i]
        QcatA[i] = s.nEA[i] * s.nB[i] / (s.nB[i] + M.KMM)
        QcatB[i] = s.nEB[i] * s.nA[i] / (s.nA[i] + M.KMM)
        QattEA.f[] = s.cytoEA[] * M.kAa
        QattEB.f[] = s.cytoEB[] * M.kBa
    end

    foreach(update, 1:n)

    #arrival is chosen uniformly between its neighbours
    rand_neighbor(i) = rand(rng, @view M.neig[i,:])

    Q = NestedQueue(
            :difA => QA * M.dA,
            :difB => QB * M.dB,
            :catA => QcatA * M.kAc,
            :catB => QcatB * M.kBc,
            :attEA => QattEA,
            :attEB => QattEB,
            :difEA => QEA * M.dEA,
            :difEB => QEB * M.dEB
        )

    println("starting simulation, $(length(Q)) events in the queue")

    t::Float64 = 0.0
    while !isempty(Q)
        (ev, i), dt = peek(Q; rng)
        t += dt
        t > T && break
        stats(t, s)
        if ev === :difA #diffusion of specie A
            j = rand_neighbor(i)
            s.nA[i] -= 1
            s.nA[j] += 1
            update(i)
            update(j)
        elseif ev === :difB #diffusion of specie B
            j = rand_neighbor(i)
            s.nB[i] -= 1
            s.nB[j] += 1
            update(i)
            update(j)
        elseif ev === :catA #B+EA->A+EA
            s.nB[i] -= 1
            s.nA[i] += 1
            update(i)
        elseif ev === :catB #A+EB->B+EB
            s.nA[i] -= 1
            s.nB[i] += 1
            update(i)
        elseif ev === :attEA #attachment of EA from cytosol
            s.cytoEA[] -= 1
            s.nEA[i] += 1
            update(i)
        elseif ev === :detEA #detachment of EA
            s.nEA[i] -= 1
            s.cytoEA[] += 1
            update(i)
        elseif ev === :attEB #attachment of EB from cytosol
            s.cytoEB[] -= 1
            s.nEB[i] += 1
            update(i)
        elseif ev === :detEB #detachment of EB
            s.nEB[i] -= 1
            s.cytoEB[] += 1
            update(i)
        elseif ev === :difEA #diffusion of EA
            j = rand_neighbor(i)
            s.nEA[i] -= 1
            s.nEA[j] += 1
            update(i)
            update(j)
        elseif ev === :difEB #diffusion of EB
            j = rand_neighbor(i)
            s.nEB[i] -= 1
            s.nEB[j] += 1
            update(i)
            update(j)
        end
    end
end
