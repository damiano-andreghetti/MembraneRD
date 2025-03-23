struct State{T}
	nA::Vector{T}
	nB::Vector{T}
	nEA::Vector{T}
	nEB::Vector{T}
	cytoEA::Base.RefValue{T}
	cytoEB::Base.RefValue{T}
end

function State(M::Model, totA::Int, totB::Int, memEA::Int, memEB::Int, cytoEA::Int, cytoEB::Int; rng = Random.default_rng())
    N = length(M)
    nA, nB, nEA, nEB = fill(0,N), fill(0,N), fill(0,N), fill(0,N)

	for _ in 1:totA
		nA[rand(rng, 1:N)]+=1
	end
	for _ in 1:totB
		nB[rand(rng, 1:N)]+=1
	end
    for _ in 1:memEA
		nEA[rand(rng, 1:N)]+=1
	end
	for _ in 1:memEB
		nEB[rand(rng, 1:N)]+=1
	end
    State(nA, nB, nEA, nEB, Ref(cytoEA), Ref(cytoEB))
end

Base.length(s::State) = length(s.nEA)