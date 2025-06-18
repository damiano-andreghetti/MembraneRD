struct Model{G}
    g::G
    posx::Vector{Float64}
    posy::Vector{Float64}
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
	kAs::Float64
	kBs::Float64
    KMM::Float64
    rho_0::Float64
end

Base.length(M::Model) = nv(M.g)
