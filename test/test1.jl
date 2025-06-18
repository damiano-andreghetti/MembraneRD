using MembraneRD, Test, StableRNGs

function build_model_state(L; rng = Random.default_rng())
    g,posx,posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    unit_length = 1.0
    theta = 0.2126944621086619#Stot/V
    Ac = unit_length^2#area associated to each cell, set unit lengthscale
    Stot = N*Ac #area of each cell is assumed to be 1, setting unit lengthscale
    V = Stot/theta
    d_timescale = 0.01#this sets unit timescale
    dA,dB,dEA,dEB = d_timescale, d_timescale, d_timescale, d_timescale
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
	kAs=0.0
	kBs=0.0

    M = Model(g, posx, posy, dA, dB, dEA, dEB, kAc, kBc, kAa, kAd, kBa, kBd, kAs,kBs,KMM, 0.0)

    totmol = N * 10
    totA, totB = floor(Int, totmol / 2), floor(Int, totmol / 2)
    totEA, totEB = floor(Int, 0.1*N), floor(Int, 0.1*N)
    memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+((theta/(kAd_th/kAa_th))*(totA/Stot))))
    memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+((theta/(kBd_th/kBa_th))*(totB/Stot))))
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    s = State(M, totA, totB, memEA, memEB, cytoEA, cytoEB; rng)
    M,s
end

@testset "reproducibility" begin
    T = 2000.0
    L = 100
    rng = StableRNG(22)
    M,s = build_model_state(L; rng)
    run_RD!(s, M, T; rng)
    @test sum(s.nA) == 37442
    @test sum(s.nB) == 62558
end

