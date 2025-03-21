module MembraneRD

using ExponentialQueues, Random

export Model, run_RD!, State, gen_hex_lattice, gen_rect_lattice

include("lattice.jl")
include("model.jl")
include("gillespie.jl")

end
