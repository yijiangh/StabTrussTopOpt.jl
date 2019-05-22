
# __precompile__()

module StabTrussTopOpt

# dependencies
using JSON
using Makie
using LinearAlgebra
using SparseArrays

# export TOProblemPar

include("common_types.jl")
include("draw_utils.jl")
include("parse_json.jl")
include("tto_define.jl")
include("solve.jl")

end # module
