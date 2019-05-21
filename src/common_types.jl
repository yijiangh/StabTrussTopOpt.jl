
mutable struct TOProblemPar
    ndim::Int          # dimension: 2 or 3
    n::Int             # nodes in ground structure
    m::Int             # nodes in ground structure
    X::Matrix{Float64} # n x ndim matrix, node coords
    T::Matrix{Int}     # m x 2 matrix, element connectivity
    # E_DOF_map::Vector{Int} # m x (2xndim) matrix, element-dof (e_id, 1..2xndim)
    B::Matrix{Float64} # m x (ndim x 2) matrix, compatability matrix
    # TODO: multiple load cases, Array{Matrix{Float64}}
    F::Matrix{Float64} # n x ndim matrix, load
    # TODO: multiple support cases, Array{Matrix{Int}}
    Bound::Matrix{Int} # n x ndim matrix, boundary, = 1 fixed, = 0 free
end

# mutable struct TOProblem
#
# end
