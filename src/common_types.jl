
mutable struct TOProblemPar
    ndim::Int          # dimension: 2 or 3
    n::Int             # num of nodes in ground structure
    m::Int             # num of elements in ground structure
    X::Matrix{Float64} # n x ndim matrix, node coords
    T::Matrix{Int}     # m x 2 matrix, element connectivity
    B::Matrix{Float64} # m x (ndim x 2) matrix, compatability matrix, i.e. gamma_i, the direction cosine, i ∈ [m]

    # TODO: multiple load cases, Array{Matrix{Float64}}
    F::Matrix{Float64} # n x ndim matrix, load

    # TODO: multiple support cases, Array{Matrix{Int}}
    Bound::Matrix{Int} # n x ndim matrix, boundary, = 1 fixed, = 0 free

    dof_map::Matrix{Int} # m x (2xndim) matrix, element-dof (e_id, 1..2xndim)
end

mutable struct TOProblem
    m::Int             # num of nodes in ground structure
    c::Vector{Float64} # m x 1 vector, coeff of liear objective function

    τ::Float64 # stability load factor, .>= 1 to ensure global stability
    E::Float64 # Young's  modulus
    σ_c::Float64 # compression stress limit
    σ_t::Float64 # tension stress limit

    # TODO: multi-load cases
    # τ::Vector{Float64}
    l::Vector{Float64} # m x 1 vector, element length
    # TODO: stored as a sparse vector
    f::Vector{Float64} # nfree x 1 vector, external force mapped in free dof indices

    Γ::SparseMatrixCSC{Float64, Int} # m x nfree sp compatability matrix
    
    Δ::Array{SparseMatrixCSC{Float64, Int}, 1} # m x sp(nfree x nfree) array, Δ_i for geometric stiffness matrix
    eK::Array{SparseMatrixCSC{Float64, Int}, 1} # m x sp(nfree x nfree), K_i for element stiffness matrix in full dof index form

    x_init::Vector{Float64} # init x value
end
