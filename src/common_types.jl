
mutable struct TOProblemPar
    ndim::Int          # dimension: 2 or 3
    n::Int             # nodes in ground structure
    m::Int             # nodes in ground structure
    X::Matrix{Float64} # n x ndim matrix, node coords
    T::Matrix{Int}     # m x 2 matrix, element connectivity
    B::Matrix{Float64} # m x (ndim x 2) matrix, compatability matrix, i.e. gamma_i, the direction cosine, i ∈ [m]

    # TODO: multiple load cases, Array{Matrix{Float64}}
    F::Matrix{Float64} # n x ndim matrix, load

    # TODO: multiple support cases, Array{Matrix{Int}}
    Bound::Matrix{Int} # n x ndim matrix, boundary, = 1 fixed, = 0 free

    dof_map::Vector{Int} # m x (2xndim) matrix, element-dof (e_id, 1..2xndim)
end

mutable struct TOProblem
    n_primal_var::Int # number of primal variables
    n_matrix_ineq:Int # number of matrix ineq
    n_linear_ineq:Int # number of linear ineq

    c::Vector{Float64} # n_primal_var x 1 vector, coeff of liear objective function

    x_init::Vector{Float64} # init x value

    τ::Vector{Float64} # stability load factor, .>= 1 to ensure global stability
    l::Vector{Float64} # m x 1 vector, element length
    Delta::Array{SparseMatrixCSC{Float64, Int},1} # m x sp(nfree x nfree), Δ_i for geometric stiffness matrix in full dof index form
    eK::Array{SparseMatrixCSC{Float64, Int},1} # m x sp(nfree x nfree), K_i for element stiffness matrix in full dof index form
end
