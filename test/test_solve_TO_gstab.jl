using StabTrussTopOpt
sto = StabTrussTopOpt
using Makie
using Test

# function example_relaxed_stab_opt()

file_dir = joinpath(pwd(), "test", "instances")
mesh_file = "tim.json"
load_file = "tim_load_support_case.json"

mesh_file = joinpath(file_dir, mesh_file)
load_file = joinpath(file_dir, load_file)

# parse ground mesh, supp, load from json
par = sto.parse_TO_problem(mesh_file, load_file)

# draw ground mesh
init_sc = Makie.Scene()
sto.draw_truss!(init_sc, par.X, par.T, par.Bound)
sto.draw_load!(init_sc, par.X, par.F)
init_sc

# define TO problem
# Steel-S235
E = 2.1e8 # kN/m^2
τ = 1.0
σ_c = 2.3e5 # kN/m^2
σ_t = 3.6e5
problem = sto.TO_define(par, E=E, τ=τ, σ_c=σ_c, σ_t=σ_t)

opt_a = sto.TO_solve_relaxed_stab(problem)

opt_a = abs.(opt_a)
# draw solution
final_sc = Scene()
sto.draw_truss!(final_sc, par.X, par.T, par.Bound, area = opt_a, line_width=1e4)
sto.draw_load!(final_sc, par.X, par.F, load_scale=1e-2)
display(final_sc)
# end

# example_relaxed_stab_opt()
