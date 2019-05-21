using StabTrussTopOpt
sto = StabTrussTopOpt
using Makie

file_dir = joinpath(pwd(), "test", "instances")
mesh_file = "tim.json"
load_file = "tim_load_support_case.json"

mesh_file = joinpath(file_dir, mesh_file)
load_file = joinpath(file_dir, load_file)

par = sto.parse_TO_problem(mesh_file, load_file)

init_sc = Makie.Scene()
sto.draw_truss!(init_sc, par.X, par.T, par.Bound)
sto.draw_load!(init_sc, par.X, par.F)
init_sc

# to_problem = TO_define(to_problem_par)

# solve_TO(to_problem)

# draw(problem)
