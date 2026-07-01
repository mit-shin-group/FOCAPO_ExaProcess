using ExaModels
using Serialization
using JuMP, Revise, Ipopt
using Printf
using HSL_jll

include("helper.jl")
include("constraint.jl")
include("JuMP_unit_model.jl")

INLET_DATA_FILE = joinpath(@__DIR__, "data", "inlet_data_100_scenarios.jls")
inlet_data = deserialize(INLET_DATA_FILE)

N_abs = 10
Q = 0.0

built = build_absorber_flash_model(inlet_data, N_abs, Q)
result = solve_jump_flash(built; verbose=true)
(; backend, status, objective, species_names, variables, diagnostics) = result

println("Number of scenarios: ", length(inlet_data.yV))
print_timing_summary("Build", built.build_timing)
print_timing_summary("Solve", result.solve_timing)

