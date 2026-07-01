using ExaModels
using MadNLP
using MadNLPGPU
using MadNLPHSL
using Printf
using CUDA
using CUDSS
using Serialization
using NLPModelsIpopt

include("helper.jl")
include("constraint.jl")
include("Exa_unit_model.jl")


INLET_DATA_FILE = joinpath(@__DIR__, "data", "inlet_data_5000_scenarios.jls")
inlet_data = deserialize(INLET_DATA_FILE)

N_abs = 10
Q = 0.0
backend = CUDABackend()

if backend isa CUDABackend
    CUDA.functional() || error("CUDA is not functional; cannot build the GPU ExaModels core.")
end

built = build_absorber_flash_model(inlet_data, N_abs, backend, Q)
result = solve_exa_flash(built; backend=backend, verbose=true)
(; backend, status, objective, species_names, variables, diagnostics) = result

println("Number of scenarios: ", length(inlet_data.yV))
print_timing_summary("Build", built.build_timing)
print_timing_summary("Solve", result.solve_timing)