module JuMPFlashBenchmark

using ExaModels
using HSL_jll
using Ipopt
using JuMP
using Printf

include(joinpath(@__DIR__, "helper.jl"))
include(joinpath(@__DIR__, "constraint.jl"))
include(joinpath(@__DIR__, "JuMP_unit_model.jl"))

end

module ExaFlashBenchmark

using CUDA
using CUDSS
using ExaModels
using MadNLP
using MadNLPHSL
using MadNLPGPU
using NLPModelsIpopt
using Printf

include(joinpath(@__DIR__, "helper.jl"))
include(joinpath(@__DIR__, "constraint.jl"))
include(joinpath(@__DIR__, "Exa_unit_model.jl"))

end

using CUDA
using Serialization

const DEFAULT_SCENARIO_COUNTS = [1, 10, 100, 1000, 5000]
const DEFAULT_OUTPUT_FILE = joinpath(@__DIR__, "benchmark_table.csv")
const N_ABS = 10
const Q = 0.0

function scenario_file(omega)
    return joinpath(@__DIR__, "data", "inlet_data_$(omega)_scenarios.jls")
end

function load_inlet_data(omega)
    file = scenario_file(omega)
    isfile(file) || error("Missing inlet data file: $file")
    return deserialize(file)
end

function elapsed_without_compilation(timing)
    compile_time = hasproperty(timing, :compile_time) ? timing.compile_time : 0.0
    return max(timing.time - compile_time, 0.0)
end

function print_case_finished(omega, label, case)
    total_s = case.build_s + case.solve_s
    println(
        "Finished |Omega| = ", omega,
        " ", label,
        " in ", round(total_s; digits=6), " s",
        " (build ", round(case.build_s; digits=6), " s, solve ", round(case.solve_s; digits=6), " s).",
    )
end

function run_jump_case(inlet_data)
    built = JuMPFlashBenchmark.build_absorber_flash_model(inlet_data, N_ABS, Q)
    result = JuMPFlashBenchmark.solve_jump_flash(built; verbose=false)
    return (;
        build_s=elapsed_without_compilation(built.build_timing),
        solve_s=elapsed_without_compilation(result.solve_timing),
        iter=result.iterations,
        status=result.status,
    )
end

function run_exa_case(inlet_data, backend)
    built = ExaFlashBenchmark.build_absorber_flash_model(inlet_data, N_ABS, backend, Q)
    result = ExaFlashBenchmark.solve_exa_flash(built; backend=backend, verbose=false)
    return (;
        build_s=elapsed_without_compilation(built.build_timing),
        solve_s=elapsed_without_compilation(result.solve_timing),
        iter=result.iterations,
        status=result.status,
    )
end

function run_benchmark_row(omega)
    println("Running |Omega| = ", omega)
    inlet_data = load_inlet_data(omega)
    actual_omega = length(inlet_data.yV)
    actual_omega == omega || error("Expected $omega scenarios in $(scenario_file(omega)); found $actual_omega.")

    jump = run_jump_case(inlet_data)
    print_case_finished(omega, "JuMP/Ipopt", jump)

    exa_cpu = run_exa_case(inlet_data, nothing)
    print_case_finished(omega, "ExaModels/MadNLP CPU", exa_cpu)

    exa_gpu = run_exa_case(inlet_data, CUDABackend())
    print_case_finished(omega, "ExaModels/MadNLP GPU", exa_gpu)

    return (;
        omega=omega,
        jump_build_s=jump.build_s,
        jump_solve_s=jump.solve_s,
        jump_iter=jump.iter,
        jump_status=jump.status,
        exa_cpu_build_s=exa_cpu.build_s,
        exa_cpu_solve_s=exa_cpu.solve_s,
        exa_cpu_iter=exa_cpu.iter,
        exa_cpu_status=exa_cpu.status,
        exa_gpu_build_s=exa_gpu.build_s,
        exa_gpu_solve_s=exa_gpu.solve_s,
        exa_gpu_iter=exa_gpu.iter,
        exa_gpu_status=exa_gpu.status,
    )
end

function csv_escape(value)
    text = string(value)
    if occursin('"', text) || occursin(',', text) || occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function write_csv(rows, output_file)
    columns = (
        :omega,
        :jump_build_s,
        :jump_solve_s,
        :jump_iter,
        :jump_status,
        :exa_cpu_build_s,
        :exa_cpu_solve_s,
        :exa_cpu_iter,
        :exa_cpu_status,
        :exa_gpu_build_s,
        :exa_gpu_solve_s,
        :exa_gpu_iter,
        :exa_gpu_status,
    )

    open(output_file, "w") do io
        println(io, join(String.(columns), ","))
        for row in rows
            println(io, join((csv_escape(getproperty(row, column)) for column in columns), ","))
        end
    end
end

CUDA.functional() || error("CUDA is not functional; cannot run ExaModels/MadNLP GPU benchmark.")

rows = [run_benchmark_row(omega) for omega in DEFAULT_SCENARIO_COUNTS]
write_csv(rows, DEFAULT_OUTPUT_FILE)
println("Wrote benchmark table to ", DEFAULT_OUTPUT_FILE)


