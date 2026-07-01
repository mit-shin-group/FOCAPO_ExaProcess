# FOCAPO_ExaProcess

Scripts to reproduce the flash benchmark results for the paper
"Harnessing GPU-Acceleration in Large-Scale Process Optimization".

The main benchmark runs JuMP/Ipopt, ExaModels/MadNLP on CPU, and
ExaModels/MadNLP on GPU for several scenario counts. The output is written to
`benchmark_table.csv`.

## Requirements

- Julia 1.12.1. The committed `Manifest.toml` was generated with this version.
- An NVIDIA GPU with a working CUDA driver.
- The HSL MA57 linear solver. The scripts use MA57 through `HSL_jll` for
  JuMP/Ipopt and through `MadNLPHSL` for MadNLP CPU solves.
- This repository, including the committed `data/inlet_data_*_scenarios.jls`
  files.

The full benchmark can take about 4 hours, depending on CPU, GPU, and solver
configuration.

## Step-by-step reproduction

1. Clone the repository and enter it.

   ```bash
   git clone <repository-url>
   cd FOCAPO_ExaProcess
   ```

2. Start Julia in the project environment.

   ```bash
   julia --project=.
   ```

3. Instantiate the exact package environment from `Manifest.toml`.

   ```julia
   using Pkg
   Pkg.instantiate()
   ```

4. Verify that CUDA is available to Julia.

   ```julia
   using CUDA
   CUDA.functional()
   ```

   This should print `true`. If it prints `false`, fix the NVIDIA driver/CUDA
   setup before running the GPU benchmark.

5. Run the full benchmark table.

   ```julia
   include("flash_benchmark_table.jl")
   ```

   The script runs the default scenario counts:

   ```julia
   [1, 10, 100, 1000, 5000]
   ```

   For each scenario count, it solves:

   - JuMP/Ipopt
   - ExaModels/MadNLP CPU
   - ExaModels/MadNLP GPU

6. Confirm the output.

   After the script finishes, it writes:

   ```text
   benchmark_table.csv
   ```

   The CSV contains one row per scenario count and these columns:

   ```text
   omega,jump_build_s,jump_solve_s,jump_iter,jump_status,exa_cpu_build_s,exa_cpu_solve_s,exa_cpu_iter,exa_cpu_status,exa_gpu_build_s,exa_gpu_solve_s,exa_gpu_iter,exa_gpu_status
   ```

   Successful ExaModels/MadNLP runs should report `SOLVE_SUCCEEDED`. JuMP/Ipopt
   runs may report `LOCALLY_SOLVED`; for some large cases, the status can depend
   on solver version, hardware, and runtime limits.

## Notes

- The benchmark script checks `CUDA.functional()` before running and stops if
  CUDA is not functional.
- `benchmark_table.csv` is overwritten each time `flash_benchmark_table.jl`
  completes.
- Timing values include machine-dependent build and solve times, so exact
  numbers are not expected to match across different systems.
