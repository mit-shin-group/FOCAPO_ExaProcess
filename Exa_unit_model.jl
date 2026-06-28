function _build_absorber_flash_model(inlet_data, N_abs = 1, backend=nothing, Q=0.0)
    helper_self_check()
    to_backend_array(x) = ExaModels.convert_array(x, backend)

    nS = length(inlet_data.yV) # number of scenarios
    nC = length(species_index) # number of species

    @assert nC == NSPECIES
    @assert NREACTIONS == 5
    @assert nS > 0 "inlet_data must contain at least one scenario."
    @assert length(inlet_data.xL) == nS "xL must have one entry per scenario."
    @assert length(inlet_data.V_in_s) == nS "V_in_s must have one entry per scenario."
    @assert length(inlet_data.L_in_s) == nS "L_in_s must have one entry per scenario."
    @assert length(inlet_data.T_L_in_s) == nS "T_L_in_s must have one entry per scenario."
    @assert length(inlet_data.T_V_in_s) == nS "T_V_in_s must have one entry per scenario."
    @assert all(length(y) == nC for y in inlet_data.yV) "Each yV scenario must contain nC species."
    @assert all(length(x) == nC for x in inlet_data.xL) "Each xL scenario must contain nC species."

    y_in_host = Float64.(permutedims(reduce(hcat, inlet_data.yV)))
    x_in_host = Float64.(permutedims(reduce(hcat, inlet_data.xL)))
    V_in_host = Float64.(inlet_data.V_in_s)
    L_in_host = Float64.(inlet_data.L_in_s)
    T_L_in_host = Float64.(inlet_data.T_L_in_s)
    T_V_in_host = Float64.(inlet_data.T_V_in_s)
    Q_host = if Q isa Number
        fill(Float64(Q), nS, N_abs)
    else
        Q_array = Float64.(Q)
        if size(Q_array) == (nS,)
            repeat(reshape(Q_array, nS, 1), 1, N_abs)
        else
            Q_array
        end
    end

    @assert size(y_in_host) == (nS, nC)
    @assert size(x_in_host) == (nS, nC)
    @assert all(length(v) == nS for v in (V_in_host, L_in_host, T_L_in_host, T_V_in_host))
    @assert size(Q_host) == (nS, N_abs)

    P_abs = 108.02 / 101.325
    density = 1050.0
    d_min = 0.1 * 0.64  # [m], absorber diameter lower bound
    d_max = 10.0 * 0.64  # [m], absorber diameter upper bound
    u_max_G = 10.0 * 0.3048  # [m/s], maximum gas superficial velocity (10 ft/s)
    abs_d_dummy_host = fill(d_min, nS, N_abs)

    species_order = sort(collect(keys(species_index)); by=spec -> species_index[spec])
    MW_VEC_host = map(spec -> MW[spec], species_order)
    @assert length(MW_VEC_host) == nC

    y_in = to_backend_array(y_in_host)
    x_in = to_backend_array(x_in_host)
    V_in = to_backend_array(V_in_host)
    L_in = to_backend_array(L_in_host)
    T_L_in = to_backend_array(T_L_in_host)
    T_V_in = to_backend_array(T_V_in_host)
    Q_in = to_backend_array(Q_host)
    abs_d_dummy_values = to_backend_array(abs_d_dummy_host)
    MW_VEC = to_backend_array(MW_VEC_host)
    stoich_backend = to_backend_array(stoich)

    stage_array(v) = to_backend_array(repeat(v, 1, N_abs))
    stage_species(v) =
        to_backend_array(repeat(reshape(v, 1, 1, length(v)), nS, N_abs, 1))
    stage_scenario_species(A) =
        to_backend_array(repeat(reshape(A, size(A, 1), 1, size(A, 2)), 1, N_abs, 1))

    species_ids = collect(1:nC)
    species_row = reshape(species_ids, 1, nC)

    abs_L_start = stage_array(L_in_host)
    abs_V_start = stage_array(V_in_host)
    abs_T_start = stage_array(T_V_in_host)
    abs_log_x_lvar = stage_species(liquid_log_lower.(species_ids))
    abs_log_x_uvar = stage_species(liquid_log_upper.(species_ids))
    abs_log_y_lvar = stage_species(vapor_log_lower.(species_ids))
    abs_log_y_uvar = stage_species(vapor_log_upper.(species_ids))
    abs_log_x_start = stage_scenario_species(liquid_log_start.(species_row, x_in_host))
    abs_log_y_start = stage_scenario_species(vapor_log_start.(species_row, y_in_host))

    @assert size(abs_L_start) == (nS, N_abs)
    @assert size(abs_V_start) == (nS, N_abs)
    @assert size(abs_T_start) == (nS, N_abs)
    @assert size(abs_log_x_lvar) == (nS, N_abs, nC)
    @assert size(abs_log_x_uvar) == (nS, N_abs, nC)
    @assert size(abs_log_y_lvar) == (nS, N_abs, nC)
    @assert size(abs_log_y_uvar) == (nS, N_abs, nC)
    @assert size(abs_log_x_start) == (nS, N_abs, nC)
    @assert size(abs_log_y_start) == (nS, N_abs, nC)

    c = ExaCore(; backend=backend, concrete=Val(true))

    @add_par(c, abs_y_in, y_in)
    @add_par(c, abs_x_in, x_in)
    @add_par(c, abs_V_in, V_in)
    @add_par(c, abs_L_in, L_in)
    @add_par(c, abs_T_L_in, T_L_in)
    @add_par(c, abs_T_V_in, T_V_in)
    @add_par(c, abs_Q, Q_in)
    @add_par(c, abs_d_dummy, abs_d_dummy_values)
    @add_par(c, abs_stoich, stoich_backend)
    @add_par(c, abs_MW, MW_VEC)

    @add_var(c, abs_L, nS, N_abs;
    lvar = 0.0,
    start = abs_L_start,
    )

    @add_var(c, abs_V, nS, N_abs;
    lvar = 0.0,
    start = abs_V_start,
    )

    @add_var(c, abs_T, nS, N_abs;
    lvar = 273.15,
    uvar = 400.0,
    start = abs_T_start,
    )

    @add_var(c, abs_log_x, nS, N_abs, nC;
    lvar = abs_log_x_lvar,
    uvar = abs_log_x_uvar,
    start = abs_log_x_start,
    )

    @add_var(c, abs_log_y, nS, N_abs, nC;
    lvar = abs_log_y_lvar,
    uvar = abs_log_y_uvar,
    start = abs_log_y_start,
    )

    @add_var(c, abs_xi, nS, N_abs, NREACTIONS;
    start = 0.0,
    )

    @add_var(c, abs_d, 1;
    lvar = 0,
    uvar = Inf,
    start = 0.64,
    )

    @add_expr(c, abs_x, exp(abs_log_x[s, n, k]) for s in 1:nS, n in 1:N_abs, k in 1:nC)

    @add_expr(c, abs_y, exp(abs_log_y[s, n, k]) for s in 1:nS, n in 1:N_abs, k in 1:nC)

    @add_expr(c, abs_x_Li, abs_L_in[s] * abs_x_in[s, k] for s in 1:nS, k in 1:nC)

    @add_expr(c, abs_y_Vi, abs_V_in[s] * abs_y_in[s, k] for s in 1:nS, k in 1:nC)

    @add_obj(c, abs_d[1])

    REACTIVE_IDX = (MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m, CO2, H2O)

    cons_abs_MB = nothing
    cons_abs_MB_1 = nothing
    cons_abs_MB_N = nothing
    cons_abs_MB_n = nothing

    if N_abs == 1
        @add_con(c, cons_abs_MB,
            cons_reactive_mb_single(abs_x, abs_y, abs_L, abs_V, abs_y_Vi, abs_x_Li, abs_stoich, abs_xi, s, k)
            for s in 1:nS, k in REACTIVE_IDX
        )
    else
        @add_con(c, cons_abs_MB_1,
            cons_reactive_mb_first(abs_x, abs_y, abs_L, abs_V, abs_y_Vi, abs_stoich, abs_xi, s, k)
            for s in 1:nS, k in REACTIVE_IDX
        )

        @add_con(c, cons_abs_MB_N,
            cons_reactive_mb_last(abs_x, abs_y, abs_L, abs_V, abs_x_Li, abs_stoich, abs_xi, s, k, N_abs)
            for s in 1:nS, k in REACTIVE_IDX
        )

        if N_abs > 2
            @add_con(c, cons_abs_MB_n,
                cons_reactive_mb_middle(abs_x, abs_y, abs_L, abs_V, abs_stoich, abs_xi, s, n, k)
            for s in 1:nS, n in 2:(N_abs - 1), k in REACTIVE_IDX
            )
        end
    end

    @add_con(c, cons_abs_inert_MB, cons_inert_mb(abs_y, abs_V, abs_y_Vi, s, n, k) for s in 1:nS, n in 1:N_abs, k in (N2, O2))

    @add_con(c, cons_abs_eq1, cons_rxne1(abs_T, abs_log_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_eq2, cons_rxne2(abs_T, abs_log_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_eq3, cons_rxne3(abs_T, abs_log_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_eq4, cons_rxne4(abs_T, abs_log_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_eq5, cons_rxne5(abs_T, abs_log_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_expr(c, abs_MEA_app, compute_stage_apparent_liquid_MEA(abs_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_expr(c, abs_H2O_app, compute_stage_apparent_liquid_H2O(abs_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_expr(c, abs_CO2_app, compute_stage_apparent_liquid_CO2(abs_x, s, n) for s in 1:nS, n in 1:N_abs)

    @add_expr(c, He_CO2_expr, get_H(abs_MEA_app[s, n], abs_H2O_app[s, n], abs_T[s, n]) / 101325.0 for s in 1:nS, n in 1:N_abs)

    @add_expr(c, P_H2O_sat_expr, get_Pvap(abs_T[s, n]) / 101325.0 for s in 1:nS, n in 1:N_abs)

    @add_expr(c, MW_avg_expr, sum(abs_x[s, n, k] * abs_MW[k] for k in 1:nC) * 1e-3 for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_H2O_VLE, cons_h2o_vle(P_abs, abs_y, abs_x, P_H2O_sat_expr, s, n) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_CO2_VLE, cons_co2_vle(P_abs, abs_y, abs_x, He_CO2_expr, MW_avg_expr, density, s, n) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_sum_x, mole_fraction_sum(abs_x, s, n, nC) for s in 1:nS, n in 1:N_abs)

    @add_con(c, cons_abs_sum_y, mole_fraction_sum(abs_y, s, n, nC) for s in 1:nS, n in 1:N_abs)

    # Unit check: abs_V [mol/s] * R [Pa*m^3/(mol*K)] * abs_T [K] / P [Pa] = [m^3/s].
    # Unit check: u_max_G [m/s] * pi * abs_d^2 / 4 [m^2] = [m^3/s].
    @add_con(c, cons_abs_superficial_velocity_limit,
        cons_superficial_velocity_limit(abs_V, abs_T, abs_d, P_abs, u_max_G, s, n)
        for s in 1:nS, n in 1:N_abs;
        lcon = -Inf,
        ucon = 0.0,
    )

    @add_expr(c, abs_app_Li_MEA, compute_inlet_apparent_liquid_MEA(abs_x_Li, s) for s in 1:nS)

    @add_expr(c, abs_app_Li_H2O, compute_inlet_apparent_liquid_H2O(abs_x_Li, s) for s in 1:nS)

    @add_expr(c, abs_app_Li_CO2, compute_inlet_apparent_liquid_CO2(abs_x_Li, s) for s in 1:nS)

    @add_expr(c, abs_H_liq_in, compute_inlet_liquid_enthalpy((abs_app_Li_MEA, abs_app_Li_H2O, abs_app_Li_CO2), abs_T_L_in, s) for s in 1:nS)

    @add_expr(c, abs_H_liq_out, compute_stage_liquid_enthalpy((abs_MEA_app, abs_H2O_app, abs_CO2_app), abs_T, abs_L, s, n) for s in 1:nS, n in 1:N_abs)

    @add_expr(c, abs_H_vap_in, compute_inlet_vapor_enthalpy(abs_y_Vi, abs_T_V_in, s) for s in 1:nS)

    @add_expr(c, abs_H_vap_out, compute_stage_vapor_enthalpy(abs_y, abs_T, abs_V, s, n) for s in 1:nS, n in 1:N_abs)

    cons_abs_energy = nothing
    cons_abs_energy_1 = nothing
    cons_abs_energy_N = nothing
    cons_abs_energy_n = nothing

    if N_abs == 1
        @add_con(c, cons_abs_energy,
            cons_energy_single(abs_H_liq_in, abs_H_vap_in, abs_Q, abs_H_liq_out, abs_H_vap_out, s)
            for s in 1:nS
        )
    else
        @add_con(c, cons_abs_energy_1, cons_energy_first(abs_H_liq_out, abs_H_vap_in, abs_Q, abs_H_vap_out, s) for s in 1:nS)

        @add_con(c, cons_abs_energy_N,
            cons_energy_last(abs_H_liq_in, abs_H_vap_out, abs_Q, abs_H_liq_out, s, N_abs)
            for s in 1:nS
        )

        if N_abs > 2
            @add_con(c, cons_abs_energy_n,
                cons_energy_middle(abs_H_liq_out, abs_H_vap_out, abs_Q, s, n)
            for s in 1:nS, n in 2:(N_abs - 1)
            )
        end
    end

    m = ExaModel(c)

    variables = (;
        abs_L,
        abs_V,
        abs_T,
        abs_log_x,
        abs_log_y,
        abs_xi,
        abs_d,
    )
    constraints = (;
        cons_abs_MB,
        cons_abs_MB_1,
        cons_abs_MB_N,
        cons_abs_MB_n,
        cons_abs_inert_MB,
        cons_abs_eq1,
        cons_abs_eq2,
        cons_abs_eq3,
        cons_abs_eq4,
        cons_abs_eq5,
        cons_abs_H2O_VLE,
        cons_abs_CO2_VLE,
        cons_abs_sum_x,
        cons_abs_sum_y,
        cons_abs_superficial_velocity_limit,
        cons_abs_energy,
        cons_abs_energy_1,
        cons_abs_energy_N,
        cons_abs_energy_n,
    )

    return (model=m, variables=variables, constraints=constraints)
end

function build_absorber_flash_model(inlet_data, N_abs = 1, backend=nothing, Q=0.0)
    build_timing = @timed _build_absorber_flash_model(inlet_data, N_abs, backend, Q)
    built = build_timing.value
    return merge(built, (; build_timing=build_timing))
end

function solve_exa_flash(built; backend=nothing, verbose=false)
    model = built.model
    variables = built.variables
    abs_L = variables.abs_L
    abs_V = variables.abs_V
    abs_T = variables.abs_T
    abs_log_x = variables.abs_log_x
    abs_log_y = variables.abs_log_y
    abs_xi = variables.abs_xi
    abs_d = variables.abs_d
    P_abs = 108.02 / 101.325
    u_max_G = 10.0 * 0.3048
    nS, N_abs = abs_L.size
    use_cuda_backend = backend isa CUDABackend

    if verbose
        println("Successfully added absorber flash constraints.")
        println("Reactive material balance constraints: ", length(species_react) * nS * N_abs)
        println("Inert vapor balance constraints: ", 2 * nS * N_abs)
        println("Reaction equilibrium constraints: ", 5 * nS * N_abs)
        println("VLE constraints: ", 2 * nS * N_abs)
        println("Mole fraction summation constraints: ", 2 * nS * N_abs)
        println("Superficial velocity constraints: ", nS * N_abs)
        println("Energy balance constraints: ", nS * N_abs)
        println(model)
    end

    print_level = verbose ? MadNLP.INFO : MadNLP.ERROR
    solve_timing = @timed begin
        if use_cuda_backend
            if verbose
                println("Using CUDABackend with CUDSS linear solver.")
            end
            result = madnlp(model; linear_solver=MadNLPGPU.CUDSSSolver, print_level=print_level, tol=1e-8)
            CUDA.synchronize()
            result
        else
            if verbose
                println("Using backend: ", backend, " with MadNLP on CPU.")
            end
            madnlp(model; linear_solver=MadNLPHSL.Ma57Solver, print_level=print_level)
        end
    end
    result = solve_timing.value

    abs_d_val = Array(solution(result, abs_d))
    abs_L_val = Array(solution(result, abs_L))
    abs_V_val = Array(solution(result, abs_V))
    abs_T_val = Array(solution(result, abs_T))
    abs_log_x_val = Array(solution(result, abs_log_x))
    abs_log_y_val = Array(solution(result, abs_log_y))
    abs_xi_val = Array(solution(result, abs_xi))
    vol_flow = abs_V_val .* R .* abs_T_val ./ (P_abs * 101325.0)
    residual = vol_flow .- u_max_G * pi * only(abs_d_val)^2 / 4.0
    d_required = sqrt(4.0 * maximum(vol_flow) / (u_max_G * pi))
    d_actual = only(abs_d_val)
    objective = result.objective
    status = result.status
    species_names = String.(first.(sort(collect(species_index), by=last)))

    if verbose
        println("status = ", status)
        println("backend:", backend)
        println("MADNLP objective:", objective)
        println("max signed g: ", maximum(residual))
        println("min signed g: ", minimum(residual))
        println("d required by max flow: ", d_required)
        println("d actual: ", d_actual)
        println("slack: ", d_actual - d_required)
    end

    return (;
        backend="ExaModels/MadNLP",
        status=string(status),
        objective=Float64(objective),
        species_names=species_names,
        variables=(;
            abs_d=abs_d_val,
            abs_L=abs_L_val,
            abs_V=abs_V_val,
            abs_T=abs_T_val,
            abs_log_x=abs_log_x_val,
            abs_log_y=abs_log_y_val,
            abs_xi=abs_xi_val,
        ),
        diagnostics=(;
            superficial_velocity_residual=residual,
            d_required=d_required,
            d_actual=d_actual,
            diameter_slack=d_actual - d_required,
        ),
        solve_timing=solve_timing,
    )
end