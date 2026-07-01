function _build_absorber_flash_model(inlet_data, N_abs=1, Q=0.0)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    set_optimizer_attribute(model, "max_iter", 1000)
    set_optimizer_attribute(model, "linear_solver", "ma57")


    #region absorber variable declaration
    S = 1:length(inlet_data.yV)
    nC = length(species_index)
    abs_Q = if Q isa Number
        fill(Float64(Q), length(S), N_abs)
    else
        Q_array = Float64.(Q)
        if size(Q_array) == (length(S),)
            repeat(reshape(Q_array, length(S), 1), 1, N_abs)
        else
            Q_array
        end
    end

    @variable(model, abs_y_in[s in S, i=1:nC] == inlet_data.yV[s][i])
    @variable(model, abs_V_in[s in S] == inlet_data.V_in_s[s]) # Vapor feed flowrate [=] mol/s
    @variable(model, abs_T_L_in[s in S] == inlet_data.T_L_in_s[s]) # Liquid feed temperature [=] K
    @variable(model, abs_T_V_in[s in S] == inlet_data.T_V_in_s[s]) # Vapor feed temperature [=] K
    @variable(model, abs_L_in[s in S] == inlet_data.L_in_s[s]) # Liquid feed flowrate [=] mol/s
    @variable(model, abs_x_in[s in S, c=1:nC] == inlet_data.xL[s][c])

    # absorber states
    P_abs = 108.02 / 101.325 # [=] atm, pressure of the absorber
    density = 1050.0 # assumed constant density of MEA solution [=] kg/m3
    u_max_G = 10.0 * 0.3048  # [m/s], maximum gas superficial velocity (10 ft/s)
    @variables(model, begin
        abs_L[s in S, n=1:N_abs] >= 0   # Liquid flowrate
        abs_V[s in S, n=1:N_abs] >= 0   # Vapor flowrate
        273.15 <= abs_T[s in S, 1:N_abs] <= 400 # Temperature [=] K, temperature limited to avoid thermal degradation (Davis, J., Rochelle, G., 2009. Thermal degradation of monoethanolamine at stripper conditions. Energy Procedia 1, 327–333.)

        abs_xi[s in S, n=1:N_abs, r=1:5]         # Extent of reaction for reaction 1-5
    end)
    @variable(model, liquid_log_lower(c) <= abs_log_x[s in S, n=1:N_abs, c=1:nC] <= liquid_log_upper(c))
    @variable(model, vapor_log_lower(c) <= abs_log_y[s in S, n=1:N_abs, c=1:nC] <= vapor_log_upper(c))
    @variable(model, abs_d[1:1] >= 0, start=0.64)

    abs_x = exp.(abs_log_x)
    abs_y = exp.(abs_log_y)

    #endregion

    @objective(model, Min, abs_d[1])

    #region ------------absorber constraints----------------
    # Initial amount of the 11 species [=] moles
    @expression(model, abs_x_Li[s in S, c=1:nC], abs_L_in[s] * abs_x_in[s, c])
    @expression(model, abs_y_Vi[s in S, c=1:nC], abs_V_in[s] * abs_y_in[s, c])

    # Apparent species amount in liquid feed
    @expression(model, abs_app_Li_MEA[s in S], compute_inlet_apparent_liquid_MEA(abs_x_Li, s))
    @expression(model, abs_app_Li_H2O[s in S], compute_inlet_apparent_liquid_H2O(abs_x_Li, s))
    @expression(model, abs_app_Li_CO2[s in S], compute_inlet_apparent_liquid_CO2(abs_x_Li, s))
    abs_app_Li = (abs_app_Li_MEA, abs_app_Li_H2O, abs_app_Li_CO2)

    cons_abs_MB = nothing
    cons_abs_MB_1 = nothing
    cons_abs_MB_N = nothing
    cons_abs_MB_n = nothing

    # Mole balance constraints
    REACTIVE_IDX = (MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m, CO2, H2O)
    if N_abs == 1
        cons_abs_MB = @constraint(model, [s in S, k in REACTIVE_IDX],
            cons_reactive_mb_single(abs_x, abs_y, abs_L, abs_V, abs_y_Vi, abs_x_Li, stoich, abs_xi, s, k) == 0)
    else
        cons_abs_MB_1 = @constraint(model, [s in S, k in REACTIVE_IDX],
            cons_reactive_mb_first(abs_x, abs_y, abs_L, abs_V, abs_y_Vi, stoich, abs_xi, s, k) == 0)

        cons_abs_MB_N = @constraint(model, [s in S, k in REACTIVE_IDX],
            cons_reactive_mb_last(abs_x, abs_y, abs_L, abs_V, abs_x_Li, stoich, abs_xi, s, k, N_abs) == 0)

        if N_abs > 2
        cons_abs_MB_n = @constraint(model, [s in S, n in 2:(N_abs - 1), k in REACTIVE_IDX],
                cons_reactive_mb_middle(abs_x, abs_y, abs_L, abs_V, stoich, abs_xi, s, n, k) == 0)
        end
    end

    # N2 and O2 mole balance constraints (inert species)
    cons_abs_inert_MB = @constraint(model, [s in S, n in 1:N_abs, k in (N2, O2)],
        cons_inert_mb(abs_y, abs_V, abs_y_Vi, s, n, k) == 0)

    # Reaction equilibrium constraints
    cons_abs_eq1 = @constraint(model, [s in S, n in 1:N_abs], cons_rxne1(abs_T, abs_log_x, s, n) == 0)
    cons_abs_eq2 = @constraint(model, [s in S, n in 1:N_abs], cons_rxne2(abs_T, abs_log_x, s, n) == 0)
    cons_abs_eq3 = @constraint(model, [s in S, n in 1:N_abs], cons_rxne3(abs_T, abs_log_x, s, n) == 0)
    cons_abs_eq4 = @constraint(model, [s in S, n in 1:N_abs], cons_rxne4(abs_T, abs_log_x, s, n) == 0)
    cons_abs_eq5 = @constraint(model, [s in S, n in 1:N_abs], cons_rxne5(abs_T, abs_log_x, s, n) == 0)

    # VLE constraints
    @expression(model, abs_MEA_app[s in S, n in 1:N_abs], compute_stage_apparent_liquid_MEA(abs_x, s, n))
    @expression(model, abs_H2O_app[s in S, n in 1:N_abs], compute_stage_apparent_liquid_H2O(abs_x, s, n))
    @expression(model, abs_CO2_app[s in S, n in 1:N_abs], compute_stage_apparent_liquid_CO2(abs_x, s, n))
    abs_app = (abs_MEA_app, abs_H2O_app, abs_CO2_app)
    He_CO2_expr = get_H.(abs_MEA_app, abs_H2O_app, abs_T) / 101325 # [=] atm·m³/mol
    P_H2O_sat_expr = get_Pvap.(abs_T) / 101325 # [=] atm
    MW_avg_expr = [sum(abs_x[s, n, species_index[c]] * MW[c] for c in keys(species_index)) * 1e-3 for s in S, n in 1:N_abs]

    cons_abs_H2O_VLE = @constraint(model, [s in S, n = 1:N_abs],
        cons_h2o_vle(P_abs, abs_y, abs_x, P_H2O_sat_expr, s, n) == 0)
    cons_abs_CO2_VLE = @constraint(model, [s in S, n = 1:N_abs],
        cons_co2_vle(P_abs, abs_y, abs_x, He_CO2_expr, MW_avg_expr, density, s, n) == 0)


    # Summation constraint (mole fraction)
    cons_abs_sum_x = @constraint(model, [s in S, n=1:N_abs], mole_fraction_sum(abs_x, s, n, nC) == 0)
    cons_abs_sum_y = @constraint(model, [s in S, n=1:N_abs], mole_fraction_sum(abs_y, s, n, nC) == 0)

    # Unit check: abs_V [mol/s] * R [Pa*m^3/(mol*K)] * abs_T [K] / P [Pa] = [m^3/s].
    # Unit check: u_max_G [m/s] * pi * abs_d^2 / 4 [m^2] = [m^3/s].
    cons_abs_superficial_velocity_limit = @constraint(model, [s in S, n=1:N_abs],
        cons_superficial_velocity_limit(abs_V, abs_T, abs_d, P_abs, u_max_G, s, n) <= 0)

    # Liquid and vapor enthalpy
    @expression(model, abs_H_liq_in[s in S], compute_inlet_liquid_enthalpy(abs_app_Li, abs_T_L_in, s))
    @expression(model, abs_H_liq_out[s in S, n in 1:N_abs], compute_stage_liquid_enthalpy(abs_app, abs_T, abs_L, s, n))
    @expression(model, abs_H_vap_in[s in S], compute_inlet_vapor_enthalpy(abs_y_Vi, abs_T_V_in, s))
    @expression(model, abs_H_vap_out[s in S, n in 1:N_abs], compute_stage_vapor_enthalpy(abs_y, abs_T, abs_V, s, n))

    cons_abs_energy = nothing
    cons_abs_energy_1 = nothing
    cons_abs_energy_N = nothing
    cons_abs_energy_n = nothing

    # Energy balance constraints
    if N_abs == 1
        cons_abs_energy = @constraint(model, [s in S],
            cons_energy_single(abs_H_liq_in, abs_H_vap_in, abs_Q, abs_H_liq_out, abs_H_vap_out, s) == 0)
    else
        cons_abs_energy_1 = @constraint(model, [s in S],
            cons_energy_first(abs_H_liq_out, abs_H_vap_in, abs_Q, abs_H_vap_out, s) == 0)

        cons_abs_energy_N = @constraint(model, [s in S],
            cons_energy_last(abs_H_liq_in, abs_H_vap_out, abs_Q, abs_H_liq_out, s, N_abs) == 0)

        if N_abs > 2
        cons_abs_energy_n = @constraint(model, [s in S, n in 2:(N_abs - 1)],
                cons_energy_middle(abs_H_liq_out, abs_H_vap_out, abs_Q, s, n) == 0)
        end
    end


    # --- Initialization ---
    # if the entering stream have a fixed value, use that, otherwise use to its start‐value
    abs_L0 = [inlet_data.L_in_s[s] for s in S] # liquid flow to absorber (Use provided data)
    abs_V0 = [check_val(abs_V_in[s]) for s in S] # vapor flow to absorber
    abs_TV0 = [check_val(abs_T_V_in[s]) for s in S] # vapor T to absorber
    init_LVT(model, abs_L, abs_V, abs_T, abs_L0, abs_V0, abs_TV0, S, N_abs) #Initialize absorber stage L, V, T

    abs_x0 = inlet_data.xL # liquid mole fraction to absorber (Use provided data)
    abs_y0 = [check_val.(abs_y_in[s, :]) for s in S]
    abs_log_x0 = [liquid_log_start(c, abs_x0[s][c]) for s in 1:length(S), n in 1:N_abs, c in 1:nC]
    abs_log_y0 = [vapor_log_start(c, abs_y0[s][c]) for s in 1:length(S), n in 1:N_abs, c in 1:nC]
    init_mole_frac(model, abs_log_x, abs_log_y, abs_log_x0, abs_log_y0, S, N_abs) #Initialize absorber stage liquid and vapor mole fractions
    for s in S, n in 1:N_abs, r in 1:NREACTIONS
        set_start_value(abs_xi[s, n, r], 0.0)
    end
    #endregion

    variables = (; abs_L, abs_V, abs_T, abs_log_x, abs_log_y, abs_xi, abs_d)
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
    return (model=model, variables=variables, constraints=constraints)
end

function build_absorber_flash_model(inlet_data, N_abs=1, Q=0.0)
    build_timing = @timed _build_absorber_flash_model(inlet_data, N_abs, Q)
    built = build_timing.value
    return merge(built, (; build_timing=build_timing))
end

function solve_jump_flash(built; verbose=false)
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

    if !verbose
        set_silent(model)
    end

    solve_timing = @timed optimize!(model)

    abs_d_val = value.(abs_d)
    abs_L_val = value.(abs_L)
    abs_V_val = value.(abs_V)
    abs_T_val = value.(abs_T)
    abs_log_x_val = value.(abs_log_x)
    abs_log_y_val = value.(abs_log_y)
    abs_xi_val = value.(abs_xi)
    vol_flow = abs_V_val .* R .* abs_T_val ./ (P_abs * 101325.0)
    residual = vol_flow .- u_max_G * pi * only(abs_d_val)^2 / 4.0
    d_required = sqrt(4.0 * maximum(vol_flow) / (u_max_G * pi))
    d_actual = only(abs_d_val)
    objective = objective_value(model)
    status = termination_status(model)
    iterations = JuMP.MOI.get(model, JuMP.MOI.BarrierIterations())
    species_names = String.(first.(sort(collect(species_index), by=last)))

    if verbose
        println("Termination status: ", status)
        println("Primal status: ", primal_status(model))
        println("Objective value (absorber diameter): ", objective)
        println("max signed g: ", maximum(residual))
        println("min signed g: ", minimum(residual))
        println("d required by max flow: ", d_required)
        println("d actual: ", d_actual)
        println("slack: ", d_actual - d_required)
    end

    return (;
        backend="JuMP/Ipopt",
        status=string(status),
        objective=Float64(objective),
        iterations=iterations,
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
            iterations=iterations,
        ),
        solve_timing=solve_timing,
    )
end
