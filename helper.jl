import JuMP
using Printf

function print_timing_summary(label, timing)
    @printf("%s elapsed time: %.6f seconds\n", label, timing.time)
    @printf("%s allocations: %.3f MiB\n", label, timing.bytes / 1024^2)
    @printf("%s GC time: %.6f seconds (%.2f%%)\n",
        label,
        timing.gctime,
        timing.time == 0.0 ? 0.0 : 100.0 * timing.gctime / timing.time,
    )

    if hasproperty(timing, :compile_time)
        @printf("%s compilation time: %.6f seconds (%.2f%%)\n",
            label,
            timing.compile_time,
            timing.time == 0.0 ? 0.0 : 100.0 * timing.compile_time / timing.time,
        )
    else
        println(label, " compilation time: unavailable in this Julia version")
    end

    if hasproperty(timing, :recompile_time)
        @printf("%s recompilation time: %.6f seconds (%.2f%%)\n",
            label,
            timing.recompile_time,
            timing.time == 0.0 ? 0.0 : 100.0 * timing.recompile_time / timing.time,
        )
    end
end


# ========== Shared species and problem data ==========
const species_index = Dict(:MEA => 1, :MEA2p => 2, :Hp => 3, :Carbamate => 4, :HCO3m => 5, :OHm => 6, :CO3_2m => 7, :CO2 => 8, :H2O => 9, :N2 => 10, :O2 => 11)
const species_react = [:MEA, :MEA2p, :Hp, :Carbamate, :HCO3m, :OHm, :CO3_2m, :CO2, :H2O]
const INERTS = [:N2, :O2]
const NONVOLATILES = [:MEA, :MEA2p, :Hp, :Carbamate, :HCO3m, :OHm, :CO3_2m]

const MW = Dict(
    :MEA => 61.08,
    :MEA2p => 62.08,
    :Hp => 1.008,
    :Carbamate => 121.11,
    :HCO3m => 61.016,
    :OHm => 17.007,
    :CO3_2m => 60.01,
    :CO2 => 44.01,
    :H2O => 18.015,
    :N2 => 28.0134,
    :O2 => 31.9988
)

# Stoichiometry matrix for reactions (9 liquid species x 5 reactions)
# Species order: MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m, CO2, H2O
# Reactions: 1: MEA protonation, 2: HCO3- dissociation, 3: H2O dissociation, 4: Carbamate formation, 5: CO2 hydration
stoich = [1 0 0 1 0;    # MEA
    -1 0 0 0 0;    # MEA2p (MEAH+)
    1 1 1 0 1;    # Hp
    0 0 0 -1 0;    # Carbamate (MEACOO-)
    0 -1 0 1 1;    # HCO3m
    0 0 1 0 0;    # OHm
    0 1 0 0 0;    # CO3_2m
    0 0 0 0 -1;    # CO2
    0 0 -1 -1 -1]    # H2O

const MEA = species_index[:MEA]
const MEA2p = species_index[:MEA2p]
const Hp = species_index[:Hp]
const Carbamate = species_index[:Carbamate]
const HCO3m = species_index[:HCO3m]
const OHm = species_index[:OHm]
const CO3_2m = species_index[:CO3_2m]
const CO2 = species_index[:CO2]
const H2O = species_index[:H2O]
const N2 = species_index[:N2]
const O2 = species_index[:O2]

const NSPECIES = length(species_index)
const NREACTIONS = size(stoich, 2)

# ========== Shared property calculations ==========
const LNK_A = (-3.038325, 216.049, 132.899, -0.52135, 231.465)
const LNK_B = (-7008.357, -12431.7, -13445.9, -2545.53, -12092.1)
const LNK_C = (0.0, -35.4819, -22.4773, 0.0, -36.7816)
const LNK_D = (-0.00313489, 0.0, 0.0, 0.0, 0.0)

function get_lnK(T)
    lnK = ntuple(
        r -> LNK_A[r] + LNK_B[r] / T + LNK_C[r] * log(T) + LNK_D[r] * T,
        NREACTIONS,
    )
    return lnK
end

function lnK(T, r)
    return LNK_A[r] + LNK_B[r] / T + LNK_C[r] * log(T) + LNK_D[r] * T
end

### Vapor pressure of Water [=] Pa (From Akula et al. Equation S30)
function get_Pvap(des_T)
    A = 72.55
    B = -7206.7
    C = -7.1385
    D = 4.05e-6

    return exp(A + B / des_T + C * log(des_T) + D * des_T^2)
end

### Henry's constant of des_CO2 (Simplified from Akula et al. Equation S31)
function get_H(x_MEA, x_H2O, T_L)
    w_MEA = (x_MEA * 61.08) / (x_MEA * 61.08 + x_H2O * 18.01528)
    w_H2O = 1.0 - w_MEA

    t = T_L - 273.15
    α_MW = 1.70981 + 0.03972 * t - 4.3e-4 * t^2 - 2.20377 * w_H2O

    He_CO2_H2O = 3.52e6 * exp(-2113.0 / T_L)
    base = 0.0289738 * exp(935.0 / T_L)
    N2O_ratio_term = exp(w_MEA * log(base))
    mixing_term = exp(w_MEA * w_H2O * α_MW)

    return He_CO2_H2O * N2O_ratio_term * mixing_term
end

# ========== Shared enthalpy calculations ==========
R = 8.314
const C_MEA = (2.6161, 3.706e-3, 3.787e-6, 0.0, 0.0)
const C_H2O = (4.2107, -1.696e-3, 2.568e-5, -1.095e-7, 3.038e-10)
const A_CO2 = (5.457, 1.045e-3, -1.157e5)
const A_H2O = (3.47, 1.45e-3, 1.21e4)
const A_N2 = (3.28, 5.93e-4, 4.0e3)
const A_O2 = (3.639, 5.06e-4, -2.27e4)

function get_Hvap(T)
    c = (5.66e4, 0.61204, -0.6257, 0.3988)
    Tc = 647.096
    Tr = T / Tc
    return c[1] * (1 - Tr)^(c[2] + c[3] * Tr + c[4] * Tr^2)
end

function H_H2O_liq(T)
    return MW[:H2O] * (
        C_H2O[1] * ((T - 273.15) - 25) +
        C_H2O[2] * ((T - 273.15)^2 - 25^2) / 2 +
        C_H2O[3] * ((T - 273.15)^3 - 25^3) / 3 +
        C_H2O[4] * ((T - 273.15)^4 - 25^4) / 4 +
        C_H2O[5] * ((T - 273.15)^5 - 25^5) / 5
    )
end

function H_MEA_liq(T)
    return MW[:MEA] * (
        C_MEA[1] * ((T - 273.15) - 25) +
        C_MEA[2] * ((T - 273.15)^2 - 25^2) / 2 +
        C_MEA[3] * ((T - 273.15)^3 - 25^3) / 3
    )
end

H_H2O_vap(T) = R * (A_H2O[1] * (T - 298.15) + A_H2O[2] * ((T^2 - 298.15^2) / 2) - A_H2O[3] * (1 / T - 1 / 298.15)) + get_Hvap(298.15)
H_CO2_vap(T) = R * (A_CO2[1] * (T - 298.15) + A_CO2[2] * ((T^2 - 298.15^2) / 2) - A_CO2[3] * (1 / T - 1 / 298.15))
H_N2_vap(T) = R * (A_N2[1] * (T - 298.15) + A_N2[2] * ((T^2 - 298.15^2) / 2) - A_N2[3] * (1 / T - 1 / 298.15))
H_O2_vap(T) = R * (A_O2[1] * (T - 298.15) + A_O2[2] * ((T^2 - 298.15^2) / 2) - A_O2[3] * (1 / T - 1 / 298.15))

# ========== Shared scalar expressions ==========
apparent_MEA(x_MEA, x_MEA2p, x_Carbamate) =
    x_MEA + x_MEA2p + x_Carbamate

apparent_H2O(x_H2O, x_Hp, x_OHm, x_HCO3m) =
    x_H2O + x_Hp + x_OHm + x_HCO3m

apparent_CO2(x_CO2, x_HCO3m, x_CO3_2m, x_Carbamate) =
    x_CO2 + x_HCO3m + x_CO3_2m + x_Carbamate

liquid_enthalpy_from_apparent(mea_amount, h2o_amount, co2_amount, T) =
    mea_amount * H_MEA_liq(T) +
    h2o_amount * H_H2O_liq(T) +
    co2_amount * (H_CO2_vap(T) - 84000.0)

vapor_enthalpy_from_amounts(h2o_amount, co2_amount, n2_amount, o2_amount, T) =
    h2o_amount * H_H2O_vap(T) +
    co2_amount * H_CO2_vap(T) +
    n2_amount * H_N2_vap(T) +
    o2_amount * H_O2_vap(T)

compute_inlet_apparent_liquid_MEA(x, s) =
    apparent_MEA(x[s, MEA], x[s, MEA2p], x[s, Carbamate])

compute_inlet_apparent_liquid_H2O(x, s) =
    apparent_H2O(x[s, H2O], x[s, Hp], x[s, OHm], x[s, HCO3m])

compute_inlet_apparent_liquid_CO2(x, s) =
    apparent_CO2(x[s, CO2], x[s, HCO3m], x[s, CO3_2m], x[s, Carbamate])

compute_inlet_apparent_liquid(x, s) = (
    compute_inlet_apparent_liquid_MEA(x, s),
    compute_inlet_apparent_liquid_H2O(x, s),
    compute_inlet_apparent_liquid_CO2(x, s),
)

compute_stage_apparent_liquid_MEA(x, s, n) =
    apparent_MEA(x[s, n, MEA], x[s, n, MEA2p], x[s, n, Carbamate])

compute_stage_apparent_liquid_H2O(x, s, n) =
    apparent_H2O(x[s, n, H2O], x[s, n, Hp], x[s, n, OHm], x[s, n, HCO3m])

compute_stage_apparent_liquid_CO2(x, s, n) =
    apparent_CO2(x[s, n, CO2], x[s, n, HCO3m], x[s, n, CO3_2m], x[s, n, Carbamate])

compute_stage_apparent_liquid(x, s, n) = (
    compute_stage_apparent_liquid_MEA(x, s, n),
    compute_stage_apparent_liquid_H2O(x, s, n),
    compute_stage_apparent_liquid_CO2(x, s, n),
)

compute_inlet_liquid_enthalpy(app_liquid, T_liquid, s) =
    liquid_enthalpy_from_apparent(app_liquid[1][s], app_liquid[2][s], app_liquid[3][s], T_liquid[s])

compute_stage_liquid_enthalpy(app_liquid, T_liquid, L_liquid, s, n) =
    liquid_enthalpy_from_apparent(
        app_liquid[1][s, n] * L_liquid[s, n],
        app_liquid[2][s, n] * L_liquid[s, n],
        app_liquid[3][s, n] * L_liquid[s, n],
        T_liquid[s, n],
    )

compute_inlet_vapor_enthalpy(y_amount, T_vapor, s) =
    vapor_enthalpy_from_amounts(y_amount[s, H2O], y_amount[s, CO2], y_amount[s, N2], y_amount[s, O2], T_vapor[s])

compute_stage_vapor_enthalpy(y_vapor, T_vapor, V_vapor, s, n) =
    vapor_enthalpy_from_amounts(
        y_vapor[s, n, H2O] * V_vapor[s, n],
        y_vapor[s, n, CO2] * V_vapor[s, n],
        y_vapor[s, n, N2] * V_vapor[s, n],
        y_vapor[s, n, O2] * V_vapor[s, n],
        T_vapor[s, n],
    )

const LOG_MOLE_FLOOR = exp(-50.0)

safe_log_mole_fraction(z) = log(max(Float64(z), LOG_MOLE_FLOOR))

liquid_log_start(k, z) =
    k in (N2, O2) ? -100.0 : safe_log_mole_fraction(z)

vapor_log_start(k, z) =
    k in (MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m) ? -100.0 : safe_log_mole_fraction(z)

liquid_log_lower(k) = k in (N2, O2) ? -100.0 : -50.0
liquid_log_upper(k) = k in (N2, O2) ? -100.0 : 0.0

vapor_log_lower(k) =
    k in (MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m) ? -100.0 : -50.0

vapor_log_upper(k) =
    k in (MEA, MEA2p, Hp, Carbamate, HCO3m, OHm, CO3_2m) ? -100.0 : 0.0

function helper_self_check()
    Ts = (298.15, 313.15, 350.0)
    for T in Ts
        @assert all(isfinite, get_lnK(T))
        @assert isfinite(lnK(T, 1))
        @assert isfinite(get_Pvap(T))
        @assert isfinite(get_H(0.1, 0.9, T))
        @assert isfinite(get_Hvap(T))
        @assert isfinite(H_H2O_liq(T))
        @assert isfinite(H_MEA_liq(T))
        @assert isfinite(H_H2O_vap(T))
        @assert isfinite(H_CO2_vap(T))
        @assert isfinite(H_N2_vap(T))
        @assert isfinite(H_O2_vap(T))
    end
    return true
end

# ========== JuMP-specific initialization helpers ==========
check_val(v) = JuMP.is_fixed(v) ? JuMP.fix_value(v) : JuMP.start_value(v)

"""Helper to initialize L, V, T for single- or multi-stage units with the provided values."""
function init_LVT(model, L, V, T, L0, V0, T0, S, N_stages)
    for s in S
        for n in 1:N_stages
            if L0 !== nothing
                JuMP.set_start_value(L[s, n], L0[s])
            end
            if V0 !== nothing
                JuMP.set_start_value(V[s, n], V0[s])
            end
            if T0 !== nothing
                JuMP.set_start_value(T[s, n], T0[s])
            end
        end
    end
end

"""Helper to initialize mole fractions for single- or multi-stage units with the provided values."""
function init_mole_frac(model, log_x, log_y, log_x0, log_y0, S, N_stages)
    for s in S
        for n in 1:N_stages
            for c in 1:length(species_index)
                if log_x0 !== nothing
                    JuMP.set_start_value(log_x[s, n, c], log_x0[s, n, c])
                end
                if log_y0 !== nothing
                    JuMP.set_start_value(log_y[s, n, c], log_y0[s, n, c])
                end
            end
        end
    end
end
