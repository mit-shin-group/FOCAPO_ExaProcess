function reaction_extent(stoich, xi, s, n, k)
    extent = sum(stoich[k, r] * xi[s, n, r] for r in 1:NREACTIONS)
    return extent
end

function cons_reactive_mb_single(x, y, L, V, y_Vi, x_Li, stoich, xi, s, k)
    residual = x[s, 1, k] * L[s, 1] +
        y[s, 1, k] * V[s, 1] -
        y_Vi[s, k] -
        x_Li[s, k] -
        reaction_extent(stoich, xi, s, 1, k)
    return residual
end

function cons_reactive_mb_first(x, y, L, V, y_Vi, stoich, xi, s, k)
    residual = x[s, 1, k] * L[s, 1] +
        y[s, 1, k] * V[s, 1] -
        y_Vi[s, k] -
        x[s, 2, k] * L[s, 2] -
        reaction_extent(stoich, xi, s, 1, k)
    return residual
end

function cons_reactive_mb_last(x, y, L, V, x_Li, stoich, xi, s, k, N_abs)
    residual = x[s, N_abs, k] * L[s, N_abs] +
        y[s, N_abs, k] * V[s, N_abs] -
        y[s, N_abs - 1, k] * V[s, N_abs - 1] -
        x_Li[s, k] -
        reaction_extent(stoich, xi, s, N_abs, k)
    return residual
end

function cons_reactive_mb_middle(x, y, L, V, stoich, xi, s, n, k)
    residual = x[s, n, k] * L[s, n] +
        y[s, n, k] * V[s, n] -
        y[s, n - 1, k] * V[s, n - 1] -
        x[s, n + 1, k] * L[s, n + 1] -
        reaction_extent(stoich, xi, s, n, k)
    return residual
end

function cons_inert_mb(y, V, y_Vi, s, n, k)
    residual = y[s, n, k] * V[s, n] - y_Vi[s, k]
    return residual
end

function cons_rxne1(T, log_x, s, n)
    residual = lnK(T[s, n], 1) +
        log_x[s, n, MEA2p] -
        log_x[s, n, MEA] -
        log_x[s, n, Hp]
    return residual
end

function cons_rxne2(T, log_x, s, n)
    residual = lnK(T[s, n], 2) +
        log_x[s, n, HCO3m] -
        log_x[s, n, CO3_2m] -
        log_x[s, n, Hp]
    return residual
end

function cons_rxne3(T, log_x, s, n)
    residual = lnK(T[s, n], 3) -
        log_x[s, n, OHm] -
        log_x[s, n, Hp]
    return residual
end

function cons_rxne4(T, log_x, s, n)
    residual = lnK(T[s, n], 4) +
        log_x[s, n, Carbamate] -
        log_x[s, n, MEA] -
        log_x[s, n, HCO3m]
    return residual
end

function cons_rxne5(T, log_x, s, n)
    residual = lnK(T[s, n], 5) +
        log_x[s, n, CO2] -
        log_x[s, n, HCO3m] -
        log_x[s, n, Hp]
    return residual
end

function cons_h2o_vle(P, y, x, P_H2O_sat, s, n)
    residual = P * y[s, n, H2O] - x[s, n, H2O] * P_H2O_sat[s, n]
    return residual
end

function cons_co2_vle(P, y, x, He_CO2, MW_avg, density, s, n)
    residual = P * y[s, n, CO2] - He_CO2[s, n] * x[s, n, CO2] * density / MW_avg[s, n]
    return residual
end

function mole_fraction_sum(z, s, n, nC)
    residual = sum(z[s, n, k] for k in 1:nC) - 1.0
    return residual
end

function cons_superficial_velocity_limit(V, T, d, P, u_max_G, s, n)
    residual = V[s, n] * R * T[s, n] / (P * 101325.0) -
        u_max_G * pi * d[1]^2 / 4.0
    return residual
end

function cons_energy_single(H_liq_in, H_vap_in, Q, H_liq_out, H_vap_out, s)
    residual = H_liq_in[s] + H_vap_in[s] + Q[s, 1] -
        H_liq_out[s, 1] - H_vap_out[s, 1]
    return residual
end

function cons_energy_first(H_liq_out, H_vap_in, Q, H_vap_out, s)
    residual = H_liq_out[s, 2] + H_vap_in[s] + Q[s, 1] -
        H_liq_out[s, 1] - H_vap_out[s, 1]
    return residual
end

function cons_energy_last(H_liq_in, H_vap_out, Q, H_liq_out, s, N_abs)
    residual = H_liq_in[s] + H_vap_out[s, N_abs - 1] + Q[s, N_abs] -
        H_liq_out[s, N_abs] - H_vap_out[s, N_abs]
    return residual
end

function cons_energy_middle(H_liq_out, H_vap_out, Q, s, n)
    residual = H_liq_out[s, n + 1] + H_vap_out[s, n - 1] + Q[s, n] -
        H_liq_out[s, n] - H_vap_out[s, n]
    return residual
end
