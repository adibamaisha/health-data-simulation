function err = cost_function0(p, years, prev_data, mu, t_alpha, w, g, N0)
    % Unpack parameters
    rho_S = p(1:4)';
    f = p(5);
    h = p(6);
    e = p(7);
    init = p(8);

    % Keep values positive
    if any(p <= 0)
        err = Inf;
        return;
    end

    try
        [t, prevalence] = model4_prevalence_only(years, rho_S, f, h, e, init, mu, t_alpha, w, g, N0);
        err = sum((prev_data - prevalence).^2);
    catch
        err = Inf;
    end
end
