function main_model4_compare  %fminsearch
    
    % Load real prevalence data
    years = (1990:2022)';
    prev_data = [5.5 6.8 8.1 9.0 9.6 9.9 10.0 9.8 9.5 9.2 ...
        8.7 8.3 7.8 7.3 6.9 6.6 6.3 6.1 5.9 5.7 ...
        5.6 5.5 5.4 5.2 5.1 4.9 4.8 4.6 4.4 4.2 ...
        4.0 3.9 3.7]';

    % Fixed parameters
    mu = 1/35;
    N0 = 22892651; % Kenya pop 1990
    nG = 4; nY = 3;

    % Transmission parameters
    pHIV = [0.0107, 0.0008, 0.0042];
    eta = [10.6, 11.0, 7.1] * 12;
    tau = [48 36 24 12; 36 24 12 6; 24 12 6 1; 12 6 1 0.25]/12;
    t_alpha = zeros(nG,nG,nY);
    for alpha = 1:nY
        for i = 1:nG
            for j = 1:nG
                t_alpha(i,j,alpha) = 1 - ((1 - pHIV(alpha))^(eta(alpha) * tau(i,j)));
            end
        end
    end
    w = [1/(2.5/12), 1/7.59, 1/2];
    g = ones(1,nG)/nG;

    %% ----------- Fminsearch (Nelder–Mead) -------------
    fprintf('\n========== Nelder–Mead (fminsearch) ===========\n');
    guess = [0.1, 0.5, 2.0, 10.0, 0.10, 1.0, 0.7, 0.055];
    opt1 = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',500,'MaxFunEvals',3000);
    [best_fmin, err_fmin] = fminsearch(@(p) cost_function0(p, years, prev_data, mu, t_alpha, w, g, N0), guess, opt1);

    rhoS_fmin = best_fmin(1:4)';
    f_fmin = best_fmin(5);
    h_fmin = best_fmin(6);
    e_fmin = best_fmin(7);
    init_fmin = best_fmin(8);

    fprintf('\nFminsearch results:\n');
    fprintf('rho_S = [%.4f %.4f %.4f %.4f]\n', rhoS_fmin);
    fprintf('f = %.4f, h = %.4f, e = %.4f, init = %.4f\n', f_fmin, h_fmin, e_fmin, init_fmin);
    fprintf('Error = %.6f\n', err_fmin);

    [t_fmin, prev_fmin] = simulate_prevalence(years, rhoS_fmin, f_fmin, h_fmin, e_fmin, init_fmin, mu, t_alpha, w, g, N0);

    %% ----------- Plot Result --------------------------
    figure;
    plot(years, prev_data, 'ko', 'LineWidth', 1.3, 'DisplayName', 'Observed'); hold on;
    plot(t_fmin, prev_fmin, 'b-', 'LineWidth', 1.6, 'DisplayName', 'Fminsearch Fit');
    xlabel('Year'); ylabel('Prevalence (%)');
    title('Fitting with fminsearch (Nelder–Mead)');
    legend('Location', 'best');
    grid on;
end

function [t, prevalence] = simulate_prevalence(years, rho_S, f, h, e, init, mu, t_alpha, w, g, N0)
    nG = 4; nY = 3;
    rho_B = rho_S;
    rho_Y = repmat(rho_S, 1, nY);
    Y0_total = N0 * init;
    Y0 = zeros(4,nY);
    Y0(:,1) = g' * Y0_total;
    S0 = N0 * g' - Y0(:,1);
    B0 = zeros(4,1);
    y0 = [S0; B0; reshape(Y0',[],1)];

    [t, y] = ode15s(@(t,y) model4_function(t, y, mu, rho_S, rho_B, rho_Y, t_alpha, w, h, f, g, e), ...
                    [min(years) max(years)], y0);

    S = y(:,1:4);
    B = y(:,5:8);
    Yp = y(:,9:end);
    total_inf = zeros(size(t));
    total_pop = zeros(size(t));
    for k = 1:length(t)
        Y = reshape(Yp(k,:)', [nY,4])';
        total_inf(k) = sum(Y(:));
        total_pop(k) = sum(S(k,:)) + sum(B(k,:)) + total_inf(k);
    end
    prevalence = 100 * (total_inf ./ total_pop);
end
