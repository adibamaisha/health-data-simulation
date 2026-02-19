function main_model4_compare_GA()
    %% Observed data
    observed_years = 1990:2022;
    observed_prevalence = [ ...
        5.5, 6.8, 8.1, 9.0, 9.6, 9.9, 10.0, 9.8, 9.5, 9.2, ...
        8.7, 8.3, 7.8, 7.3, 6.9, 6.6, 6.3, 6.1, 5.9, 5.7, ...
        5.6, 5.5, 5.4, 5.2, 5.1, 4.9, 4.8, 4.6, 4.4, 4.2, ...
        4.0, 3.9, 3.7];

    %% GA setup
    nVars = 6; % only 6 parameters
    lb = [0.001, 2.0, 0.0, 0.0, 0.0, 0.0]; % lower bounds
    ub = [0.1, 5.0, 1.0, 20.0, 1.0, 0.05]; % upper bounds

    repeatThreshold = 50;  % stop if best cost repeats this many generations
    epsilon = 1e-4;        % absolute change threshold for early stop
    useParallel = true;

    %  initial guess 
   x0 = [  0.0115    4.1219    0.5246    9.9073    0.9933    0.0466];   % [rho_0, beta, f, h, e, init]

    % Generate diverse initial population around x0
    populationSize = 200;
    perturb = 0.2; % 20% variation
    initialPop = repmat(x0, populationSize, 1) .* (1 + perturb*(2*rand(populationSize,nVars)-1));

    opts_ga = optimoptions('ga', ...
        'Display', 'iter', ...
        'PopulationSize', populationSize, ...
        'MaxGenerations', 500, ...
        'UseParallel', useParallel, ...
        'InitialPopulationMatrix', initialPop, ...
        'OutputFcn', @(options,state,flag) gaEarlyStop(options,state,flag,epsilon,repeatThreshold));

    %% Run GA
    [best_params, best_cost] = ga(@(params) cost_function0(params, observed_years, observed_prevalence), ...
                                  nVars, [], [], [], [], lb, ub, [], opts_ga);

    %% Display results
    fprintf('\n---------------------------\n');
    fprintf('Best Parameters Found by GA:\n');
    disp(best_params);
    fprintf('Best Cost = %.4f\n', best_cost);
    fprintf('---------------------------\n');

    %% Simulate model
    [calendar_years, prevalence_model] = model4_prevalence_only(best_params, 2022);

    %% Plot
    figure;
    plot(observed_years, observed_prevalence, 'ro-', 'LineWidth', 2); hold on;
    plot(calendar_years, prevalence_model, 'b-', 'LineWidth', 2);
    xline(2022, '--k', 'LineWidth', 1.2);
    xlabel('Year'); ylabel('HIV Prevalence (%)');
    legend('Observed', 'Model (Best Fit)', 'Projection End');
    title('HIV Prevalence in Kenya (GA Calibration)');
    grid on;
end

%% ----------------- GA Early Stop Function -----------------
function [state, options, optchanged] = gaEarlyStop(options,state,flag,epsilon,repeatThreshold)
    optchanged = false;
    persistent prevBest repeatCount

    switch flag
        case 'init'
            prevBest = Inf;
            repeatCount = 0;

        case 'iter'
            currentBest = state.Best(end);

            % Compute relative improvement for monitoring
            rel_improv = abs(prevBest - currentBest) / max(prevBest, eps);

            % Count repeated best cost
            if abs(currentBest - prevBest) < epsilon
                repeatCount = repeatCount + 1;
            else
                repeatCount = 0;  % reset counter if cost changes
            end

            % Stop if repeated best cost exceeded
            if repeatCount >= repeatThreshold
                fprintf('\nGA early stop: abs change %.6f, rel improvement %.6f, repeated best cost = %d\n', ...
                        abs(prevBest - currentBest), rel_improv, repeatCount);
                state.StopFlag = 'y';
            end

            prevBest = currentBest;

        case 'done'
            clear prevBest repeatCount;
    end
end
