function dydt = model4_function(t, y, mu, rho_S, rho_B, rho_Y, ...
                                 t_alpha, w, h, f, g, e)
                             
    nG = 4;  % Number of risk groups
    nY = 3;  % Number of HIV stages

    S = y(1:nG);
    B = y(nG+1:2*nG);
    Y = reshape(y(2*nG+1:end), [nY, nG])';

  
    N = sum(S + B + sum(Y, 2));

    % Mixing matrix G
    numerators = rho_S(:).*S + rho_B(:).*B + sum(rho_Y .* Y, 2);
    
    denom_G = sum(numerators);
   
    G = zeros(nG);
    delta=eye(4);
    for i = 1:nG
        for j = 1:nG
            G(i,j) = e*delta(i,j) + ((1 - e)*numerators(j)/denom_G);
        end
    end
    
    
    
    % Force of infection
    
    denom_lambda = (rho_S(:).*S) +(rho_B(:).*B) + sum(rho_Y .* Y,2);
    
    lambda = zeros(nG,1);
   
   
    for i = 1:nG
        for j = 1:nG
            lambda(i) = lambda(i) + (rho_S(i) * G(i,j) * ...
                        sum((squeeze(t_alpha(i,j,:))) .* (rho_Y(j,:) .* Y(j,:))') / denom_lambda(j));
        end
    end
   


    inf_S = lambda .* S;
    inf_B = h * lambda .* B;

  
    birth_S = (1 - f) * mu * N * g';
    birth_B = f * mu * N * g';

    % dS and dB
    dS = birth_S - inf_S - mu * S;
    dB = birth_B - inf_B - mu * B;

    % Infected stages
    dY = zeros(nG, nY);
    for i = 1:nG
        dY(i,1) = inf_S(i) + inf_B(i) - (mu + w(1)) * Y(i,1);
        dY(i,2) = w(1) * Y(i,1) - (mu + w(2)) * Y(i,2);
        dY(i,3) = w(2) * Y(i,2) - (mu + w(3)) * Y(i,3);

    end

    dydt = [dS; dB; reshape(dY', [], 1)];
 
end
