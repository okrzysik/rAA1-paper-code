
% Given eigenvalues m, compute the maximum eigenvalue of the associated
% NEPv1
function lambda = NEPv_lambda_max(m)

    n = numel(m);
    lambda = 0;
    for i = 1:n
        for j = 1:n
            if j == i
                continue
            else
                lambda_temp = lambda_eps_max(m(i), m(j));
                if lambda_temp > lambda
                    lambda = lambda_temp;
                    i_max = i;
                    j_max = j;
                end
            end
        end
    end
end

% The nonlinear eigenvalue lambda maximized over epsilon
function lambda = lambda_eps_max(M1, M2)
    lambda = (abs( M1 .* M2 .* (M2 - M1) ) ./ ( abs(M1.*(M1 - 1)) + abs(M2.*(M2 - 1)) )).^2; 
end