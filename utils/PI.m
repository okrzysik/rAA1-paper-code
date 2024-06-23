% Picard iteration u_{k+1} = q(u_k), with residual r_k = u_k - q(u_k)
% NOTE: r is the residual vector for the second last iterate, and similarly 
% rnorms(end) = ||r|| 
% rtol is a relative residual halting tolerance and is optional
function [rnorms, r, u] = PI(q, u0, maxiter, rtol)
    rnorms = [];
    u_current = u0;
    for iter = 1:maxiter
        u_new   = q(u_current);
        q_current = u_new;
        r_current = u_current - q_current;
        rnorms    = [rnorms; norm(r_current)];

        u_current   = u_new;

        % Halt if relative residual tolerance met
        if iter > 1 && nargin == 4
            if rnorms(end)/rnorms(1) < rtol
                break
            end
        end
    end
    u = u_current;
    r = r_current;
end