% rAA(1) for the iteration u_{k+1} = q(u_k), with residual r_k = u_k - q(u_k)
% NOTE: r is the residual vector for the second last iterate, and similarly 
% rnorms(end) = ||r|| 
% rtol is a relative residual halting tolerance and is optional
function [rnorms, r] = rAA1(q, u0, maxiter, rtol)

    rnorms    = [];
    u_old     = u0;
    u_current = u0; % Place holder
    u_new     = u0; % Place holder

    % We have u_old and u_current, and from these we compute u_new.
    for iter = 1:maxiter 
    
        % Apply underlying fixed-point iteration
        if mod(iter, 2) == 1
            % current <- new
            u_current = u_new;
    
            % u_new = q(u_current)
            u_new = q(u_current);
            % We need to save q_current
            q_current = u_new;
    
            % Evaluate residual for u_current.
            % Note: r_current = u_current - q(u_current)
            r_current = u_current - q_current;
    
        % Take regular AA(1) step
        else
            % old <- current 
            u_old     = u_current; 
            q_old     = q_current;
            r_old     = r_current;
    
            % current <- new
            u_current = u_new;
    
            % Evaluate q(u_current)
            q_current = q(u_current);
            % Evaluate residual for u_current
            % Note: r_current = u_current - q(u_current)
            r_current = u_current - q_current;
    
            % Solve LSQ problem for beta coefficient.
            beta  = -(r_current - r_old) \ r_current;
    
            % Compute new AA iterate.
            u_new = q_current + beta * ( q_current - q_old );
        end

        r = r_current;
    
        rnorms = [rnorms; norm(r_current)];

        % Halt if relative residual tolerance met
        if iter > 1 && nargin == 4
            if rnorms(end)/rnorms(1) < rtol
                break
            end
        end
    end
end