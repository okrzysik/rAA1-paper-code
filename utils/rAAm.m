% rAA(m) for the iteration u_{k+1} = q(u_k), with residual r_k = u_k - q(u_k)
% 
% R is a matrix with columns the residual vectors for all iterates, except
% the last. rnorms is the vector whose elements are the norm of the columns 
% in this matrix.
% 
% Note that this is not a memory-efficient implementation.
%
% If u_{k+1} = q(u_k) is the underlying fixed-point iteration then the AA
% iteration used here is for k>=0
%   u_{k+1} = q(u_k) + sum_{i=1}^{m_k} beta_k^i * ( q(u_k) - q(u_{k-i}) ),
% where m_k = mod(k, m+1) is the memory parameter.

function [rnorms, R] = rAAm(q, u0, maxiter, m)

    rnorms = [];

    u0 = u0(:);
    [n, ~] = size(u0);

    U = zeros(n, maxiter+1);
    Q = zeros(n, maxiter);
    R = zeros(n, maxiter);

    U(:, 1) = u0;
    
    % We have u_k and we're going to compute u_k+1
    for k = 1:maxiter
    
        % Fixed-point function applied to previous point
        Q(:, k)   = q(U(:, k)); % q(u_k)

        % Residual in previous point
        R(:, k) = U(:, k) - Q(:, k);
        rnorms = [rnorms; norm(R(:, k))];

        % Window size. Note that k-1 here not k since we loop from k = 1
        % and not k = 0 due to MATLAB's 1's based indexing.
        mk = mod(k-1, m+1);

        % Just do fixed point
        if mk == 0    
            U(:, k+1) = Q(:, k);    % u_{k+1} = q(u_k)
            
        % Actually do an AA step where we add in a weighted sum of q
        % applied to past iterates. The LSQ problem for the coefficients
        % beta is B*beta = -r_k
        else
    
            % Assemble system matrix
            B = zeros(n, mk);
            for i = 1:mk
                B(:, i) = R(:, k) - R(:, k-i);
            end

            % Solve linear system
            beta = -B \ R(:, k);
            
            % Compute updated solution
            U(:, k+1) = Q(:, k);    % u_{k+1} = q(u_k)
            for i = 1:mk            % Add in the sum
                U(:, k+1) = U(:, k+1) + beta(i) * ( Q(:, k) - Q(:, k-i) );
            end

        end

    end

    % % Compute q and r at the last point
    % Q(:, k+1) = q(U(:, k+1)); % q(u_k+1)
    % R(:, k+1) = U(:, k+1) - Q(:, k+1);
    % rnorms = [rnorms; norm(R(:, k+1))];
end