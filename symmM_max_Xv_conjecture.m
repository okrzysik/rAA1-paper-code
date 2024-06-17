% Numerically "verify" the conjecture that the max of ||X(v)*v|| / ||v||
% over all v is achieved by the nonlinear eigenvector v of X(v) that
% maximizes the associated nonlinear eigenvalue lambda over those
% characterized by the theorem in the paper.
%
% What this code does is randomly generate a symmetric matrix M from which
% X(v) can be computed. Then, for a given number of randomly chosen vectors
% v the quantity in question ||X(v)*v|| / ||v|| is computed, and it is
% compared to the conjectured upper bound.
%
% There are two types of strategy available for choosing v.
%   1: v is randomly generated
%   2: some vector z is chosen at random, it is then iterated as z <-
%       X(z)*z for some number of iterations, then v <- z.
%
% The second strategy seems to produce values of ||X(v)*v|| / ||v|| that 
% are much closer to the closer to the conjectured maximum. 
%
% NOTE: The distribution of possible values of the quantity in question
% seems sensitive to exactly how M is chosen. E.g., if M is a symmetrized
% version of randn, then the distribution of values seems to be closely
% clustered (perhaps indicating that the convergence factor doesn't depend
% strongly on the initial iterate), but if chosen based on rand then the
% distribution is much less clustered. 


tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear

save_fig = ~true;

rng(3) % Fix seed so that test matrices M are the same for both v 
% strategies and so that same v are chosen initially.

v_strategy = 1;
%v_strategy = 2; num_intermediate_iters = 10;


% Number of matrices M to test.
num_M_tests       = 16;
% Number of v to test per matrix M.
num_v_tests_per_M = 10000;



figure(v_strategy)
clf
hold on

% Test the conjecture on a bunch of randomly chosen matrices M.
for M_idx = 1:num_M_tests

    % Size of the matrix M. 
    n = 5*M_idx;

    % % Dense M, will have norm less than or equal to one
    % M = rand(n, n); 
    % M = M + M';
    % M = M / ((1+rand(1))*norm(M)); 

    % Diagonal M, will norm less than or equal to one
    M = diag(rand(n, 1)) - diag(rand(n, 1));

    % Matrix A
    I = eye(n);
    A = I - M;

    % Form the matrix function X(v).
    alpha = @(v) (A*v)' * v / norm(A*v)^2;
    R     = @(v) M*(I - alpha(v)*A) * v;
    X     = @(v) M*(I - alpha(R(v))*A) * M*(I - alpha(v)*A);

    % Compute the conjectured upper bound
    [ratio_max_conj, v_max_conj] = lambda_max(M);
    ratio_max_conj = ratio_max_conj^0.25;

    % Compute the quantity in question over a bunch of different v.
    ratio = [];
    for v_idx = 1:num_v_tests_per_M

        if v_strategy == 1
            v = rand(n, 1) - rand(n, 1);
            v_new = X(v)*v;

        elseif v_strategy == 2
            v = randn(n, 1);

            for iter = 1:num_intermediate_iters
                if iter > 1; v = v_new; end
                v_new = X(v)*v;
            end
        end

        ratio_temp = ( norm(v_new, 2) / norm(v, 2) )^0.25;

        % Check if conjecture has been violated
        if ratio_temp > ratio_max_conj
            [ratio_temp, ratio_max_conj]
            if ratio_temp > 0.1
                fprintf('Conjecture seemingly violated...\n')
            else
                fprintf('Conjecture seemingly violated, but w/ small conv factor...\n')
            end
        end

        ratio = [ratio; ratio_temp];
    end

    figure(v_strategy)
    plot(M_idx*ones(size(ratio)), ratio, 'gx')
    plot(M_idx, ratio_max_conj, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm')

    % Sanity check: Does the worst-case conjectured vector achieve the
    % worst-case conjectured bound when plugged into the expression.
    %plot(M_idx, (norm(X(v_max_conj)*v_max_conj, 2) / norm(v_max_conj, 2) ).^0.25, 'r*')

    plot(M_idx, norm(M)^0.25, 'b>', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
   

    m = eig(M);
    figure(2*v_strategy+1)
    plot(M_idx*ones(size(m)), m, 'o')
    hold on

end

% Pretty up plot.
figure(v_strategy)
box on
xlim([1 num_M_tests])
ylim([0 1])
xlabel('$\ell$: test index')
ylabel('$\big( \Vert X_{\ell}(\mathbf{v}) \mathbf{v} \Vert / \Vert \mathbf{v} \Vert \big)^{1/4}$')
title(sprintf('$\\mathbf{v}$ strategy = %d', v_strategy))

% Save plot
if save_fig
    fig_name = sprintf('./figures/max-Xv-conjecture-vstrat%d-numv%d', v_strategy, num_v_tests_per_M);
    figure_saver(gcf, fig_name, false);
end

figure(2*v_strategy+1)
box on
title('Eigenvalues $\{ m_i \}_{i = 1}^{n(\ell)}$ of $M_{\ell}$')
xlabel('$\ell$: test index')
%ylim([-1, 1])
xlim([1 num_M_tests])

% Save plot
if save_fig
    fig_name = sprintf('./figures/max-Xv-conjecture-eigM', v_strategy, num_v_tests_per_M);
    figure_saver(gcf, fig_name, false);
end







% Compute the maximum nonlinear eigenvalue of X(v) over the set of those
% characterized in the paper. Also compute the associated maximum nonlinear
% eigenvector.
function [lambda, v_max] = lambda_max(M)

    [n, ~] = size(M);
    [V, D] = eig(M);
    m = diag(D);
    
    i_max = 0;
    j_max = 0;
    
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
    
    eps_max = m(j_max)*(m(i_max) - 1) / ( m(i_max)*(m(j_max) - 1) );
    eps_max = sqrt(abs(eps_max));
    v_max = V(:, i_max) + eps_max * V(:, j_max);

end

% The nonlinear eigenvalue lambda maximized over epsilon
function lambda = lambda_eps_max(M1, M2)
    lambda = (abs( M1 .* M2 .* (M2 - M1) ) ./ ( abs(M1.*(M1 - 1)) + abs(M2.*(M2 - 1)) )).^2; 
end
