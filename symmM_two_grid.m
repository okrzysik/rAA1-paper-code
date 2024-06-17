% Uses rAA(1) to accelerate a two-grid algorithm for solving a 2D Poisson
% problem. 
%
% The code iteratively solves the linear system A*x = f by repeated calls 
% to a multigrid function like
%   x <- multigrid(x, b) until convergence.
% I.e., the Picard iteration function is q(x) = multigrid(x, b)

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear all
close all
rng(1)

test_matrix = 3; % Just for labelling things.
save_M_eigs = ~true; % Save eigenvalues of M? For plotting.

save_fig = ~true;


% Add to plots the convergence factors computed from theory. This requires
% assembling the multigrid iteration matrix.
compute_theory_bound = true;


num_x0  = 30;  % Number of different x0 to try
maxiter = 300; % Number of iterations

%% Problem and multigrid parameters
nx = 16; % There are (n-1)^2 interior grid points.
n = (nx-1)^2;
h = 1/nx; % Step size.
options.omega = 4/5; % Weighted Jacobi scheme's weighting factor.
num_cycles = 20;

A = matPoisson(nx); % Coefficient matrix of the matrix system Au = f.
A_coarse = matPoisson(nx/2);
f = zeros(n, 1); % RHS of linear system f.
%f = randn(n, 1); % RHS of linear system f.
x = randn(n, 1); % Use a random initial guess for the approximate solution v.

% Interpolation
P = interpolate(nx);
% Restriction
R = 0.25*P';


%% Run the basic two-grid algorithm 
% See what the convergence factor is, and check if this convergence factor 
% is consistent with that from the Picard iteration.
rnorms = [];
% Run two-grid cycles
fprintf('Basic two-grid method:\n')
for iter = 1:num_cycles

    x = two_grid(x, f, A, A_coarse, P, R, options);

    rnorms = [rnorms; norm(f - A*x)];

    if iter == 1
        fprintf('||r||_%d = %.2e\n', iter-1, rnorms(iter))
    else
        fprintf('||r||_%d = %.1e, ||r||_%d/||r||_%d = %.2e\n', ...
            iter-1, rnorms(end), iter, iter-1, rnorms(end)/rnorms(end-1))
    end
end
fprintf('\n\n')



%% Picard iteration
% For this problem, the iteration function M is one two-grid cycle 
q = @(v) two_grid(v, f, A, A_coarse, P, R, options);

% Iteration numbers for plotting.
k = (1:maxiter-1)';
plot_skip = 5;
k_plot = k(1:plot_skip:end);

% Loop over different x0 vectors
for x0_test = 1:num_x0
    x0 = rand(n, 1) - rand(n, 1);

    [rnorms_rAA, r] = rAA1(q, x0, maxiter);
    rnorms_PI       = PI(q,   x0, maxiter);
   
    % Root-average conv factor
    figure(1)
    plot(k_plot, rnorms_rAA(2:plot_skip:end).^(1./k_plot), '-*', ...
        'MarkerSize', 8)
    hold on
    plot(k_plot, rnorms_PI(2:plot_skip:end).^(1./k_plot),  '-<', ...
        'MarkerSize', 8)
    

    % Some things RE the NEPv
    % % Decompose r into the eigenvectors of M.
    % if maxiter == 0
    %     r = x0;
    % end
    % r     = r / norm(r);
    % gamma = M_evecs \ r;
    % figure(3)
    % semilogy(x0_test * ones(size(gamma)), sort(abs(gamma)), '-x', ...
    %     'MarkerSize', 8, 'MarkerFaceColor', 'auto')
    % hold on

    % Compute some measure of how close the residual is to solving the NEPv1.
    % z = X(r)*r;
    % theta = acos( z' * r / (norm(z)*norm(r)) )
    % 
    % z = (X(r)*r) ./ r;
    % [mean(z), std(z)]
end


%% Assemble two-grid error propagation matrix
% See Trottenberg et al. p. 40 for error propagator formula.
if compute_theory_bound

    % Fine-grid matrix
    A1 = A;
    % Smoothing matrix.
    [~, S] = wJacobi(A1, f, f, options);
    % Interpolation matrix
    P = interpolate(nx);
    % Restriction matrix
    R = 0.25*P';
    % Coarse-grid matrix
    A2 = matPoisson(nx/2);
    
    % Coarse-grid correction matrix
    K = eye( (nx-1)^2 ) - P * (A2 \ R) * A1;
    % Two-grid matrix
    M = S * K * S;

    % Compute orthogonal eigenvectors of M
    [M_evecs, M_eigs] = schur(M);
    M_eigs = eig(M_eigs);
    assert(max(abs(imag(M_eigs))) < 1e-10, 'M eigs have significant imaginary component...')
    M_eigs = real(M_eigs);

    if save_M_eigs
        save(sprintf('./data/symmM_eigs_test%d.mat', test_matrix), 'M_eigs')
    end

    % Compute eigenvalues of M. There's small round off error here so
    % eigenvalues come out with small imaginary part.
    m = sort(real(eig(M)));
    
    % Sanity check: Do M and the q in fact have the same action?
    e = randn((nx-1)^2, 1);
    assert(norm(M*e - q(e), inf) < 1e-10, 'M and q do not have the same action?')
    
    rho_rAA1_theory = NEPv_lambda_max(m)^0.25;
    rho_PI_theory   = max(abs(m));

    figure(1)
    hold on
    plot(k, rho_rAA1_theory*ones(size(k)), '--b', 'LineWidth', 6)
    plot(k, rho_PI_theory*ones(size(k))  , ':k', 'LineWidth', 6)

    % Matrix A
    I = eye((nx-1)^2);
    A = I - M;
    
    % Form the matrix function X(v).
    alpha = @(v) (A*v)' * v / norm(A*v)^2;
    R     = @(v) M*(I - alpha(v)*A) * v;
    X     = @(v) M*(I - alpha(R(v))*A) * M*(I - alpha(v)*A);
end



% Tidy up each plot.
figure(1)
xlabel('$k$')
ylabel('$\varrho_k(\mathbf{r}_0) := \Vert \mathbf{r}_k \Vert^{1/k}$')
ylim([0.85*rho_rAA1_theory, 1.1*rho_PI_theory])
title('test = 3')
box on

if save_fig
    fig_name = sprintf('./figures/rhok-symmM-MG-nx%d', nx);
    figure_saver(gcf, fig_name, false);
end

% Eigenvalue plot
figure(2)
plot(m, 'rx')
title('Eigenvalues $\{ m_i \}$ of $M$')
xlabel('$i$')
axis tight
box on

if save_fig
    fig_name = sprintf('./figures/eig-symmM-MG-nx%d', nx);
    figure_saver(gcf, fig_name, false);
end


% % Plot coefficients of residual in eigenvector basis.
% figure(3)
% title(sprintf('$k = %d$', maxiter))
% ylim([1e-2, 1])
% xlabel('$\mathbf{r}_0$ index')
% ylabel('$| V^{-1} \widehat{\mathbf{r}}_k (\mathbf{r}_0) |$')








