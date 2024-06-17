% Uses rAA(1) to accelerate an underlying fixed-point iteration with
% iteration function q(x) = M*x with M a symmetric matrix.
%
% There are two different types of test matrix M, each with differently
% structured eigenvalues. Simply choose which "test_matrix" below that you
% would like.
%
% Note that the distribution of convergence factors obtained can be
% sensitive to how exctly the iterates are initialized. E.g, using randn
% vs. rand.
%
% For the rAA(1) iteration, for each x0, the final residual vector will be 
% decomposed into the basis of eigenvectors of M and the magnitude of the
% coefficients plotted.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear
close all
rng(5)

test_matrix = 1; % Choose me
test_matrix = 2; % Or choose me


save_fig = ~true;
save_M_eigs = ~true; % Save eigenvalues of M? For overlapying plots of eigenvalues of different M

n = 32; % Dimension of system/number of spatial DOFs

num_x0  = 30; % Number of different x0 to try
maxiter = 300; % Number of iterations


% Uncomment the below maxiter lines to re-produce the rk-decomposition
% plots from the supplementary materials.
save_rk_decomposition_fig = ~true; % Save decomposition of of rk into eigenvectors of M?
% maxiter = 1;
% maxiter = 61;
% maxiter = 121;
% maxiter = 181;


% Dense M, will have norm less than or equal to one
if test_matrix == 1
    M = rand(n, n); 
    M = M + M';
    M = M / ((1+rand(1))*norm(M)); 
elseif test_matrix == 2
    M = randn(n, n); 
    M = M + M';
    M = M / ((1+rand(1))*norm(M)); 
end


% Compute orthogonal eigenvectors of M
[M_evecs, M_eigs] = schur(M);
M_eigs = diag(M_eigs);

if save_M_eigs
    save(sprintf('./data/symmM_eigs_test%d.mat', test_matrix), 'M_eigs')
end

% RHS vector
b = zeros(n, 1);

% Iteration function
q = @(x) M*x + b;

% Iteration numbers for plotting.
k = (1:maxiter-1)';
% Plotting junk
plot_skip = 5;
k_plot = k(1:plot_skip:end);
rho_rAA1_min = 1;

% Loop over different starting iterate
for x0_test = 1:num_x0

    % Starting iterate
    x0 = rand(n, 1) - rand(n, 1);

    [rnorms_rAA, r] = rAA1(q, x0, maxiter);
    rnorms_PI       = PI(q, x0, maxiter);
   
    % Root-average conv factor
    figure(1)
    plot(k_plot, rnorms_rAA(2:plot_skip:end).^(1./k_plot), '-*', ...
        'MarkerSize', 9)
    hold on
    plot(k_plot, rnorms_PI(2:plot_skip:end).^(1./k_plot),  '-<', ...
        'MarkerSize', 9)
    
    rho_rAA1_min = min(rho_rAA1_min, rnorms_rAA(end));


    % Decompose r into the eigenvectors of M.
    if maxiter == 0
        r = r0;
    end
    r     = r / norm(r);
    gamma = M_evecs \ r;
    figure(3)
    semilogy(x0_test * ones(size(gamma)), sort(abs(gamma)), '-x', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'auto')
    hold on
end

if maxiter > 0
rho_rAA1 = NEPv_lambda_max(M_eigs)^0.25;
rho_PI   = max(abs(M_eigs));

% Tidy up each plot.
figure(1)
plot(k, rho_rAA1*ones(size(k)), '--b', 'LineWidth', 6)
plot(k, rho_PI*ones(size(k))  , ':k', 'LineWidth', 6)
xlabel('$k$')
ylabel('$\varrho_k(\mathbf{r}_0) := \Vert \mathbf{r}_k \Vert^{1/k}$')
%ylim([0.95*rho_rAA1_min, 1.1*rho_PI])
title(sprintf('test = %d', test_matrix))


% Set axis limits so plots are easier to read
if test_matrix == 1
    ylim([0.1 0.65])
end

if test_matrix == 2
    ylim([0.58 0.725])
end

if save_fig
    fig_name = sprintf('./figures/rhok-symmM-rand-test%d', test_matrix);
    figure_saver(gcf, fig_name, false);
end


% Eigenvalue plot
figure(2)
plot(eig(M), 'rx')
title('Eigenvalues $\{ m_i \}$ of $M$')
xlabel('$i$')
axis tight

if save_fig
    fig_name = sprintf('./figures/eig-symmM-rand-test%d', test_matrix);
    figure_saver(gcf, fig_name, false);
end

end


% Plot coefficients of residual in eigenvector basis.
figure(3)
title(sprintf('$k = %d$', maxiter-1)) % Note maxiter-1, since the AA code returns the residual for the second last iterate.
ylim([1e-10, 1])
xlabel('$\mathbf{r}_0$ index')
ylabel('$| V^{-1} \widehat{\mathbf{r}}_k (\mathbf{r}_0) |$')

if save_rk_decomposition_fig
    fig_name = sprintf('./figures/asymtotic-rk-symmM-rand-test%d-k%d', test_matrix, maxiter-1);
    figure_saver(gcf, fig_name, false);
end