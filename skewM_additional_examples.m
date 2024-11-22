% For skew-symmetric matrices M, compare theoretical rAA(1) and PI 
% convergence factor to numerically computed convergence factors.
%
% Solve the linear system A*x = b, where A = I - M with M skew symmetric.
%
% In this test a sequence of randomly generated matrices M are considered.
% In all cases, it seems that the asymptotic convergence factor of rAA(1)
% is independent of the initial iterate.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear
close all
rng(1)

save_fig = ~true;

% Try several different size M.
n_array = [10, 14, 22]; % Dimension of M.
eig_markers = {'r>', 'bx', 'go'}; % Markers for plotting corresponding eigenvalues

num_x0  = 30; % Number of different x0 to try
maxiter = 300; % Number of iterations

% Iteration numbers for plotting.
k = (0:maxiter-1)';

plot_skip = 5;
k_plot = k(2:plot_skip:end);


for n_idx = 1:numel(n_array)
    n = n_array(n_idx);

    % Take M as the skew-symmetric part of a randomly generated matrix.
    M = rand(n, n);
    M = M' - M;
    a = 2/3;
    b = 2;
    scale = a + (b-a)*rand(1); % generate random # between (a,b)
    M = M / (scale*norm(M)); % Force M to have norm in interval (1/b, 1/a)

    m_max = max(abs(eig(full(M))));

    % RHS vector
    b = zeros(n, 1); % This should be zero to avoid round off errors
    
    % Iteration function
    q = @(x) M*x + b;
    
    plot_skip = 5;
    k_plot = k(1:plot_skip:end);
    
    for x0_test = 1:num_x0
        % Starting iterate
        x0 = rand(n, 1) - rand(n, 1);

        rnorms_rAA = rAA1(q, x0, maxiter);
        rnorms_PI  = PI(q, x0, maxiter);
            
        % Root-average conv factor
        figure(n_idx)
        plot(k_plot, rnorms_rAA(2:plot_skip:end).^(1./k_plot), '-*', ...
            'MarkerSize', 8)
        hold on
        plot(k_plot, rnorms_PI(2:plot_skip:end).^(1./k_plot),  '-<', ...
            'MarkerSize', 8)
        
        % Residual norms
        figure(2)
        semilogy(k, rnorms_rAA, '-*')
        hold on
        semilogy(k, rnorms_PI,  '-<')
        
    end

    rho_rAA1 = rAA1_conv_fac(m_max);
    rho_PI   = m_max;
    
    % Tidy up each plot.
    figure(n_idx)
    plot(k, rho_rAA1*ones(size(k)), '--b', 'LineWidth', 6)
    plot(k, rho_PI*ones(size(k))  , ':k', 'LineWidth', 6)
    xlabel('$k$')
    ylabel('$\varrho_k(\mathbf{r}_0) := \Vert \mathbf{r}_k \Vert^{1/k}$')
    title(sprintf('test = %d', n_idx))
    ylim([0.9*rho_rAA1, 1.1*rho_PI])
    xlim([min(k), max(k)])

    if save_fig
        fig_name = sprintf('./figures/skew-M-add-test%d-rhok', n_idx);
        figure_saver(gcf, fig_name, false);
    end

    % Plot eigenvalues
    figure(101)
    m_abs = abs(eig(M));
    i = linspace(0, 1, numel(m_abs));
    plot(i, sort(m_abs), eig_markers{n_idx}, 'DisplayName', sprintf('test = %d', n_idx))
    hold on
end

% Tidy up eigenvalue plot.
axis tight
title('Absolute eigenvalues $\{|m_i|\}_{i = 1}^n$ of $M$')
lh = legend();
lh.set('location', 'best')
ax = gca;
ax.XTickLabel = {};
ax.XTick = [];
box on
xlabel('$k$', 'Color', 0.99*[1 1 1])
ylabel('$\varrho_k(\mathbf{r}_0) := \Vert \mathbf{r}_k \Vert^{1/k}$', 'Color', 0.99*[1 1 1])

if save_fig
    fig_name = sprintf('./figures/skew-M-add-test-eigs');
    figure_saver(gcf, fig_name, false);
end

% Worst-case conv fact for rAA(1).
function rho = rAA1_conv_fac(m)
    rho = m ./ (( 1 + m.^2 ).^0.25);
end