% For skew-symmetric matrices M, compare theoretical rAA(1) and PI 
% convergence factor to numerically computed convergence factors.
%
% Solve the linear system A*x = b, where A = I - M with M skew symmetric.
% We consider linear systems of this form that arise from the 
% discretization of the linear advection problem u_t = u_x. If u_x is
% discretized with a central discretization, say, L, then the PDE becomes
% an ODE system du/dt = L*u, where L*u \approx u_x. Then, applying any DIRK
% time integration results in linear systems of the form (I - dt*L)*x = b,
% where dt is the time-step size. If the DIRK method is backward Euler,
% then the linear system is (I - dt*L)*u^{n+1} = u^{n}, with u^{n} an
% approximation to u at time t_n, but more generally, the quantity to be
% solved for could be any stage vector in a DIRK method.
% 
% So, we identify the matrix M as M = dt*L. For L we consider both
% 2nd- and 6th-order spatial discretizations. Note that for skew-symmetric
% spatial discretizations one often has to use implicit time integration
% for stability since explicit methods aren't well suited to handle
% eigenvalues along the imaginary axis.
%
% In summary, we have something like:
% PDE: u_t = u_x  -> ODE:  du/dt = L*u  ->  Linear system:  (I - dt*L)*u^{n+1} = u^{n}. 

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 
close all
rng(1)


save_fig = ~true;

% Order of accuracy of the spatial discretization L
space_disc_order = 2; n = 90; % Dimension of system/number of spatial DOFs
%space_disc_order = 6; n = 30; % Dimension of system/number of spatial DOFs


num_r0 = 30; % Number of different r0 to try
maxiter = 200; % Number of iterations


% Get matrix M
CFL   = 0.7; % dt = h*CFL.
h     = 1/n;
dt    = CFL * h;
L     = spatial_disc(n, space_disc_order); %full(L)
M     = dt*L;
m_max = max(abs(eig(full(M))));


% RHS vector
b = zeros(n, 1);

% Iteration function
q = @(x) M*x + b;

% Iteration numbers for plotting.
k = (0:maxiter-1)';

plot_skip = 5;
k_plot = k(2:plot_skip:end);

for r0_test = 1:num_r0    
    % Starting iterate
    x0 = rand(n, 1) - rand(n, 1);

    rnorms_rAA = rAA1(q, x0, maxiter);
    rnorms_PI  = PI(q, x0, maxiter);
   
    % Root-average conv factor
    figure(1)
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
    
    % % Ratio of relative residual norms
    % figure(3)
    % plot(k(1:(end-1)/2), rnorms_rAA(3:2:end)./rnorms_rAA(1:2:end-2), '-*')
    % hold on
    % plot(k(1:(end-1)/2), rnorms_PI(3:2:end)./rnorms_PI(1:2:end-2),  '-<')
end

rho_rAA1 = rAA1_conv_fac(m_max);
rho_PI   = m_max;

% Tidy up each plot.
figure(1)
plot(k, rho_rAA1*ones(size(k)), '--b', 'LineWidth', 6)
plot(k, rho_PI*ones(size(k))  , ':k', 'LineWidth', 6)
xlabel('$k$')
ylabel('$\varrho_k(\mathbf{r}_0) := \Vert \mathbf{r}_k \Vert^{1/k}$')
ylim([0.9*rho_rAA1, 1.1*rho_PI])
%xlim([min(k), max(k)])
title(sprintf('$M = \\delta t L_{%d},\\, |m_{\\max}|=%.2f$', space_disc_order, m_max))

if save_fig
    fig_name = sprintf('./figures/rhok-skewM-n%d-L%d', n, space_disc_order);
    figure_saver(gcf, fig_name, false);
end


figure(2)
xlabel('$k$')
ylabel('$\Vert \mathbf{r}_k \Vert$')
xlim([min(k), max(k)])

% figure(3)
% plot(k(1:(end-1)/2), rho_rAA1^2*ones(size(k(1:(end-1)/2))), '--k', 'LineWidth', 4)
% plot(k(1:(end-1)/2), rho_PI^2*ones(size(k(1:(end-1)/2)))  , ':b', 'LineWidth', 4)
% xlabel('$k$')
% ylabel('$\Vert \mathbf{r}_{2k+2} \Vert / \Vert \mathbf{r}_{2k} \Vert$')


% Worst-case conv fact for rAA(1).
function rho = rAA1_conv_fac(m)
    rho = m ./ (( 1 + m.^2 ).^0.25);
end


% Skew symmetric spatial disc of d/dx with n points. 
function L = spatial_disc(n, order)
    e = ones(n, 1);
    h = 1/n;

    if order == 2
        L = spdiags([-e 0*e e], -1:1, n, n) / (2*h);
        % Set up periodic BCs.
        L(1, end) = -L(1, 2);
        L(end, 1) = -L(end, end-1);

    elseif order == 6
        assert(n >= 9, 'Disc implementation does not work for fewer than 9 points')

        L = spdiags([-1*e 9*e -45*e 0*e 45*e -9*e 1*e], -3:3, n, n) / (60*h);
        % Set up periodic BCs.
        % Top RH corner
        E = full(L(4:6, 1:3));
        L(1:3, end-2:end) = E;

        % Bottom LH corner
        E = full(L(end-5:end-3, end-2:end));
        L(end-2:end, 1:3) = E;

    else

        error('Only 2nd- and 6th-order spatial discretizations available')
    end
end

