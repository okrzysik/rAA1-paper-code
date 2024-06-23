% For skew-symmetric matrices M, compare theoretical rAA(1) and PI 
% convergence factor to numerically computed convergence factors.
%
% Solve the linear system A*x = b, where A = I - M with M skew symmetric.
% We consider linear systems of this form that arise from the 
% discretization of the linear advection problem u_t = (alpha*u)_x. If 
% (alpha*u)_x is discretized with a skew-symmetric discretization, say, D, 
% then the PDE becomes an ODE system du/dt = L*u, where L*u \approx 
% (alpha*u)_x. Then, applying any DIRK time integration or Crank-Nicolson 
% results in linear systems of the form (I - a*dt*L)*x = b, where dt is the 
% time-step size and a is some constant.
% 
% We can then solve this linear system by the Picard iteration
% x_{k+1} = (a*dt*L)*x_{k} + b.
% So, we identify the fixed-point iteration matrix M as M = a*dt*L. For L 
% we consider 
%   1. Central finite differences of 2nd and 6th order
%   2. Fourier collocation spectral method
%
% In summary, we have something like:
% PDE: u_t = (alpha*u)_x  -> ODE:  du/dt = L*u  ->  Linear system:  (I - a*dt*L)*u^{n+1} = b^{n}. 
%
% For info on the spectral method, see Hesthaven (2017) section 13.2.2

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 
close all
rng(1)


save_fig = ~true;


%% Choose your spatial discretization: Uncomment either block of code
% Fourier collocation method
disc = 'spectral';
n = 31; % Dimension of linear system
CFL = 0.85;

% % Central finite difference
% disc = 'finite-difference';
% FD_order = 2; % Order of accuracy of finite difference
% % FD_order = 6;
% n = 64; % Dimension of linear system
% CFL = 1.7;


num_r0 = 30; % Number of different r0 to try
maxiter = 300; % Number of iterations


% Solve u_t = (alpha*u)_x. Choose some alpha > 0 such that diag(alpha) is
% SPD such that its square root exists.
alpha = @(x) 0.5*(1 + cos(x).^2);
%alpha = @(x) ones(size(x));

% Equispaced grid on [0, 2pi). 
if strcmp(disc, 'finite-difference')
    x = linspace(0, 2*pi, n+1)';
    x(end) = [];
    h = x(2) - x(1); % Grid spacing.
    D = spatial_disc(n, h, FD_order);
    
    
% Fourier collocation method. This means we force the residual to vanish
% pointwise at 2*N + 1 collocation points, resulting in 2*N + 1 eqations in 
% 2*N + 1 unknowns. If u(x, t) is the exact PDE solution, then
%   u(x, t) approx u_h(x, t) = sum_{j = 0}^{2N} u_h(x_j, t) * h_j(x), 
% where h_j(x) is the Lagrange polynomials defined by x_j. 
% The unknowns we solve for are u_h(x,j t) at all points x_j. Note that the
% collocation points are equispaced and equal to the interpolation points.
elseif strcmp(disc, 'spectral')

    assert(mod(n, 2) == 1, 'n must be odd because n = 2*N + 1')
    % n = 2*N + 1 -> N = (n - 1)/2
    N = (n - 1)/2;
    x = 2*pi/(2*N+1)*(0:2*N)';
    h = x(2) - x(1); % Grid spacing.
    D = FourierD(N, 1);

end

%% Get spatial discretization matrix L.
dt = CFL*h; % Time-step size.

% Diagonal matrix of wave-speed evaluated at FD mesh points or interpolation points.
Lambda = diag(alpha(x));
L      = D*Lambda;


%% Sanity check: Time-step with CN and check the solution looks right in the eyeball norm
t = (0:dt:2*pi)';
U = zeros(numel(x), numel(t));
U(:, 1) = cos(x); % Initial condition

% Crank-Nicolson applied to u'(t) = L*u does 
% (u_{n+1} - u_n)/dt = 0.5*( L*u_{n+1} + L*u_n ), so that
% [I - dt*0.5*L]*u_{n+1} = [I + dt*0.5*L]*u_n
% or Y*u_{n+1} = Z*u_{n}.
Y = eye(n) - 0.5*dt*L;
Z = eye(n) + 0.5*dt*L;

for step = 2:numel(t)
    U(:, step) = Y \ (Z*U(:, step-1));
end

figure(11)
[X, T] = meshgrid(x, t);
mesh(X, T, U')
xlabel('$x$')
ylabel('$t$')
zlabel('$u(x, t)$')
title('Sanity check: PDE solution')


%% Get matrix M
M     = 0.5*dt*L;
M_eigs = eig(full(M));
assert(max(abs(real(M_eigs))) < 1e-10, 'M eigs are not imaginary...')
m_max = max(abs(M_eigs)); % Spectral radius

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
xlim([min(k), max(k)])
if strcmp(disc, 'spectral')
    title(sprintf('$\\textrm{spectral},\\, \\rho(M)=%.2f$', m_max))
elseif strcmp(disc, 'finite-difference')
    title(sprintf('$\\textrm{finite-difference},\\, \\rho(M)=%.2f$', m_max))
end
box on

if save_fig
    if strcmp(disc, 'spectral')
        fig_name = sprintf('./figures/rhok-skewM-spec-n%d-N%d', n, N);
    elseif strcmp(disc, 'finite-difference')
        fig_name = sprintf('./figures/rhok-skewM-FD-n%d-L%d', n, FD_order);
    end
    figure_saver(gcf, fig_name, false);
end

figure(2)
xlabel('$k$')
ylabel('$\Vert \mathbf{r}_k \Vert$')
xlim([min(k), max(k)])


% Worst-case conv fact for rAA(1).
function rho = rAA1_conv_fac(m)
    rho = m ./ (( 1 + m.^2 ).^0.25);
end


% Skew symmetric spatial disc of d/dx with n points spaced by h.
function D = spatial_disc(n, h, order)
    e = ones(n, 1);

    if order == 2
        D = spdiags([-e 0*e e], -1:1, n, n) / (2*h);
        % Set up periodic BCs.
        D(1, end) = -D(1, 2);
        D(end, 1) = -D(end, end-1);

    elseif order == 6
        assert(n >= 9, 'Disc implementation does not work for fewer than 9 points')

        D = spdiags([-1*e 9*e -45*e 0*e 45*e -9*e 1*e], -3:3, n, n) / (60*h);
        % Set up periodic BCs.
        % Top RH corner
        E = full(D(4:6, 1:3));
        D(1:3, end-2:end) = E;

        % Bottom LH corner
        E = full(D(end-5:end-3, end-2:end));
        D(end-2:end, 1:3) = E;

    else

        error('Only 2nd- and 6th-order spatial discretizations available')
    end
end


% This function is from the MATLAB code in Hesthaven's book.
function [D] = FourierD(N,q)
    % function [D] = FourierD(N);
    % Purpose: Initialize q'th order Fourier differentiation matrix
    column = [0 (-1).^(1:2*N)./(2*sin(pi/(2*N+1)*(1:2*N)))]';
    D = toeplitz(column, column([1 (2*N+1):-1:2]))^q;
end

