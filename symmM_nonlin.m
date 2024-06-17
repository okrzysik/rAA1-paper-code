% Solve the nonlinear elliptic BVP:
%   -(u_xx + u_yy) + gamma*u*exp(u) = f.
%
% This example comes from chapter 6 of the multigrid tutorial by Briggs et
% al. (see specifics from p. 102).
%
% We consider solving the problem with an inexact Newton method wherein the
% linearized system at each Newton iteration is approximately solved with a
% single two-grid iteration.
%
% We consider accelerating this inexact Newton iteration using rAA(1).
%
% There are two examples. 
%   1. f = 0. In this case we can iterate the method for as long as we want
%   in which case we can clearly see the improvement of rAA(1). On the
%   other hand, the solution in this case is trivial, u = 0. As such, when 
%   evaluated at the fixed point, the Jacobian is (almost) equal to the 
%   standard discretized Laplacian, which is kind of boring since its these
%   spectrum of this Jacobian that determines the asymptotic convergence
%   factor.
%
%   2. f = some prescribed function. In this case the solution u \neq 0,
%   so we run into numerical precision issues and can only iterate until we
%   get about 16 digits of accuracy in the residual. On the other hand, the
%   Jacobian evaluated at the fixed point is more interesting than just the
%   standard discrete Laplacian.
%   Note that in this case how tightly the convergence factor gets to the
%   theoretical predictions is a bit sensitive on exactly how u0 is chosen
%   and also how f is chosen.
%
%
% Further details:
%
% Problem is posed on unit square [0,1] \times [0,1] with zero Dirichlet
% BCs on all sides. There are nx-1 interior points/unknown per direction.
%
% The two-grid method is fairly standard using linear interpolation. It
% uses weighted Jacobi relaxation which seems to be fairly adequate. A
% Galerkin coarse-grid operator is used.
%
% As gamma increases so does the nonlinearity, but notice also that the
% problem becomes more and more diagonally dominant so it actually gets
% easier to solve. I tend to find that its hardest to solve for smaller
% values of gamma, e.g., ~1.
%
% NOTE: We solve a system of equations A(u) = 0. In Newton's method, A(u)
% is the residual, but this is not the residual we compute in AA, since in
% AA we write Newton's method as the fixed-point iteration u_new =
% q(u_old), and the AA residual is r_old = u_old - q(u_old) != A(u_old).
%
% Suppose T(u_k) is the two-grid error propagation matrix for the Jacobian 
% J_k of A(u_k), then Newton step takes the form u_{k+1} = u_k -
% T(u_k)*A(u_k), so that the iteration function is q(u) = u - T(u)*A(u). We
% need to evaluate the Jacobian of this q function and then compute its
% eigenvalues. 


tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear all
close all
rng(1)


save_fig = ~true;
compute_theory_bound = true;
num_different_u0 = 10; % Number of different initial iterates to test

% Problem parameters
maxiter = 40;
rtol    = 1e-16; % Halt when residual reduced by this amount
gamma   = 5;
nx      = 32;


n = (nx - 1)^2; % Total # of unknowns
h = 1/nx;
params.gamma = gamma;
params.nx    = nx;
params.h     = h;

x = (h:h:1-h)'; 
y = x;
[X, Y] = meshgrid(x, y);

% Pick a RHS f.
f = PDE_rhs(X, Y, params); f_is_zero = false;
f = rand(n, 1) - rand(n, 1); f_is_zero = false;
%f = zeros((nx-1)^2, 1);    f_is_zero = true;
params.f = f;

% Initial iterate for exact Newon method
u0 = zeros((nx-1)^2, 1);

% Interpolation
params.P = interpolate(nx);
% Restriction
params.R = 0.25*params.P';


%% Get exact solution so error can be computed
% If non-trivial RHS, compute exact solution with exact Newton method
if ~f_is_zero
    fprintf('Exact Newton\n')
    params.jacobian_solver = 'direct';
    rnorms = [];
    u_old = u0;
    for iter = 1:maxiter 
    
        [u_new, r_old] = newton_step(u_old, params);
        rnorms = [rnorms; norm(r_old)];
        if iter == 1
            fprintf('||r||_%d = %.2e\n', iter-1, rnorms(iter))
        else
            fprintf('||r||_%d = %.2e, ||r||_%d / ||r||_%d = %.2e\n', ...
                iter-1, rnorms(end), iter, iter-1, rnorms(end)/rnorms(end-1))
        end
        
        if rnorms(end)/rnorms(1) < rtol
            break
        end
        u_old = u_new;
    end

    u_exact = u_new;

    % figure(1)
    % U_EXACT = reshape(u_exact, [nx-1, nx-1]);
    % mesh(X, Y, U_EXACT);
    % 
    % figure(2)
    % U_EXACT = reshape(PDE_sol(X, Y, params), [nx-1, nx-1]);
    % mesh(X, Y, U_EXACT);
    % pause

% If f is zero then exact solution is too.
else
    u_exact = zeros(size(f));

end

% Picard and rAA(1) iterations. These both use an inexact Newton iteration,
% so ensure jacobian_solver is a two-grid method and not an exact solve.
params.jacobian_solver = 'two-grid';

% Iteration function
q = @(u) newton_step(u, params);

% Record max number of iterations taken over all starting iterates, for plotting.
maxiter_actual = 0;

% Loop over different u0
for u0_idx = 1:num_different_u0

    % Initial iterate
    u0 = 1e-1*(rand(n, 1) - rand(n, 1));

    
    %% Picard iteration
    rnorms_PI = PI(q, u0, maxiter, rtol);
    
    numiter_actual = numel(rnorms_PI);
    maxiter_actual = max(maxiter_actual, numiter_actual);
    
    figure(1)
    % Note the indexing. We don't raise ||r_0|| to the power 1/0, but
    % instead start ||r_k||^{1/k} from k=1, noting that rnorms(1) is
    % ||r_0||
    plot(rnorms_PI(2:end).^(1./(1:numiter_actual-1)'), '-<')
    hold on
    ylim([0, 1])
    
    figure(2)
    semilogy(0:numiter_actual-1, rnorms_PI/rnorms_PI(1), '-<')
    hold on


    %% rAA(1) iteration
    rnorms_AA = rAA1(q, u0, maxiter, rtol);
    
    numiter_actual = numel(rnorms_AA);
    maxiter_actual = max(maxiter_actual, numiter_actual);

    figure(1)
    plot(rnorms_AA(2:end).^(1./(1:numiter_actual-1)'), '-*')
    hold on
    ylim([0, 1])
    
    figure(2)
    semilogy(0:numiter_actual-1, rnorms_AA/rnorms_AA(1), '-*')
    hold on

end


%% Pretty up and save figures
figure(1)
ylabel('$\Vert \mathbf{r}_k \Vert^{1/k}$')
xlabel('$k$')
xlim([1, maxiter_actual-1])
ylim([0.1 0.4])

figure(2)
ylabel('$\Vert \mathbf{r}_k \Vert / \Vert \mathbf{r}_0 \Vert$')
xlabel('$k$')
axis tight
xlim([0, maxiter_actual])


if save_fig
    fig_name = sprintf('./figures/rk-newton-nx%d', nx);
    figure_saver(gcf, fig_name, false);
end



%% Compute theoretical values for convergence factors.
if compute_theory_bound
% We need to evaluate the Jacobian of q at the fixed point, i.e., the exact
% solution u_exact. To approximate this Jacobian we use a finite difference and
% thus only require the action of q. I.e., the directional derivative of q
% at the point u in the direction d is approximately (q(u + mu*d) - q(u))/mu
% for small parameters mu. To get the jth column of this Jacobian we take d
% as the jth canonical unit vector.

I       = eye((nx-1)^2); % Columns of I are the canonical unit vectors we 
% need to compute the derivative of q in the direction of.
q_prime = zeros((nx-1)^2, (nx-1)^2);

mu = 1e-6; % Finite difference step size.

%params.jacobian_solver = 'direct';

q_u = u_exact; % q at the fixed point is just u_exact.
for col = 1:(nx-1)^2
    e = I(:, col);

    % Evaluate q(u + mu*e)
    [q_umue, ~]     = newton_step(u_exact + mu*e, params);

    % Insert column into q_prime.
    q_prime(:, col) =  (q_umue - q_u)/mu;
end

M = q_prime;

% Jacobian is not symmetric, but does it have real eigenvalues?
m = eig(M);

assert(max(abs(imag(m))) < 1e-10, 'M eigs have significant imaginary component...')
m = real(m);

rho_rAA1_theory = NEPv_lambda_max(m)^0.25;
rho_PI_theory   = max(abs(m));

% Overlay theoretical convergence factors
figure(1)
hold on
k = 1:maxiter_actual;
plot(k, rho_rAA1_theory*ones(size(k)), '--b', 'LineWidth', 6)
plot(k, rho_PI_theory*ones(size(k))  , ':k', 'LineWidth', 6)

end


if save_fig
    fig_name = sprintf('./figures/rhok-newton-nx%d', nx);
    figure_saver(gcf, fig_name, false);
end


% One step of Newton's method. Given an approximation u_old, update it to
% u_new. r_old is the residual for the input.
function [u_new, r_old] = newton_step(u_old, params)
    r_old    = residual(u_old, params);

    e_lin    = zeros(size(r_old));
    J        = jacobian(u_old, params);

    if strcmp(params.jacobian_solver, 'direct')
        e_lin = J \ -r_old;

    elseif strcmp(params.jacobian_solver, 'two-grid')
        
        J_coarse = params.R * J * params.P;
        e_lin = zeros(size(r_old)); % Initial guess at e_lin
        %two_grid(x, f, A, A_coarse, P, R, options);
        options.omega = 4/5; % Jacobi weight
        e_lin = two_grid(e_lin, -r_old, J, J_coarse, params.P, params.R, options);
    end

    u_new    = u_old + e_lin;
end


% Assembles the Jacobian J at the current point u. This is based in the
% description on p. 104 of Briggs et al.
% J is block tridiagonal with diagonal blocks B on its block sub and super
% diagonals and tridigonal blocks J_i, i = 1,...,nx-1 along its block
% diagonal.
%
% MATLAB's sparse command is used to assemble the matrix, and this just
% requires (row,col,data) triplets for every non-zero.
function J = jacobian(u, params)
    % Unpack parameters
    nx    = params.nx;
    h     = params.h;
    gamma = params.gamma;

    u     = reshape(u, [nx-1, nx-1]);

    % subdiagonal has nx-2 B blocks, so does super diagonal. B is just
    % diagonal, so has nx-1 nnz.
    % main diagonal has nx-1 blocks, each tridiagonal.
    nnz = (nx-1)*(nx-2) + ((nx-2) + (nx-1) + (nx-2))*(nx - 1) + (nx-1)*(nx-2);

    ii = [];
    jj = [];
    vv = [];

    % Local non-zero data for B matrices
    B_ii = (1:nx-1)';
    B_jj = (1:nx-1)';
    B_vv = -1/h^2 * ones(nx-1, 1);

    % Local non-zero data for J_i matrices.
    % Each J_i block is tridiagonal. Organize the data so that the diagonal
    % data for J_i is stored in the first nx-1 places so that it's simple to
    % access and ammend in the loop below.
    J_ii = [(1:nx-1)';      (2:nx-1)';              (1:nx-2)'];
    J_jj = [(1:nx-1)';      (1:nx-2)';              (2:nx-1)'];
    % To begin with, insert dummy data into the first nx-1 entries of J_vv
    J_vv = [zeros(nx-1, 1); -1/h^2 * ones(nx-2, 1); -1/h^2 * ones(nx-2, 1)];

    % Loop over all block rows and in each, add the row col and data
    % entries for the blocks into those arrays for J.
    for block_row_idx = 1:nx-1

        % Fix up diagonal J data.
        J_vv(1:nx-1) = 4/h^2 + gamma*( 1 + u(block_row_idx, :) ).*exp( u(block_row_idx, :) );

        % One less than index of first row in current block row
        global_row_offset = (block_row_idx-1)*(nx-1);
        
        % First block row: No B on block sub diagonal
        if block_row_idx == 1

            ii = [ii; J_ii + global_row_offset; B_ii + global_row_offset];
            jj = [jj; J_jj + global_row_offset; B_jj + global_row_offset + (nx-1)];
            vv = [vv; J_vv; B_vv];

        % Last block row: No B on block super diagonal
        elseif block_row_idx == nx-1

            ii = [ii; B_ii + global_row_offset;          J_ii + global_row_offset];
            jj = [jj; B_jj + global_row_offset - (nx-1); J_jj + global_row_offset];
            vv = [vv; B_vv; J_vv];

        % General case where there are three blocks per row.
        else

            ii = [ii; B_ii + global_row_offset;          J_ii + global_row_offset; B_ii + global_row_offset];
            jj = [jj; B_jj + global_row_offset - (nx-1); J_jj + global_row_offset; B_jj + global_row_offset + (nx-1)];
            vv = [vv; B_vv; J_vv; B_vv];

        end

    end

    J = sparse(ii, jj, vv, (nx-1)^2, (nx-1)^2, nnz);
end

% Nonlinear residual corresponding to discretized PDE. 
function r = residual(u, params)
    
    % Unpack parameters
    nx    = params.nx;
    h     = params.h;
    gamma = params.gamma;
    f     = params.f;

    % Reshape vectors into 2D arrays so indexing is easier.
    u = reshape(u, [nx-1, nx-1]);
    f = reshape(f, [nx-1, nx-1]);
    r = zeros(size(u));


    % Compute the residual in each of the cases. Note all zero Dirichlet
    % BCs are used. We have special cases along each boundary and at each
    % corner since in these cases at least one of the points we connect to
    % is on a boundary.
    
    % Bottom-left corner
    ii = 1; jj = 1;
    r(ii, jj) = (4*u(ii, jj) ...
                - 0          - u(ii+1, jj) - 0           - u(ii, jj+1) ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Bottom-right corner
    ii = nx-1; jj = 1;
    r(ii, jj) = (4*u(ii, jj) ...
                - u(ii-1, jj) - 0          - 0            - u(ii, jj+1) ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Top-left corner
    ii = 1; jj = nx-1;
    r(ii, jj) = (4*u(ii, jj) ...
                - 0           - u(ii+1, jj) - u(ii, jj-1) - 0 ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Top-right corner
    ii = nx-1; jj = nx-1;
    r(ii, jj) = (4*u(ii, jj) ...
                - u(ii-1, jj) - 0           - u(ii, jj-1) - 0 ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Bottom boundary
    ii = (2:nx-2)'; jj = 1;
    r(ii, jj) = (4*u(ii, jj) ...
                - u(ii-1, jj) - u(ii+1, jj) - 0          - u(ii, jj+1) ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Left boundary
    ii = 1; jj = (2:nx-2)';
    r(ii, jj) = (4*u(ii, jj) ...
                - 0          - u(ii+1, jj) - u(ii, jj-1) - u(ii, jj+1) ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Right boundary
    ii = nx-1; jj = (2:nx-2)';
    r(ii, jj) = (4*u(ii, jj) ...
                - u(ii-1, jj) - 0          - u(ii, jj-1) - u(ii, jj+1) ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));

    % Top boundary
    ii = (2:nx-2)'; jj = nx-1;
    r(ii, jj) = (4*u(ii, jj) ...
                - u(ii-1, jj) - u(ii+1, jj) - u(ii, jj-1) - 0 ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));
    
    % Standard case: No boundary data required.
    ii = (2:nx-2)'; jj = (2:nx-2)';
    r(ii, jj) = (4*u(ii, jj) ...
                - u(ii-1, jj) - u(ii+1, jj) - u(ii, jj-1) - u(ii, jj+1) ...
                )/h^2 ...
                + gamma*u(ii, jj).*exp(u(ii, jj));
            
    r = r - f;
    % Flatten back into vector.
    r = r(:);
end


% RHS of PDE given on p. 103 of Briggs et al.
function f = PDE_rhs(X, Y, params)
    gamma = params.gamma;
    f = 2*( (X - X.^2) + (Y - Y.^2) ) + gamma*( X - X.^2 ).*(Y - Y.^2).*exp( (X - X.^2).*(Y - Y.^2) );
    f = f(:);
end

% Exact PDE solution for the above RHS.
function u = PDE_sol(X, Y, params)
    u = (X - X.^2) .* (Y - Y.^2);
    u = u(:);
end