% Here we apply rAA(m) accelerate an underlying fixed-point iteration
% matrix using a symmetric matrix M. 
% In particular, we are interested in approximate periodicity of the
% resulting residuals.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear
close all
rng(1)

test_matrix = 1; % Choose me
%test_matrix = 2; % Or choose me

save_fig = ~true;

n = 14; % Dimension of matrix M
m_array = [1, 3, 5, 7]; % Restart parameter

num_x0  = 4; % Number of different x0 to try
maxiter = 300; % Number of iterations


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

% Starting iterate
X0 = zeros(n, num_x0);
for i = 1:num_x0
    X0(:, i) = rand(n, 1) - rand(n, 1);
end

% RHS vector
b = zeros(n, 1);

% Iteration function
q = @(x) M*x + b;

for m_idx = 1:numel(m_array)
        
    
    m = m_array(m_idx);

    % Loop over different starting iterate
    for x0_test = 1:num_x0
    
        x0 = X0(:, x0_test);
    
        [rnorms, R] = rAAm(q, x0, maxiter, m);
    
        % k_skip = 1:2*(m+1):maxiter;
        % 
        % cosangle_vec = [];
        % for idx = 1:numel(k_skip)-1
        %     u = R(:, k_skip(idx));   u = u / norm(u);
        %     v = R(:, k_skip(idx+1)); v = v / norm(v);
        %     cos_angle = u' * v;
        % 
        %     cosangle_vec = [cosangle_vec; cos_angle];
        % end
        % 
        % figure(2*m_idx)
        % semilogy(k_skip(1:end-1), abs(1-cosangle_vec), '-o')
        % hold on

        k_skip = 2*(m+1);
    
        cosangle_vec = [];
        for idx = 1:maxiter-k_skip
            u = R(:, idx);   u = u / norm(u);
            v = R(:, idx + k_skip); v = v / norm(v);
            cos_angle = u' * v;

            cosangle_vec = [cosangle_vec; cos_angle];
        end
    
        figure(2*m_idx)
        semilogy(1:maxiter-k_skip, abs(1-cosangle_vec), '-x')
        hold on

        % Just a sanity check that the method is converging.
        figure(2*m_idx+1)
        semilogy(1:maxiter, rnorms, '-x')
        hold on
    end
    
    % Tidy up each plot.
    figure(2*m_idx)
    xlabel('$k$')
    title(sprintf('$m = %d$', m))
    set(gca, 'xminorgrid','off','yminorgrid','off', 'xgrid','on','ygrid','on')
    if save_fig
        fig_name = sprintf('./figures/rAAm-periodic-m%d', m);
        figure_saver(gcf, fig_name, false);
    end

    % Tidy up each plot.
    figure(2*m_idx+1)
    xlabel('$k$')
    title(sprintf('residuals $m = %d$', m))
end