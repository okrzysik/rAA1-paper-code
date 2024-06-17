% For skew symmetric M, plot the convergence factor of rAA(1) vs. Picard
% iteration. Also plot the number of iterations required to reach a given
% tolerance of 10^{-nu} given each convergence factor.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear all
close all

%%%%%%%%%%%%%%%%%
save_fig = true;
%%%%%%%%%%%%%%%%%


m_max_PI = linspace(0, 1);
m_max_AA = linspace(0, 1.272);

% Worst-case root convergence factors for each method
rho_worst_PI = m_max_PI;
rho_worst_AA = m_max_AA ./ (1 + m_max_AA.^2).^0.25;

mycols    = {'k', 'r', 'b', [0, 0.5, 0], [0.75, 0, 0.95], [0.9290, 0.6940, 0.1250], [0.3010 0.7450 0.9330]};


%% Convergence factor plot
figure(1)
plot(m_max_AA, rho_worst_AA, '--*', ...
    'LineWidth', 1, ...
    'Color', mycols{5}, ...
    'DisplayName', 'rAA(1)')
hold on
plot(m_max_PI, rho_worst_PI, '-<', ...
    'LineWidth', 1, ...
    'Color', mycols{6}, ...
    'DisplayName', 'PI')
xlabel('$\rho(M)$')
ylabel('$\varrho^{\mathrm{worst}}$')
axis tight
lh = legend();
lh.set('Location', 'Best')


if save_fig
    fig_name = sprintf('./figures/skewM_rho_worst');
    figure_saver(gcf, fig_name, false);
end


%% Number of iterations plot
figure(2)
nu_array = (1:2:7);
for nu_idx = 1:numel(nu_array)
    
    nu = nu_array(nu_idx);

    plot(m_max_AA, -nu ./ log10(rho_worst_AA), ...
        '--*', ...
        'color', mycols{nu_idx}, ...
        'DisplayName', sprintf('$\\nu_* = %d$', nu))
    hold on
    plot(m_max_PI, -nu ./ log10(rho_worst_PI), ...
        '-<', ...
        'color', mycols{nu_idx}, ...
        'HandleVisibility', 'off')
end

axis tight
ylim([0, 50])
xlim([0, m_max_AA(end)])

lh = legend();
lh.set('Location', 'Best')
xlabel('$\rho(M)$')
ylabel('$k_*$')
%title('$k_* = -\nu_* / \log_{10}(\varrho^{\mathrm{worst}}(M))$')
title('$\Vert \mathbf{e}_{k_*} \Vert / \Vert \mathbf{e}_0 \Vert = 10^{-\nu_*}$')

if save_fig
    fig_name = sprintf('./figures/skewM_num_iters');
    figure_saver(gcf, fig_name, false);
end