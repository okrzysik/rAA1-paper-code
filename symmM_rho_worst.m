% plots the worst-case root-linear convergence factor for rAA(1) on the
% space of all 2x2 symmetric matrices M and plot this relative to the
% convergence factor of the underlying Picard iteration.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear 

save_figs = ~true;


% m_min = -3;
% m_max =  3;

m_min = -1;
m_max =  1;

nu = 1e-3;
m1 = linspace(m_min+nu, m_max-nu, 1e3);
m2 = m1;

[M1, M2] = meshgrid(m1, m2);



%% AA(1)
rho_worst_AA = abs( M1 .* M2 .* (M2 - M1) ) ./ ( abs(M1.*(M1 - 1)) + abs(M2.*(M2 - 1)) ); 
rho_worst_AA = sqrt(rho_worst_AA);

% Mask values larger than 1 to NaN so they don't appear in the plot.
rho_worst_AA(rho_worst_AA > 1) = NaN;

figure
contourf(M1, M2, rho_worst_AA, 10)
xlabel('$m_1$')
ylabel('$m_2$')
title('$\varrho^{\mathrm{worst}}_{\mathrm{AA}}$')
caxis manual
caxis([0 1]);
colorbar

yticks([m_min:(m_max - m_min)/4:m_max])
yticks([m_min:(m_max - m_min)/4:m_max])
axis([m_min, m_max, m_min, m_max])

if save_figs
    figure_saver(gcf, './figures/rho-worst-AA', ~true)
end


%% Picard iteration
rho_worst_PI = max(abs(M1), abs(M2));

% Mask values larger than 1 to NaN so they don't appear in the plot.
rho_worst_PI(rho_worst_PI > 1) = NaN;

figure
contourf(M1, M2, rho_worst_PI, 10)
xlabel('$m_1$')
ylabel('$m_2$')
title('$\varrho^{\mathrm{worst}}_{\mathrm{PI}}$')
caxis manual
caxis([0 1]);
colorbar

yticks([m_min:(m_max - m_min)/4:m_max])
yticks([m_min:(m_max - m_min)/4:m_max])
axis([m_min, m_max, m_min, m_max])


if save_figs
    figure_saver(gcf, './figures/rho-worst-PI', ~true)
end


%% Their ratio
figure
contourf(M1, M2, rho_worst_AA ./ rho_worst_PI, 10)
hold on
plot(m1, -m2, 'k--', 'LineWidth', 3)
xlabel('$m_1$')
ylabel('$m_2$')
title('$\varrho^{\mathrm{worst}}_{\mathrm{AA}} / \varrho^{\mathrm{worst}}_{\mathrm{PI}}$')
caxis manual
caxis([0 1]);
colorbar

yticks([m_min:(m_max - m_min)/4:m_max])
yticks([m_min:(m_max - m_min)/4:m_max])
axis([m_min, m_max, m_min, m_max])

if save_figs
    figure_saver(gcf, './figures/rho-worst-ratio', ~true)
end