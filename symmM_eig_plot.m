% Overlay on a single plot the three sets of eigenvalues saved in ./data
clc
clear all
close all


save_fig = ~true;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

load('./data/symmM_eigs_test1.mat')
i = linspace(0, 1, numel(M_eigs));
plot(i, sort(M_eigs), 'r>', 'DisplayName', '$\mathrm{test = 1}$')

hold on
load('./data/symmM_eigs_test2.mat')
i = linspace(0, 1, numel(M_eigs));
plot(i, sort(M_eigs), 'bx', 'DisplayName', '$\mathrm{test = 2}$')


hold on
load('./data/symmM_eigs_test3.mat')
i = linspace(0, 1, numel(M_eigs));
plot(i, sort(M_eigs), 'go', 'DisplayName', '$\mathrm{test = 3}$')

lh = legend();
lh.set('Location', 'Best', 'FontSize', 22)

axis tight
title('Eigenvalues $\{m_i\}_{i = 1}^n$ of $M$')
ax = gca;
ax.XTickLabel = {};
ax.XTick = [];
box on
xlabel('$k$', 'Color', 0.99*[1 1 1])
ylabel('$\varrho_k(\mathbf{r}_0) := \Vert \mathbf{r}_k \Vert^{1/k}$', 'Color', 0.99*[1 1 1])
% ax.YColor = 0.999*[1 1 1];
% ax.XColor = 0.999*[1 1 1];



if save_fig
    fig_name = sprintf('./figures/eig-symmM');
    figure_saver(gcf, fig_name, false);
end
