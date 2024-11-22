% plots cross-sections of the r0-dependent r-linear convergence factor 
% for 2x2 symmetric matrices M with eigenvalues m1, and m2. The convergence 
% factor is plotted as a function of epsilon which parametrizes the initial 
% residual, r0 = v1 + epsilon v2, and cross sections are shown for several 
% values of m1, m2.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clc
clear

save_fig = ~true;
fig_name = './figures/two_by_two_r0_conv_cs';

rng(4)

num_examples = 9;

epsilon = logspace(-2, 2, 50);

styles = {'-o',':x', '-.+', '--*', '-s', ':d', '-.v', '-->'};
ax = axes(); 
ax.LineStyleOrder = styles; 

for kk = 1:num_examples
    m1 = round(randn(1), 1, 'significant');
    m2 = round(randn(1), 1, 'significant');

    while m1 == 1
        m1 = round(randn(1), 1, 'significant');
    end

    while m2 == 1
        m2 = round(randn(1), 1, 'significant');
    end

    semilogx(epsilon, conv_fac(m1, m2, epsilon), ...
        styles{mod(kk, length(styles))+1}, ...
        'DisplayName', sprintf('$(%.1f,%.1f)$', m1, m2), ...
        'LineWidth', 2)
        
    hold on
end

box on
lh = legend();
lh.set('Location', 'NorthOutside')
%lh.set('Location', 'EastOutside')
%lh.NumColumns = num_examples/3;
lh.NumColumns = num_examples;
lh.Orientation = 'horizontal';
xlabel('$\varepsilon$')
ylabel('$\varrho_{\mathrm{AA}}(\mathbf{r}_0)$')
xlim([-0.25, epsilon(end)])

function rho = conv_fac(m1, m2, eps)
    rho = (m2 - m1).^2 .* ...
          (m1.*m2).^2 ./ ...
          ( (m1 - 1).^2 + eps.^2*(m2 - 1).^2 ) ./ ... 
          ( m1.^2 + m2.^2 ./ eps.^2 );

    rho = rho.^0.25;
end

if save_fig
    figure_saver(gcf, fig_name, false);
end