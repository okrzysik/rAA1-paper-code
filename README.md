This repository contains the MATLAB codes used for the numerical examples in the paper "Asymptotic convergence of restarted Anderson acceleration for certain normal linear systems" by
Oliver A. Krzysik, Hans De Sterck, and Adam Smith.

To ensure all of the examples work as intended, you should first run 'startup.m' to add subdirectories to your MATLAB path, and to configure plotting properties.

Each figure from the paper can be generated by running the scripts as listed below.

* Figure 1:
    Top row: 'symmM_two_by_two_r0_conv_cs.m' 
    Bottom row: 'symmM_rho_worst.m'

* Figure 2:
    Test 1 and test 2 convergence factors: 'skewM_rand_examples.m' 
    Test 3 convergence factor: 'symmM_two_grid.m'
    Eigenvalue plot: 'symmM_eig_plot.m'

* Figure 3:
    'symmM_nonlin.m'

* Figure 4:
    'skewM_conv_plot.m'

* Figure 5:
    'skewM_advection_test.m'
    Left panel: Uncomment the finite-difference example from the code
    Right panel: Uncomment the spectral example from the code

* Supplementary material Figure 1:
    'symmM_rand_examples_rAAm.m'

* Supplementary material Figure 2:
    'symmM_rand_examples.m': To generate the different plots in the figure, ensure test\_matrix = 1, and set the maxiter parameter to one more than the values indicated in the titles of the plots.

* Supplementary material Figure 3:
    'max_Xv_conjecture.m' 

* Supplementary material Figure 4:
    'skew_M_additional_examples.m' 



Additional information.
* 'utils/' has scripts implementing basic Picard iteration, the rAA(1) iteration, scripts relating to saving plots, and bits and pieces required by the multigrid examples

* 'data/' has .mat files storing the eigenvalues of the matrices M used in Figure 2

* 'figures/' has copies of all plots from the paper as .png files
