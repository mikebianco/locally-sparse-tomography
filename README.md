# Locally sparse tomography
This is a Matlab implementation of the locally sparse travel time tomography (LST) algorithm for seismic and acoustic data as presented in the paper:
M.J. Bianco and P. Gerstoft, "Travel time tomography with adaptive dictionaries," IEEE Trans. Computational Imaging, Vol. 4, No. 4, 2018.

The other tomography methods in the paper, conventional and total variation-regularized tomography, are also implemented for comparison. Figures similar to those in the paper (Bianco and Gerstoft 2018) are generated in a series of synthetic tests. These tests include cases with and without noise.

The simulations are run from 'main_LST.m'. This script is configured for a noise-free simulation (noise can be added in line 51, for example). The LST toolbox is self-contained and is ready-to-run.
