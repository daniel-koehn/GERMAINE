# GERMAINE

2D acoustic FDFD modelling, FWI and RTM code, which I developed together with Denise De Nil.

The Frequency Domain Finite-Difference (FDFD) Code GERMAINE solves the 2D isotropic acoustic/TE-mode wave equation using a 9-point mixed-grid FD stencil with PML absorbing boundary conditions according to 

- Hustedt, B., Operto, S., and Virieux, J. (2004) Mixed-grid and staggered-grid finite difference methods for frequency domain acoustic wave modelling. Geophysical Journal International, 157:1269–1296.
- Operto, S., Virieux, J., Ribodetti, A. and Anderson, J.A. (2009) Finite-difference frequency-domain modeling of viscoacoustic wave propagation in 2D tilted transversely isotropic (TTI) media. Geophysics 74(5):T75-T95.

The forward wavefield is calculated via a LU-decompostion and forward/backward substitution using UMFPACK, which is part of the sparse matrix library SuiteSparse:

http://faculty.cse.tamu.edu/davis/suitesparse.html

The code is parallelized with MPI using a very simple shot-frequency parallelization. The FWI code is based on the adjoint state-method with quasi-Newton l-BFGS optimization (Nocedal & Wright 2006).

GERMAINE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2.0 of the License only.

GERMAINE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License in LICENSE.md for more details.

If you show modelling/inversion results in a paper or presentation please refer to:

Köhn, D., De Nil, D. and Rabbel, W. (2017) Tutorial: Introduction to frequency domain modelling and FWI of georadar data with GERMAINE, DOI:10.13140/RG.2.2.29354.03523

https://danielkoehnsite.wordpress.com/blog/tutorials/fdfd-modelling-and-fwi-of-georadar-data/

Daniel Koehn and Denise De Nil
