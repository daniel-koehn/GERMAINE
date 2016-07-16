# GERMAINE

2D acoustic FDFD modelling and FWI code, which I developed together with Denise De Nil.

The 2D Frequency Domain Finite-Difference (FDFD) Code GERMAINE solves the 2D Helmholtz equation using a 9-point FD stencil with PML absorbing boundary conditions according to 

- I. Singer, E. Turkel, 2004, A perfectly matched layer for the Helmholtz equation in a semi-infinite strip. Journal of Computational Physics, 201(2), 439-465.
- Z. Chen, D. Cheng, W. Feng, H. Yang, 2013, An optimal 9-point finite difference scheme for the Helmholtz equation with PML, Int. J. Numer. Anal. Model., 10, 389-410. 

The forward wavefield is calculated via a LU-decompostion and forward/backward substitution using UMFPACK, which is part of the sparse matrix library SuiteSparse:

http://faculty.cse.tamu.edu/davis/suitesparse.html

The code is parallelized with MPI using a very simple shot parallelization. The FWI code is based on the adjoint state-method with quasi-Newton l-BFGS optimization (Nocedal & Wright 2006).

GERMAINE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2.0 of the License only.

GERMAINE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License in LICENSE.md for more details.

If you show modelling/inversion results in a paper or presentation please give a reference to the following papers:


Daniel Koehn and Denise De Nil
