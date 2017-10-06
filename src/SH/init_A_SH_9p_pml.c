/*----------------------------------------------------------------------------------------------------------
 *  Initiate impedance matrix for 2D Helmholtz equation with PMLs 
 *  using a 9p-mixed grid according to 
 * 
 *  Hustedt, B., Operto, S., and Virieux, J. (2004). Mixed-grid and staggered-grid finite difference
 *  methods for frequency domain acoustic wave modelling. Geophysical Journal International, 157:1269–1296.
 *
 *  Operto, S., Virieux, J., Ribodetti, A. and Anderson, J.A. (2009) Finite-difference frequency-domain modeling 
 *  of viscoacoustic wave propagation in 2D tilted transversely isotropic (TTI) media . Geophysics 74(5):T75-T95.
 *
 *  D. Koehn
 *  Kiel, 16.08.2017
 *  --------------------------------------------------------------------------------------------------------*/

#include "fd.h"

void init_A_SH_9p_pml(struct PML_AC *PML_AC, struct matSH *matSH, struct waveAC *waveAC){

	extern int NX, NY, NONZERO;
        extern float DH, S;	

	/* local variables */
	int i, j, k;
        complex float tmp, Omega2;
        //complex float tmpA, tmpB, tmpC;
	complex float dyp, dym, dxp, dxm, dx0, dy0;
        float a, b, c, d, e, idh2, damp;
	//SuiteSparse_long count;
	int count, ishift;

	damp = 0.0;

        /* define FD parameters */
	a = 0.6287326;	
	b = 0.3712667;
	c = 1 - a - b;
	b = 0.25 * b;
        c = 0.25 * c;
        d = 0.4382634;
	e = 1 - d;

        idh2 = 1.0/(DH*DH);

	/* assemble impedance matrix */        

        k=0;        /* index of main diagonal  */	
        count=0;    /* count non-zero elements */

	/* set squared complex angular frequency*/
	Omega2 = cpowf(2.0 * M_PI * ((*waveAC).freq + S * I),2.0);

	init_mat_SH(waveAC,matSH);

        for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		    /* average shear modulus */
		    av_mu(matSH,i,j);

		    /* store PML parameters */
		    dy0 = (*PML_AC).dampyr[j] + (*PML_AC).dampyi[j] * I;
		    dx0 = (*PML_AC).dampxr[i] + (*PML_AC).dampxi[i] * I;	

		    dyp = (*PML_AC).dampyhr[j] + (*PML_AC).dampyhi[j] * I;
		    dym = (*PML_AC).dampyhr[j-1] + (*PML_AC).dampyhi[j-1] * I;

		    dxp = (*PML_AC).dampxhr[i] + (*PML_AC).dampxhi[i] * I;
		    dxm = (*PML_AC).dampxhr[i-1] + (*PML_AC).dampxhi[i-1] * I;

		    /* NW gridpoint */
		    if((i > 1) && (j > 1)){		

           	       tmp = c * (Omega2 * (*matSH).rho[j-1][i-1])
                           + d * ( 0.25 * idh2 * dy0 * (*matSH).bmm * dym
                                  +0.25 * idh2 * dx0 * (*matSH).bmm * dxm);
                
                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - 1 - NX;
                       count++;

		    }   

		    /* N gridpoint */
		    if(j > 1){		       

           	       tmp = b * (Omega2 * (*matSH).rho[j-1][i])
                	   + d * ( 0.25 * idh2 * dy0 * ((*matSH).bpm*dym+(*matSH).bmm*dym)
                		  +0.25 * idh2 * dx0 * (-(*matSH).bpm*dxp-(*matSH).bmm*dxm))
                	   + e * (dy0 * idh2 * (*matSH).bm0 * dym);

		       tmp += - damp * tmp;

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - NX;
                       count++;

		    }  

		    /* NE gridpoint */
		    if((i < NX) && (j > 1)){

          	       tmp = c * (Omega2 * (*matSH).rho[j-1][i+1])
                           + d * ( 0.25 * idh2 * dy0 * (*matSH).bpm * dym
                                  +0.25 * idh2 * dx0 * (*matSH).bpm * dxp);                

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1 - NX;
                       count++;

		    }  

		    /* W gridpoint */
		    if(i > 1){

		       tmp = b * (Omega2 * (*matSH).rho[j][i-1])
                	   + d * ( 0.25 * idh2 * dy0 * (-(*matSH).bmp*dyp-(*matSH).bmm*dym)
                                  +0.25 * idh2 * dx0 * ((*matSH).bmp*dxm+(*matSH).bmm*dxm))
                           + e * (dx0 * idh2 *(*matSH).bm0 * dxm);

	               tmp += - damp * tmp;

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k-1;
                       count++;

		    }

		    /* central gridpoint */

		    tmp = a * (Omega2 * (*matSH).rho[j][i]) 
                        + d * (-0.25 * idh2 * dy0 * ((*matSH).bmp*dyp + (*matSH).bpm*dym + (*matSH).bpp*dyp + (*matSH).bmm*dym)
                               -0.25 * idh2 * dx0 * ((*matSH).bmp*dxm + (*matSH).bpm*dxp + (*matSH).bpp*dxp + (*matSH).bmm*dxm))
                        + e * (-dy0 * idh2 * (dyp*(*matSH).b0p+dym*(*matSH).b0m) - dx0 * idh2 * (dxp*(*matSH).bp0+dxm*(*matSH).bm0));

		    tmp += damp * 4.0 * tmp;  

                    (*waveAC).Ar[count] = creal(tmp); 
                    (*waveAC).Ai[count] = cimag(tmp);
                    (*waveAC).irow[count] = k;
                    (*waveAC).icol[count] = k;
                    count++;

		    /* E gridpoint */
		    if(i < NX){

 		       tmp = b * (Omega2 * (*matSH).rho[j][i+1])
                           + d * ( 0.25 * idh2 * dy0 * (-(*matSH).bpm*dym-(*matSH).bpp*dyp)
                                  +0.25 * idh2 * dx0 * ((*matSH).bpm*dxp+(*matSH).bpp*dxp))
                           + e * (dx0 * idh2 * (*matSH).bp0 * dxp);

	               tmp += - damp * tmp;

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1;
                       count++;

		    }

		    /* SW gridpoint */
		    if( (j < NY) && (i > 1) ){

           	       tmp = c * (Omega2 * (*matSH).rho[j+1][i-1])
                           + d * ( 0.25 * idh2 * dy0 * (*matSH).bmp * dyp
                                  +0.25 * idh2 * dx0 * (*matSH).bmp * dxm);

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - 1 + NX;
                       count++;

		    }

		    /* S gridpoint */
		    if(j < NY){

                       tmp = b * (Omega2 * (*matSH).rho[j+1][i])
                           + d * (0.25 * idh2 * dy0 * ((*matSH).bmp*dyp+(*matSH).bpp*dyp)
                                  + 0.25 * idh2 * dx0 * (-(*matSH).bmp*dxm-(*matSH).bpp*dxp))
                	   + e * (dy0 * idh2 * (*matSH).b0p * dyp);

		       tmp += -damp * tmp; 

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + NX;
                       count++;

		    }     
            
		    /* SE gridpoint */
		    if((i < NX) && (j < NY)) {

                       tmp = c * (Omega2 * (*matSH).rho[j+1][i+1])
                	   + d * ( 0.25 * idh2 * dy0 * (*matSH).bpp * dyp
               			  +0.25 * idh2 * dx0 * (*matSH).bpp * dxp);                

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1 + NX;
                       count++;

		    }

                    k++;                       

		}
	}

}




