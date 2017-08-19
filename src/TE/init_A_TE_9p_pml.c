/*----------------------------------------------------------------------------------------------------------
 *  Initiate impedance matrix for 2D TE-mode problem with PMLs 
 *  using a 9p-mixed grid according to 
 * 
 *  Hustedt, B., Operto, S., and Virieux, J. (2004). Mixed-grid and staggered-grid finite difference
 *  methods for frequency domain acoustic wave modelling. Geophysical Journal International, 157:1269–1296.
 *
 *  Operto, S., Virieux, J., Ribodetti, A. and Anderson, J.A. (2009) Finite-difference frequency-domain modeling 
 *  of viscoacoustic wave propagation in 2D tilted transversely isotropic (TTI) media . Geophysics 74(5):T75-T95.
 *
 *  Lavoué, F., Brossier, R., Métivier, L., Garambois, S. & Virieux, J. (2014) Two-dimensional permittivity 
 *  and conductivity imaging by full waveform inversion of multioffset GPR data: a frequency-domain quasi-Newton 
 *  approach, Geophysical Journal International, 197, 248-268.
 *
 *  D. Koehn
 *  Kiel, 19.08.2017
 *  --------------------------------------------------------------------------------------------------------*/

#include "fd.h"

void init_A_TE_9p_pml(struct PML_AC *PML_AC, struct matTE *matTE, struct waveAC *waveAC){

	extern int NX, NY, NONZERO;
        extern float DH, S;	

	/* local variables */
	int i, j, k;
        complex float tmp, Omega2, Omega, epse;
        //complex float tmpA, tmpB, tmpC;
	complex float dyp, dym, dxp, dxm, dx0, dy0;
        float a, b, c, d, e, idh2, damp, mu0;
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
	Omega = 2.0 * M_PI * ((*waveAC).freq - S * I);
	Omega2 = cpowf(2.0 * M_PI * ((*waveAC).freq - S * I),2.0);

	/* define magnetic permeability mu0 */
	mu0 = 4.0 * M_PI * 1e-7;

        for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		    /* store PML parameters */
		    dy0 = (*PML_AC).dampyr[j] + (*PML_AC).dampyi[j] * I;
		    dx0 = (*PML_AC).dampxr[i] + (*PML_AC).dampxi[i] * I;	

		    dyp = (*PML_AC).dampyhr[j] + (*PML_AC).dampyhi[j] * I;
		    dym = (*PML_AC).dampyhr[j-1] + (*PML_AC).dampyhi[j-1] * I;

		    dxp = (*PML_AC).dampxhr[i] + (*PML_AC).dampxhi[i] * I;
		    dxm = (*PML_AC).dampxhr[i-1] + (*PML_AC).dampxhi[i-1] * I;

		    /* NW gridpoint */
		    if((i > 1) && (j > 1)){		

		       epse = (*matTE).epsilon[j-1][i-1] + ((*matTE).sigma[j-1][i-1] / Omega) * I;

           	       tmp = c * (Omega2 * mu0 * epse)
                           + d * ( 0.25 * idh2 * dy0 * dym
                                  +0.25 * idh2 * dx0 * dxm);
                
                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - 1 - NX;
                       count++;

		    }   

		    /* N gridpoint */
		    if(j > 1){		       

		       epse = (*matTE).epsilon[j-1][i] + ((*matTE).sigma[j-1][i] / Omega) * I;

           	       tmp = b * (Omega2 * mu0 * epse)
                	   + d * ( 0.25 * idh2 * dy0 * (2.0*dym)
                		  +0.25 * idh2 * dx0 * (-dxp-dxm))
                	   + e * (dy0 * idh2 * dym);

		       tmp += - damp * tmp;

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - NX;
                       count++;

		    }  

		    /* NE gridpoint */
		    if((i < NX) && (j > 1)){

		       epse = (*matTE).epsilon[j-1][i+1] + ((*matTE).sigma[j-1][i+1] / Omega) * I;

          	       tmp = c * (Omega2 * mu0 * epse)
                           + d * ( 0.25 * idh2 * dy0 * dym
                                  +0.25 * idh2 * dx0 * dxp);                

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1 - NX;
                       count++;

		    }  

		    /* W gridpoint */
		    if(i > 1){

		       epse = (*matTE).epsilon[j][i-1] + ((*matTE).sigma[j][i-1] / Omega) * I;

		       tmp = b * (Omega2 * mu0 * epse)
                	   + d * ( 0.25 * idh2 * dy0 * (-dyp-dym)
                                  +0.25 * idh2 * dx0 * (2.0*dxm))
                           + e * (dx0 * idh2 * dxm);

	               tmp += - damp * tmp;

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k-1;
                       count++;

		    }

		    /* central gridpoint */

		    epse = (*matTE).epsilon[j][i] + ((*matTE).sigma[j][i] / Omega) * I;

		    tmp = a * (Omega2 * mu0 * epse) 
                        + d * (-0.25 * idh2 * dy0 * (dyp + dym + dyp + dym)
                               -0.25 * idh2 * dx0 * (dxm + dxp + dxp + dxm))
                        + e * (-dy0 * idh2 * (dyp + dym) - dx0 * idh2 * (dxp+dxm));

		    tmp += damp * 4.0 * tmp;  

                    (*waveAC).Ar[count] = creal(tmp); 
                    (*waveAC).Ai[count] = cimag(tmp);
                    (*waveAC).irow[count] = k;
                    (*waveAC).icol[count] = k;
                    count++;

		    /* E gridpoint */
		    if(i < NX){

		       epse = (*matTE).epsilon[j][i+1] + ((*matTE).sigma[j][i+1] / Omega) * I;

 		       tmp = b * (Omega2 * mu0 * epse)
                           + d * ( 0.25 * idh2 * dy0 * (-dym-dyp)
                                  +0.25 * idh2 * dx0 * (dxp+dxp))
                           + e * (dx0 * idh2 * dxp);

	               tmp += - damp * tmp;

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1;
                       count++;

		    }

		    /* SW gridpoint */
		    if( (j < NY) && (i > 1) ){

		       epse = (*matTE).epsilon[j+1][i-1] + ((*matTE).sigma[j+1][i-1] / Omega) * I;

           	       tmp = c * (Omega2 * mu0 * epse)
                           + d * ( 0.25 * idh2 * dy0 * dyp
                                  +0.25 * idh2 * dx0 * dxm);

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - 1 + NX;
                       count++;

		    }

		    /* S gridpoint */
		    if(j < NY){

		       epse = (*matTE).epsilon[j+1][i] + ((*matTE).sigma[j+1][i] / Omega) * I;

                       tmp = b * (Omega2 * mu0 * epse)
                           + d * (0.25 * idh2 * dy0 * (dyp+dyp)
                                  + 0.25 * idh2 * dx0 * (-dxm-dxp))
                	   + e * (dy0 * idh2 * dyp);

		       tmp += -damp * tmp; 

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + NX;
                       count++;

		    }     
            
		    /* SE gridpoint */
		    if((i < NX) && (j < NY)) {

		       epse = (*matTE).epsilon[j+1][i+1] + ((*matTE).sigma[j+1][i+1] / Omega) * I;

                       tmp = c * (Omega2 * mu0 * epse)
                	   + d * ( 0.25 * idh2 * dy0 * dyp
               			  +0.25 * idh2 * dx0 * dxp);                

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




