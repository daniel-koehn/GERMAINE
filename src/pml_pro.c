/*-------------------------------------------------------------------------------
 *  Define damping profiles within CFS-PML boundary condition according to 
 * 
 *  Z. Chen, D. Cheng, W. Feng, H. Yang, 2013, An optimal 9-point finite 
 *  difference scheme for the Helmholtz equation with PML, Int. J. Numer. 
 *  Anal. Model., 10, 389-410.
 *  
 *  and
 *  
 *  W. Zhang, Y. Shen, 2013, Unsplit complex frequency-shifted PML 
 *  implementation using auxiliary differential equations for seismic 
 *  wave modeling. Geophysics, 75(4), T141–T154.
 *  
 *  D. Koehn
 *  Kiel, 20.09.2016
 *  -----------------------------------------------------------------------------*/

#include "fd.h"

void pml_pro(struct PML_AC *PML_AC, struct waveAC *waveAC){

	extern int NX, NY, NPML, FREE_SURF;
        extern float DH, A0_PML, OMEGA_PML;
	extern float A0_PML, OMEGA_PML, PML_VEL, PML_BETA0;
	extern float PML_POWD, PML_POWB, PML_POWA;
	extern char SNAP_FILE[STRING_SIZE];	

	/* local variables */
	int i, j, k;
        complex float tmp, sx, sy, ci;
        float lPML, d0, a0;
	char filename[STRING_SIZE];

        /* define complex i */
        ci = 0.0 + 1.0 *I;

	/* define thickness of PML [m] */
	lPML = NPML * DH;

	/* define PML damping parameters */
        OMEGA_PML = 2.0 * M_PI * (*waveAC).freq;

	/* define polynomial degrees */
	PML_POWD = 2.0;
	PML_POWA = 1.0;
	PML_POWB = 2.0;

	/* define d0 */
	d0 = - log(A0_PML) * PML_VEL*(PML_POWD + 1.0) / (2.0 * lPML);
	
	/* define a0 */
	a0 = M_PI * (*waveAC).freq;

        /* set free surface to zero */
        FREE_SURF = 0;

        /* calculate PML damping profiles sigma_x and sigma_y in x- and y-direction */
	for (i=1;i<=NX;i++){  
	    
	    (*PML_AC).d_x[i] = 0.0;
	    (*PML_AC).a_x[i] = 0.0;
	    (*PML_AC).b_x[i] = 1.0;
	    
	    /* define damping profile at left PML boundary */
	    if(i <= NPML){
	      (*PML_AC).d_x[i] = d0 * pow((NPML-i)*DH/lPML,PML_POWD);
	      (*PML_AC).a_x[i] = a0 * (1.0 - pow((NPML-i)*DH/lPML,PML_POWA));
	      (*PML_AC).b_x[i] = 1.0 + (PML_BETA0 - 1.0) * pow((NPML-i)*DH/lPML,PML_POWB);
	    }
	    
	    /* define damping profile at right PML boundary */
	    if(i >= NX - NPML){
	      (*PML_AC).d_x[i] = d0 * pow((DH*((NX-NPML)-i))/lPML,PML_POWD);
	      (*PML_AC).a_x[i] = a0 * (1.0 - pow((DH*((NX-NPML)-i))/lPML,PML_POWA));
	      (*PML_AC).b_x[i] = 1.0 + (PML_BETA0 - 1.0) * pow((DH*((NX-NPML)-i))/lPML,PML_POWB);
	    }
	    
	}

	for (j=1;j<=NY;j++){
	    
	    (*PML_AC).d_y[j] = 0.0;
	    (*PML_AC).a_y[j] = 0.0;	    
	    (*PML_AC).b_y[j] = 1.0;

	    /* define damping profile at top PML boundary */
	    if((j <= NPML)&&(FREE_SURF==0)){
	      (*PML_AC).d_y[j] = d0 * pow((NPML-j)*DH/lPML,PML_POWD);
	      (*PML_AC).a_y[j] = a0 * (1.0 - pow((NPML-j)*DH/lPML,PML_POWA));
	      (*PML_AC).b_y[j] = 1.0 + (PML_BETA0 - 1.0) * pow((NPML-j)*DH/lPML,PML_POWB);
	    }
	    
	    /* define damping profile at bottom PML boundary */
	    if(j >= NY - NPML){
	      (*PML_AC).d_y[j] = d0 * pow((DH*((NY-NPML)-j))/lPML,PML_POWD);
	      (*PML_AC).a_y[j] = a0 * (1.0 - pow((DH*((NY-NPML)-j))/lPML,PML_POWA));
	      (*PML_AC).b_y[j] = 1.0 + (PML_BETA0 - 1.0) * pow((DH*((NY-NPML)-j))/lPML,PML_POWB);
	    }
	    
	}

	/* calculate sx and sy */
	for (j=1;j<=NY;j++){
	    for (i=1;i<=NX;i++){  
		
		sx = (*PML_AC).b_x[i] + ((*PML_AC).d_x[i]/((*PML_AC).a_x[i] + OMEGA_PML*ci));
		sy = (*PML_AC).b_y[j] + ((*PML_AC).d_y[j]/((*PML_AC).a_y[j] + OMEGA_PML*ci));

		/* sx = 1.0 - ((*PML_AC).d_x[i] * ci / OMEGA_PML);
		sy = 1.0 - ((*PML_AC).d_y[j] * ci / OMEGA_PML); */

                /* calculate real and imaginary parts of factors A, B and C */
		(*PML_AC).Ar[j][i] = creal(sy/sx);
		(*PML_AC).Ai[j][i] = cimag(sy/sx);

		(*PML_AC).Br[j][i] = creal(sx/sy);
		(*PML_AC).Bi[j][i] = cimag(sx/sy);

		(*PML_AC).Cr[j][i] = creal(sx*sy);
		(*PML_AC).Ci[j][i] = cimag(sx*sy);

	    }
	}

	/* calculate arithmetic average values of A and B in x- and y-direction */

        /* initiate arrays */
	store_mat((*PML_AC).Ar,(*PML_AC).Axr,NX,NY);
	store_mat((*PML_AC).Ai,(*PML_AC).Axi,NX,NY);

	store_mat((*PML_AC).Br,(*PML_AC).Byr,NX,NY);
	store_mat((*PML_AC).Bi,(*PML_AC).Byi,NX,NY);

        /* arithmetic averages */
	for (j=2;j<=NY-1;j++){
	    for (i=2;i<=NX-1;i++){  
		
		(*PML_AC).Axr[j][i] = ((*PML_AC).Ar[j][i] + (*PML_AC).Ar[j][i+1]) / 2.0;        
		(*PML_AC).Axi[j][i] = ((*PML_AC).Ai[j][i] + (*PML_AC).Ai[j][i+1]) / 2.0;        

		(*PML_AC).Byr[j][i] = ((*PML_AC).Br[j][i] + (*PML_AC).Br[j+1][i]) / 2.0;
		(*PML_AC).Byi[j][i] = ((*PML_AC).Bi[j][i] + (*PML_AC).Bi[j+1][i]) / 2.0;		

	    }
	}

   /*sprintf(filename,"%s_Axr.bin",SNAP_FILE);
   writemod(filename,(*PML_AC).Ar,3);*/
       
}




