/*-------------------------------------------------------------------------------
 *  Define damping profiles within PML boundary condition according to 
 * 
 *  Z. Chen, D. Cheng, W. Feng, H. Yang, 2013, An optimal 9-point finite 
 *  difference scheme for the Helmholtz equation with PML, Int. J. Numer. 
 *  Anal. Model., 10, 389-410.
 *
 *  D. Koehn
 *  Kiel, 19.06.2016
 *  -----------------------------------------------------------------------------*/

#include "fd.h"

void pml_pro(struct PML_AC *PML_AC, struct waveAC *waveAC){

	extern int NX, NY, NPML, FREE_SURF;
        extern float DH, A0_PML, OMEGA_PML;
	extern char SNAP_FILE[STRING_SIZE];	

	/* local variables */
	int i, j, k;
        complex float tmp, sx, sy, ci;
        float lPML;
	char filename[STRING_SIZE];

        /* define complex i */
        ci = 0.0 + 1.0 *I;

	/* define thickness of PML [m] */
	lPML = NPML * DH;

	/* define PML damping parameters */
        OMEGA_PML = 2.0 * M_PI * (*waveAC).freq;

        /* set free surface to zero */
        FREE_SURF = 0;

        /* calculate PML damping profiles sigma_x and sigma_y in x- and y-direction */
	for (i=1;i<=NX;i++){  
	    
	    (*PML_AC).sigma_x[i] = 0.0;
	    
	    /* define damping profile at left PML boundary */
	    if(i < NPML){
	      (*PML_AC).sigma_x[i] = OMEGA_PML * A0_PML * pow((NPML-i)*DH/lPML,2.0);
	    }
	    
	    /* define damping profile at right PML boundary */
	    if(i >= NX - NPML){
	      (*PML_AC).sigma_x[i] = OMEGA_PML * A0_PML * pow((DH*((NX-NPML+1)-i))/lPML,2.0);
	    }
	    
	}

	for (j=1;j<=NY;j++){
	    
	    (*PML_AC).sigma_y[j] = 0.0;
	    
	    /* define damping profile at top PML boundary */
	    if((j < NPML)&&(FREE_SURF!=1)){
	      (*PML_AC).sigma_y[j] = OMEGA_PML * A0_PML * pow((NPML-j)*DH/lPML,2.0);
	    }
	    
	    /* define damping profile at bottom PML boundary */
	    if(j >= NY - NPML){
	      (*PML_AC).sigma_y[j] = OMEGA_PML * A0_PML * pow((DH*((NY-NPML+1)-j))/lPML,2.0);
	    }
	    
	}

	/* calculate sx and sy */
	for (j=1;j<=NY;j++){
	    for (i=1;i<=NX;i++){  
		
		sx = 1.0 - ((*PML_AC).sigma_x[i]*ci/OMEGA_PML);
		sy = 1.0 - ((*PML_AC).sigma_y[j]*ci/OMEGA_PML);

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
		
		(*PML_AC).Axr[j][i] = ((*PML_AC).Ar[j][i] + (*PML_AC).Ar[j][i+1])/2.0;        
		(*PML_AC).Axi[j][i] = ((*PML_AC).Ai[j][i] + (*PML_AC).Ai[j][i+1])/2.0;        

		(*PML_AC).Byr[j][i] = ((*PML_AC).Br[j][i] + (*PML_AC).Br[j+1][i])/2.0;
		(*PML_AC).Byi[j][i] = ((*PML_AC).Bi[j][i] + (*PML_AC).Bi[j+1][i])/2.0;		

	    }
	}

   /*sprintf(filename,"%s_Axr.bin",SNAP_FILE);
   writemod(filename,(*PML_AC).Ar,3);*/
       
}




