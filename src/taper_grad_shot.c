/*
 * taper gradient with a gaussian frame to damp inversion artefacts near the sources and receivers
 sws == 1 vertical taper (for tomography geometry)
 sws == 2 horizontal taper (for reflection geometry)
 sws == 3 local radial taper at the source and receiver positions
 sws == 4  
 */

#include "fd.h"

void taper_grad_shot(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int ishot)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY;
	extern FILE *FP;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2, VTON, SRTON;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1;

	/*extern int GRADT1, GRADT2, GRADT3, GRADT4;*/
	float amp, a, *window, grad_tap, **waveconvtmp;
	char modfile[STRING_SIZE];
	
	extern float SRTRADIUS;
	extern int SRTSHAPE, FILTSIZE;
        float **taper_coeff, taper;
	float R, frh, f0, x, y, r;
	
	taper_coeff= matrix(1,NY,1,NX);
	
	R = SRTRADIUS;
	frh = 0.1;
	f0 = 1e-3;
	a = - (log(frh-f0)-log(1.0-f0))/log(2.0);

        /*****************************/
        /* Taper at source positions */
        /*****************************/
	n=ishot;
	for (iy=1;iy<=NY;iy++){
             for (ix=1;ix<=NX;ix++){
	
	          /* calculate global coordinates */
		  x = ix * DH;
		  y = iy * DH; 
		  r = sqrt(pow((x-srcpos[1][n]),2.0)+pow((y-srcpos[2][n]),2.0));
		  taper = 1.0;
		  
		  if(r<=R){
		    taper = f0 + (1.0 - f0) * pow((0.5 - 0.5 *cos(PI*r/R)),a);		 
		  }
		  
		  taper_coeff[iy][ix] = taper;
		  waveconv[iy][ix] *= taper;
		       
             }
	}
    
        /*****************************/
        /* Taper at receiver positions */
        /*****************************/
	for (n=1;n<=ntr;n++) {
	     for (iy=1;iy<=NY;iy++){
                  for (ix=1;ix<=NX;ix++){
	
		       x = ix * DH;
		       y = iy * DH; 
		       r = sqrt(pow((x-(recpos[1][n]*DH)),2.0)+pow((y-(recpos[2][n]*DH)),2.0));
		       taper = 1.0;
		  
		       if(r<=R){
		         taper = f0 + (1.0 - f0) * pow((0.5 - 0.5 *cos(PI*r/R)),a);
		       }
		  
		           taper_coeff[iy][ix] = taper;
		           waveconv[iy][ix] *= taper;
		      
                  }
	     }
	}
	
	
	/*
        sprintf(modfile,"taper_coeff_%i.bin",ishot);
        writemod(modfile,taper_coeff,3); */
	
	free_matrix(taper_coeff,1,NY,1,NX);
	
}



