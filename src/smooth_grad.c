/*
 * Smoothing model with a median frame to damp inversion artefacts near the sources and receivers
 * M. Schaefer - more or less copied from spat_filt.c - May 2011
 * Kiel, 12.12.2015 - modified for RAJZEL (D. Koehn)
 */

#include "fd.h"

void smooth_grad(float ** waveconv)
{

/* global variables */
extern float DH;
extern int FREE_SURF, NX, NY;
extern char JACOBIAN[STRING_SIZE];
extern int FILT_SIZE_GRAD;

/* local variables */
int i, j, ii, jj;
int i1, j1, filtsize, hfs;

float **model_tmp, **filterpart, grad, normgauss, smooth_meter;
	
	
if (FILT_SIZE_GRAD==0)	return;
if (!(FILT_SIZE_GRAD % 2)) {
if (FILT_SIZE_GRAD > 0)	FILT_SIZE_GRAD += 1;
else			FILT_SIZE_GRAD -= 1;
}

hfs = abs(FILT_SIZE_GRAD)/2;
printf("\n hfs: %d \n",hfs);
model_tmp = matrix(-hfs+1,NY+hfs,-hfs+1,NX+hfs);

filterpart=matrix(1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));

/* load merged model */
for (i=1;i<=NX;i++){
	for (j=1;j<=NY;j++){
	      	model_tmp[j][i]=waveconv[j][i];
      	}	
}

/* apply 2D-Gaussian filter on vp model */

/* extrapolate array */
/* left/right boundary */
for (j=1;j<=NY;j++){
     for (i=-hfs+1;i<=0;i++){
          model_tmp[j][i] = model_tmp[j][1];
     }

     for (i=NX+1;i<=NX+hfs;i++){
	  model_tmp[j][i] = model_tmp[j][NX];
     }

}

/* top/bottom boundary incl. corners */
for (j=-hfs+1;j<=0;j++){

      for (i=-hfs+1;i<=NX+hfs;i++){
	model_tmp[j][i] = model_tmp[1][i];
      }

}
for (j=NY+1;j<=NY+hfs;j++){

      for (i=-hfs+1;i<=NX+hfs;i++){
	model_tmp[j][i] = model_tmp[NY][i];
      }

}

/* filter */
for (j=1;j<=NY;j++){
     for (i=1;i<=NX;i++){

	  /* create a filtersize x filtersize matrix */
	  for (ii=-hfs;ii<=hfs;ii++){
	       for (jj=-hfs;jj<=hfs;jj++){

			      /*if ((ii+hfs+1)<(-hfs+1)) err(" (ii+hfs+1)<(-hfs+1) ! ");
			      if ((ii+hfs+1)>(hfs+NXG)) err(" (ii+hfs+1)>(hfs+NXG) ! ");*/
			      filterpart[jj+hfs+1][ii+hfs+1] = model_tmp[j+jj][i+ii];

	       }
	  }

	  /* filter */
	  waveconv[j][i] = median2d(filterpart,abs(FILT_SIZE_GRAD),abs(FILT_SIZE_GRAD));				

      }
}
      
free_matrix(model_tmp,-hfs+1,NY+hfs,-hfs+1,NX+hfs);
free_matrix(filterpart,1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
			
}/* end of application condition for the smoothing */
