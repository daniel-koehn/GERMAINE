/*
 * Application of Gaussian filter to gradient
 * 
 * D. Koehn
 * Kiel, 12.12.2015
 */

#include "fd.h"

void gauss_filt(float ** waveconv){

/* global variables */
extern float DH, VREF, FC_low;
extern int FREE_SURF, NX, NY, MYID;
extern char JACOBIAN[STRING_SIZE];
extern float FILT_SIZE_GRAD, FILT_SIZE_GRAD1;

/* local variables */
int i, j, ii, jj;
int i1, j1, filtsize, hfsx, hfsy;
float **model_tmp, **kernel, grad, normgauss, smooth_meter;
float conv;
float r, sigmax, sigmay, sx, sy, sum=0.0;
float vref, lam, frac_lam_x, frac_lam_y, xlam, ylam;

/* calculate wavelength of reference velocity */
lam = VREF / FC_low;

/* define filter size as fraction of reference velocity wavelength */
frac_lam_x = FILT_SIZE_GRAD * lam;
frac_lam_y = FILT_SIZE_GRAD1 * lam;

/* calculate filter size */
xlam = 3.0 * frac_lam_x;
ylam = 3.0 * frac_lam_y;

hfsx = (int)(xlam/DH);
hfsy = (int)(ylam/DH);

sx = frac_lam_x * frac_lam_x;
sy = frac_lam_y * frac_lam_y;

if(MYID==0){
   printf("\n hfsx: %d \n",hfsx);
   printf("\n hfsy: %d \n",hfsy);
}

model_tmp = matrix(-hfsy+1,NY+hfsy,-hfsx+1,NX+hfsx);
kernel=matrix(1,2*hfsy+1,1,2*hfsx+1);

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
     for (i=-hfsx+1;i<=0;i++){
          model_tmp[j][i] = model_tmp[j][1];
     }

     for (i=NX+1;i<=NX+hfsx;i++){
	  model_tmp[j][i] = model_tmp[j][NX];
     }

}

/* top/bottom boundary incl. corners */
for (j=-hfsy+1;j<=0;j++){

      for (i=-hfsx+1;i<=NX+hfsx;i++){
	model_tmp[j][i] = model_tmp[1][i];
      }

}
for (j=NY+1;j<=NY+hfsy;j++){

      for (i=-hfsx+1;i<=NX+hfsx;i++){
	model_tmp[j][i] = model_tmp[NY][i];
      }

}

/* create filter kernel */
for (ii=-hfsx;ii<=hfsx;ii++){
     for (jj=-hfsy;jj<=hfsy;jj++){
         
	 kernel[jj+hfsy+1][ii+hfsx+1] = exp(-((ii*ii*DH*DH)/sx) - ((jj*jj*DH*DH)/sy));
         sum += kernel[jj+hfsy+1][ii+hfsx+1];

     }
}

/* normalize kernel */
for (i=1;i<=2*hfsx;i++){
     for (j=1;j<=2*hfsy;j++){
         
         kernel[j][i] /= sum;

     }
}

/* apply Gaussian filter to gradient */
for (j=1;j<=NY;j++){
     for (i=1;i<=NX;i++){

          conv = 0.0;
          /* loop over kernel*/
	  for (ii=-hfsx;ii<=hfsx;ii++){
	       for (jj=-hfsy;jj<=hfsy;jj++){

	            conv += model_tmp[j+jj][i+ii]*kernel[jj+hfsy+1][ii+hfsx+1];				

               }
          }

          /* output of filtered gradient */
          waveconv[j][i] = conv;

      }
}
      
free_matrix(model_tmp,-hfsy+1,NY+hfsy,-hfsx+1,NX+hfsx);
free_matrix(kernel,1,2*hfsy+1,1,2*hfsx+1);
			
}

