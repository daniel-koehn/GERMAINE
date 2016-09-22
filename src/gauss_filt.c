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
int i, j, ii, jj, filt_size_x, filt_size_y;
int i1, j1, filtsize, hfsx, hfsy;
float **model_tmp, **kernel, grad, normgauss, smooth_meter;
float conv;
float r, sigmax, sigmay, sx, sy, sum=0.0;
float vref, lam;

/* calculate wavelength of reference velocity */
lam = VREF / FC_low;

/* define filter size as fraction of reference velocity wavelength */
filt_size_x = round((FILT_SIZE_GRAD * lam)/DH);
filt_size_y = round((FILT_SIZE_GRAD1 * lam)/DH);

if (filt_size_x==0)	return;
if (!(filt_size_x % 2)) {
    if (filt_size_x > 0)	
        filt_size_x += 1;
    else			
        filt_size_x -= 1;
}

if (filt_size_y==0)	return;
if (!(filt_size_y % 2)) {
    if (filt_size_y > 0)	
        filt_size_y += 1;
    else			
        filt_size_y -= 1;
}


hfsx = abs(filt_size_x)/2;
sigmax = hfsx/2;
sx = 2.0 * sigmax *sigmax;

hfsy = abs(filt_size_y)/2;
sigmay = hfsy/2;
sy = 2.0 * sigmay *sigmay;

if(MYID==0){
   printf("\n hfsx: %d \n",hfsx);
   printf("\n hfsy: %d \n",hfsy);
}

model_tmp = matrix(-hfsy+1,NY+hfsy,-hfsx+1,NX+hfsx);
kernel=matrix(1,abs(filt_size_y),1,abs(filt_size_x));

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
         
         kernel[jj+hfsy+1][ii+hfsx+1] = exp(-((ii*ii)/sx) - ((jj*jj)/sy));
         sum += kernel[jj+hfsy+1][ii+hfsx+1];

     }
}

/* normalize kernel */
for (i=1;i<=filt_size_x;i++){
     for (j=1;j<=filt_size_y;j++){
         
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
free_matrix(kernel,1,abs(filt_size_y),1,abs(filt_size_x));
			
}

