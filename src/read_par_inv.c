/*------------------------------------------------------------------------
 *   Read FWI-workflow-parameters from Stdin                          
 *
 *  D. Koehn
 *  Kiel, the 8th of August 2013
 *  ----------------------------------------------------------------------*/

/* reading FWI-workflow parameters for DENISE */

#include "fd.h"

void read_par_inv(FILE *fp,int nstage,int stagemax){

/* declaration of global variables */
extern int MYID, INVMAT, PHYSICS;
extern int SPATFILTER, NF, MISFIT;
extern float PRO, FC_low, FC_high, VREF, S, BETA_MAT1, BETA_MAT2;
extern float FILT_SIZE_GRAD, FILT_SIZE_GRAD1, LAMBDA_1, LAMBDA_2;
extern float MAT1_NORM, MAT2_NORM, MAT1_NORM0, MAT2_NORM0;

/* definition of local variables */
int i;
char str [80];

fscanf(fp,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s",str,str,str,str,str,str,str,str,str,str,str,str,str,str);
for (i=1;i<=nstage;i++){
      fscanf(fp,"%f%f%f%i%f%i%f%f%f%i%f%f%f%f",&PRO,&FC_low,&FC_high,&NF,&S,&SPATFILTER,&FILT_SIZE_GRAD,&FILT_SIZE_GRAD1,&VREF,&MISFIT,&BETA_MAT1,&BETA_MAT2,&LAMBDA_1,&LAMBDA_2);
}

fclose(fp);

/* calculate scaling factors on all MPI processes*/
if(BETA_MAT1>0.0){MAT1_NORM = MAT1_NORM0 * BETA_MAT1;}
if(BETA_MAT2>0.0){MAT2_NORM = MAT2_NORM0 * BETA_MAT2;}

if((MYID==0)&&(INVMAT==1)){

   printf("=========================================== \n");
   printf("       GERMAINE-stage %d of %d \n",nstage,stagemax);
   printf("=========================================== \n");
   
   printf("\n\n");
   printf(" Invert NF = %d frequencies between FC_low = %e Hz and FC_high = %e Hz.\n", NF, FC_low, FC_high);
   printf(" Laplace constant S = %e 1/s\n", S);
   printf(" Smoothing (spatial filtering) of the gradients: \n ");
   if(SPATFILTER){
     printf(" \tSPATFILTER=%d: Gradients are smoothed.\n",SPATFILTER);
   }

   if(SPATFILTER==1){
      printf(" \tGaussian filter with FILT_SIZE_GRAD = %f and FILT_SIZE_GRAD1 = %f fraction of average wavelength\n",FILT_SIZE_GRAD,FILT_SIZE_GRAD1); 
      printf(" \tusing reference velocity VREF = %f\n",VREF); 
   }
   printf("\n\n");
   printf("Misfit function: \n ");
   if(MISFIT==1){
      printf(" \t L2-norm MISFIT = %d\n",MISFIT); 
   }
   if(MISFIT==2){
      printf(" \t logarithmic L2-norm MISFIT = %d\n",MISFIT); 
   }
   if(MISFIT==3){
      printf(" \t phase only logarithmic L2-norm MISFIT = %d\n",MISFIT); 
   }
   printf("\n\n");

   printf(" Tikhonov regularization: \n ");
   if(PHYSICS==1){   
       printf(" BETA_VP = %e \t BETA_RHO = %e \n", BETA_MAT1, BETA_MAT2);
       printf(" VP_NORM = %e \t RHO_NORM = %e \n", MAT1_NORM, MAT2_NORM);
       printf(" LAMBDA_VP = %e \t LAMBDA_RHO = %e \n", LAMBDA_1, LAMBDA_2);
   }

   if(PHYSICS==4){   
       printf(" BETA_SIGMA = %e \t BETA_EPSILON = %e \n", BETA_MAT1, BETA_MAT2);
       printf(" SIGMA_NORM = %e \t EPSILON_NORM = %e \n", MAT1_NORM, MAT2_NORM);
       printf(" LAMBDA_SIGMA = %e \t LAMBDA_EPSILON = %e \n", LAMBDA_1, LAMBDA_2);
   }

}

}
