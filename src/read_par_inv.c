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
extern int MYID, INVMAT;
extern int SPATFILTER, NF, MISFIT;
extern float PRO, FC_low, FC_high, VREF;
extern float FILT_SIZE_GRAD, FILT_SIZE_GRAD1;

/* definition of local variables */
int i;
char str [80];

fscanf(fp,"%s%s%s%s%s%s%s%s%s",str,str,str,str,str,str,str,str,str);
for (i=1;i<=nstage;i++){
      fscanf(fp,"%f%f%f%i%i%f%f%f%i",&PRO,&FC_low,&FC_high,&NF,&SPATFILTER,&FILT_SIZE_GRAD,&FILT_SIZE_GRAD1,&VREF,&MISFIT);
}

fclose(fp);

if((MYID==0)&&(INVMAT==1)){
   printf("=========================================== \n");
   printf("       GERMAINE-stage %d of %d \n",nstage,stagemax);
   printf("=========================================== \n");
   
   printf("\n\n");
   printf(" Invert NF = %d frequencies between FC_low = %e Hz and FC_high = %e Hz.\n", NF, FC_low, FC_high);
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

}

}
