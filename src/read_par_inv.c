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
extern int SPATFILTER, NF;
extern float PRO, WD_DAMP, WD_DAMP1;
extern float FC_low, FC_high;
extern int FILT_SIZE_GRAD, FILT_SIZE_GRAD1;

/* definition of local variables */
int i;
char str [80];

fscanf(fp,"%s%s%s%s%s%s%s%s%s",str,str,str,str,str,str,str,str,str);
for (i=1;i<=nstage;i++){
      fscanf(fp,"%f%f%f%i%i%f%f%i%i",&PRO,&FC_low,&FC_high,&NF,&SPATFILTER,&WD_DAMP,&WD_DAMP1,&FILT_SIZE_GRAD,&FILT_SIZE_GRAD1);
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
      printf(" \t(Wavenumber domain damping with WD_DAMP = %e)\n",WD_DAMP); 
   }
   if(SPATFILTER==2){
      printf(" \t(Smooth2 damping with WD_DAMP = %e, WD_DAMP1 = %e)\n",WD_DAMP,WD_DAMP1); 
   }
   if(SPATFILTER==3){
      printf(" \t(Median filter with FILT_SIZE_GRAD = %d)\n",FILT_SIZE_GRAD); 
   }
   if(SPATFILTER==4){
      printf(" \t(Gaussian filter with FILT_SIZE_GRAD = %d and FILT_SIZE_GRAD1 = %d)\n",FILT_SIZE_GRAD,FILT_SIZE_GRAD1); 
   }
   if(SPATFILTER==0){
      printf(" \tSPATFILTER=%d: Gradients are not smoothed.\n",SPATFILTER);
   }

}

}
