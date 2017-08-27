/*------------------------------------------------------------------------
 * Step length estimation by parabolic line search
 *
 * D. Koehn
 * Kiel, 22.08.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float parabolicls_TE(struct fwiTE *fwiTE, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matTE *matTE, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage, float alpha, float L2){

extern int MIN_ITER,STEPMAX, NX, NY, MYID;
extern char JACOBIAN[STRING_SIZE];
extern float EPS_SCALE, SCALEFAC, MAT1_NORM, MAT2_NORM;

float opteps_vp, L2_tmp;
int h, i, j, n, ishot;

/* Variables for step length calculation */
int step1, step2, step3, itests, iteste, stepmax, countstep;
int itest;
float scalefac, tmp, eps_scale;
float maxgrad_sigmar, maxsigmar, maxgrad_epsilonr, maxepsilonr, *L2t, *epst1;

L2t = vector(1,3);
epst1 = vector(1,3);

/* normalize material parameters */
scale_grad((*matTE).sigma,1.0/MAT1_NORM,(*matTE).sigmar,NX,NY);
scale_grad((*matTE).epsilon,1.0/MAT2_NORM,(*matTE).epsilonr,NX,NY);

scalefac = SCALEFAC;  /* scale factor for the step length */
stepmax  = STEPMAX;   /* number of maximum misfit calculations/steplength 2/3*/ 

step1=0;
step2=0;

/* start with first guess for step length alpha */
maxgrad_sigmar = maximum_m((*fwiTE).Hgrad_sigma,NX,NY);
     maxsigmar = maximum_m((*matTE).sigmar,NX,NY);
        alpha = EPS_SCALE * maxsigmar/maxgrad_sigmar;

countstep=0; /* count number of forward calculations */
L2t[1] = L2;

itests=2;
iteste=2;

if(MYID==0){
  printf("\n========================================================================== \n");
  printf("\n ***** Estimate optimum step length by parabolic line search   **********  \n");
  printf("\n========================================================================== \n\n"); 
}

while((step2!=1)||(step1!=1)){

      for(itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

          /* copy sigmar -> sigma_old */
          /* copy epsilonr -> epsilon_old */
	  store_mat((*matTE).sigmar,(*fwiTE).sigma_old,NX,NY);
	  store_mat((*matTE).epsilonr,(*fwiTE).epsilon_old,NX,NY);

          /* test sigmar and epsilonr-update */
	  calc_mat_change_wolfe_multi_para((*fwiTE).Hgrad_sigma,(*matTE).sigmar,(*fwiTE).sigma_old,alpha,1);
	  calc_mat_change_wolfe_multi_para((*fwiTE).Hgrad_epsilon,(*matTE).epsilonr,(*fwiTE).epsilon_old,alpha,2);

          /* convert sigmar -> sigma and epsilonr -> epsilon to perform forward modelling */
	  scale_grad((*matTE).sigmar,MAT1_NORM,(*matTE).sigma,NX,NY);
	  scale_grad((*matTE).epsilonr,MAT2_NORM,(*matTE).epsilon,NX,NY);

          L2_tmp = obj_TE(fwiTE,waveAC,PML_AC,matTE,srcpos,nshots,recpos,ntr,iter,nstage);

	  /* Tikhonov regularization (cost function) */
	  L2t[itest] = Tikhonov_cost_TE(fwiTE,matTE,L2_tmp,iter);

	  /* copy sigma_old -> sigmar */
	  /* copy epsilon_old -> epsilonr */
	  store_mat((*fwiTE).sigma_old,(*matTE).sigmar,NX,NY);
	  store_mat((*fwiTE).epsilon_old,(*matTE).epsilonr,NX,NY);
	     
      } /* end of L2 test */

      /* Did not found a step size which reduces the misfit function */
      if((step1==0)&&(L2t[1]<=L2t[2])){

        alpha = alpha/scalefac; 
        countstep++;

      }

      /* Found a step size with L2t[2] < L2t[3] */
      if((step1==1)&&(L2t[2]<L2t[3])){

        epst1[3]=alpha;
        step2=1;

      }

      /* Could not found a step size with L2t[2] < L2t[3]*/
      if((step1==1)&&(L2t[2]>=L2t[3])){

         epst1[3]=alpha;

         /* increase step length to find  a larger misfit function than L2t[2]*/
         alpha = alpha + (alpha/scalefac);
         countstep++;                       
      }         

      /* found a step size which reduces the misfit function */
      if((step1==0)&&(L2t[1]>L2t[2])){

         epst1[2]=alpha; 
         step1=1;
         iteste=3;
         itests=3;
         countstep=0;
	 
         /* find a second step length with a larger misfit function than L2t[2]*/
         alpha = alpha + (alpha/scalefac);

      }

      step3=0;

      if((step1==0)&&(countstep>stepmax)){

         if(MYID==0){printf(" Steplength estimation failed!");} 
         step3=1;
	 eps_scale = 0.0;
         break;

      }

      if((step1==1)&&(countstep>stepmax)){

         if(MYID==0){
            printf("Could not found a proper 3rd step length which brackets the minimum\n");}
            step1=1;
            step2=1;

      }

if(MYID==0){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t alpha = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e\n",iteste,itests,step1,step2,alpha,countstep,stepmax,scalefac,L2t[1],L2t[2],L2t[3]);}

} /* end of while loop */

if(step1==1){ /* only find an optimal step length if step1==1 */

   /* calculate optimal step length epsilon for Vp*/
   if(MYID==0){
      printf("===================================================== \n");
      printf("    calculate optimal step length epsilon for Vp      \n");
      printf("===================================================== \n");
   }

   opteps_vp = calc_opt_step(L2t,epst1,1);
   eps_scale = opteps_vp;

}

/* convert sigmar -> sigma and epsilonr -> epsilon to perform forward modelling */
scale_grad((*matTE).sigmar,MAT1_NORM,(*matTE).sigma,NX,NY);
scale_grad((*matTE).epsilonr,MAT2_NORM,(*matTE).epsilon,NX,NY);

free_vector(L2t,1,3);
free_vector(epst1,1,3);

return eps_scale;

}

