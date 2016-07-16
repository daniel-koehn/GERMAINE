/*------------------------------------------------------------------------
 * Step length estimation by parabolic line search
 *
 * D. Koehn
 * Kiel, 11.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float step_length_est(float ** waveconv, float ** Vp, float **S, int iter, int partest, 
float ** srcpos, int ntr, int ** recpos, int itest, int *step1, int *step3, int nfstart, int nshot1, 
int nshot2, int nsrc, float * epst1, float * L2t){

extern int MIN_ITER,STEPMAX, NX, NY, MYID;
extern char JACOBIAN[STRING_SIZE];
extern float EPS_SCALE, SCALEFAC;

float opteps_vp, *Tmod, *Tobs, *Tres, **Vpnp1, **TT, **Snp1;
int h, i, j, n, ishot;

/* Variables for step length calculation */
int step2, itests, iteste, stepmax, countstep;
float scalefac, eps_scale, L2sum, tmp, L2;

Vpnp1 =  matrix(1,NY,1,NX);
Tmod = vector(1,ntr);
Tobs = vector(1,ntr);
Tres = vector(1,ntr);
TT =  matrix(1,NY,1,NX);
Snp1 =  matrix(1,NY,1,NX);

scalefac = SCALEFAC;  /* scale factor for the step length */
stepmax  = STEPMAX;   /* number of maximum misfit calculations/steplength 2/3*/ 

*step1=0;
step2=0;

/* start with first guess for step length alpha */
eps_scale=EPS_SCALE; /* maximum model change = 1% of the maximum model value */
countstep=0; /* count number of forward calculations */

itests=2;
iteste=2;

while((step2!=1)||(*step1!=1)){

      for(itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

          /* test Vp-update */
	  tmp=calc_mat_change(waveconv,Vp,Vpnp1,iter,eps_scale,1,nfstart);

          /* calculate slowness S from test Vp model*/
          calc_S(Vpnp1,Snp1);
 
          L2 = 0.0;
          for (ishot=nshot1;ishot<nshot2;ishot++){         

		  /* solve forward problem with Eikonal solver */
		  eikonal(Snp1,TT,Tmod,srcpos,ishot,nsrc,recpos,ntr);

                  /* calculate traveltime residuals at receiver positions */
 	          L2+=calc_FA_res(Tmod,Tobs,Tres,ntr,ishot);

          }

          /* assemble objective function from all MPI processes */
	  L2sum = 0.0;
          MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	  L2t[itest]= 0.5 * L2sum;
	     
      } /* end of L2 test */

      /* Did not found a step size which reduces the misfit function */
      if((*step1==0)&&(L2t[1]<=L2t[2])){

        eps_scale = eps_scale/scalefac; 
        countstep++;

      }

      /* Found a step size with L2t[2] < L2t[3]*/
      if((*step1==1)&&(L2t[2]<L2t[3])){

        epst1[3]=eps_scale;
        step2=1;

      }

      /* Could not found a step size with L2t[2] < L2t[3]*/
      if((*step1==1)&&(L2t[2]>=L2t[3])){

         epst1[3]=eps_scale;

         /* increase step length to find  a larger misfit function than L2t[2]*/
         eps_scale = eps_scale + (eps_scale/scalefac);
         countstep++;                       
      }         

      /* found a step size which reduces the misfit function */
      if((*step1==0)&&(L2t[1]>L2t[2])){

         epst1[2]=eps_scale; 
         *step1=1;
         iteste=3;
         itests=3;
         countstep=0;
         /* find a second step length with a larger misfit function than L2t[2]*/
         eps_scale = eps_scale + (eps_scale/scalefac);

      }

      *step3=0;

      if((*step1==0)&&(countstep>stepmax)){

         if(MYID==0){
            printf(" Steplength estimation failed!");} 
         *step3=1;
         break;

      }

      if((*step1==1)&&(countstep>stepmax)){

         if(MYID==0){
            printf("Could not found a proper 3rd step length which brackets the minimum\n");}
            *step1=1;
            step2=1;

      }

if(MYID==0){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e\n",iteste,itests,*step1,step2,eps_scale,countstep,stepmax,scalefac,L2t[1],L2t[2],L2t[3]);}

} /* end of while loop */

if(*step1==1){ /* only find an optimal step length if step1==1 */

   /* calculate optimal step length epsilon for Vp and Vs*/
   if(MYID==0){
      printf("===================================================== \n");
      printf("    calculate optimal step length epsilon for Vp      \n");
      printf("===================================================== \n");
   }

   opteps_vp = calc_opt_step(L2t,epst1,1);
   eps_scale = opteps_vp;

}

free_matrix(Vpnp1,1,NY,1,NX);
free_matrix(Snp1,1,NY,1,NX);
free_matrix(TT,1,NY,1,NX);
free_vector(Tmod,1,ntr);
free_vector(Tobs,1,ntr);
free_vector(Tres,1,ntr);

return eps_scale;
}

