/*------------------------------------------------------------------------------------
 *  Calculate step length which satisfies the Wolfe condition 
 *  from T. van Leeuwen: https://github.com/tleeuwen/Penalty-Method/
 *  (adapted from http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3)
 *  
 *  D. Koehn
 *  Kiel, 23.08.2017
 *  ----------------------------------------------------------------------------------*/

#include "fd.h"

float wolfels_TE(struct fwiTE *fwiTE, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matTE *matTE, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage, float alpha, float L2){

	/* declaration of global variables */
        extern int NX, NY, GRAD_METHOD, STEPMAX, MYID;
        extern float EPS_SCALE, C1, C2, SCALEFAC;
	extern float MAT1_NORM, MAT2_NORM; 

        /* declaration of local variables */
        float normg, mu, nu, ft, ** Snp1;
        float ** gt_sigma, ** gt_epsilon, g0s0, gts0, maxgrad, maxsigma, tmp;     
	int lsiter, done;

        gt_sigma =  matrix(1,NY,1,NX);
        gt_epsilon =  matrix(1,NY,1,NX);

        /* copy -> gt */
	store_mat((*fwiTE).grad_sigma,gt_sigma,NX,NY);
	store_mat((*fwiTE).grad_epsilon,gt_epsilon,NX,NY);

	/* normalize material parameters */
	scale_grad((*matTE).sigma,1.0/MAT1_NORM,(*matTE).sigmar,NX,NY);
	scale_grad((*matTE).epsilon,1.0/MAT2_NORM,(*matTE).epsilonr,NX,NY);

        /* calculate initial step length */
        if(iter==1){

           /*normg = norm_matrix((*fwiTE).grad_sigma,NX,NY);
           normg += norm_matrix((*fwiTE).grad_epsilon,NX,NY);
           alpha = 1.0/normg;*/

           /*maxgrad = maximum_m((*fwiTE).Hgrad_sigma,NX,NY);
           maxsigma = maximum_m((*matTE).sigmar,NX,NY);

           alpha = EPS_SCALE * maxsigma/maxgrad;*/

	   alpha = 1.0;

        }else{

           alpha = 1.0;

        }

        lsiter = 0;
	done = 0;
	mu = 0.0;
        nu = 1e30;        

        if(MYID==0){
           printf(" Estimate step length by Wolfe line search \n");
           printf(" ----------------------------------------- \n");
        }

	while(done!=1){
    
              if(lsiter < STEPMAX){

                 if(MYID==0){
                    printf("gradient evaluation %d ... \n",lsiter);
                    printf("alpha = %e ... \n",alpha);
                 }

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

		 ft = grad_obj_TE(fwiTE,waveAC,PML_AC,matTE,srcpos,nshots,recpos,ntr,iter,nstage);
		 
		 /*descent((*fwiTE).grad_sigma,gt_sigma);
 		 descent((*fwiTE).grad_epsilon,gt_epsilon);*/
	  	 store_mat((*fwiTE).grad_sigma,gt_sigma,NX,NY);
	  	 store_mat((*fwiTE).grad_epsilon,gt_epsilon,NX,NY);

	  	 /* copy sigma_old -> sigmar */
	  	 /* copy epsilon_old -> epsilonr */
	  	 store_mat((*fwiTE).sigma_old,(*matTE).sigmar,NX,NY);
	  	 store_mat((*fwiTE).epsilon_old,(*matTE).epsilonr,NX,NY);

                 lsiter++;

              }else{

                 alpha = 0.0;
                 break;

              }
           
              gts0 = dotp_matrix(gt_sigma,(*fwiTE).Hgrad_sigma,NX,NY);
              gts0+= dotp_matrix(gt_epsilon,(*fwiTE).Hgrad_epsilon,NX,NY);

              g0s0 = dotp_matrix((*fwiTE).grad_sigma,(*fwiTE).Hgrad_sigma,NX,NY);
              g0s0+= dotp_matrix((*fwiTE).grad_epsilon,(*fwiTE).Hgrad_epsilon,NX,NY);

              if(MYID==0){
                 printf("Wolfe Conditions \n");
                 printf("---------------- \n");
                 printf("ft = %e \t <= L2  = %e \n",ft,L2);
                 printf("ft = %e \t <= L2 + C1*alpha*g0s0 = %e \n",ft,L2 + C1*alpha*g0s0);
                 printf("gts0 = %e \t >= C2*g0s0 = %e \n",gts0,C2*g0s0);
              }
 
              if (ft > (L2 + C1*alpha*g0s0)){
                  nu = alpha;
	          alpha = (nu + mu)/SCALEFAC;
              }else if(gts0 < (C2*g0s0)){
                  mu = alpha;
		  if(nu==1e30){
                     alpha = SCALEFAC*alpha;
                  }else{
		     alpha = (nu + mu)/SCALEFAC;
		  }
              }else{
                  done = 1;
              }

        }

        if(MYID==0){
           printf("final alpha = %e \n",alpha);
        }

	/* convert sigmar -> sigma and epsilonr -> epsilon to perform forward modelling */
	scale_grad((*matTE).sigmar,MAT1_NORM,(*matTE).sigma,NX,NY);
	scale_grad((*matTE).epsilonr,MAT2_NORM,(*matTE).epsilon,NX,NY);

        free_matrix(gt_sigma,1,NY,1,NX);
        free_matrix(gt_epsilon,1,NY,1,NX);

return alpha;
                	    
}
