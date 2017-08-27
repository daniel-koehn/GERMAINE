/*------------------------------------------------------------------------
 * Contribution of Tikhonov regularization term to cost function
 * 
 * Daniel Koehn
 * Kiel, 26.08.2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

float Tikhonov_cost_TE(struct fwiTE *fwiTE, struct matTE *matTE, float L2, int iter){

	extern int NX, NY, IDX, IDY, EST_HYPER, MYID;
	extern float MAT1_NORM, MAT2_NORM, BETA_MAT1, BETA_MAT2;

	int i, j, h;
	float laplace_sigma, laplace_epsilon, tikh_sigma, tikh_epsilon;
	float l2_tikhonov;


	/* calculate normalized material parameters */
	scale_grad((*matTE).sigma,1.0/MAT1_NORM,(*matTE).sigmar,NX,NY);
	scale_grad((*matTE).epsilon,1.0/MAT2_NORM,(*matTE).epsilonr,NX,NY);

	tikh_sigma = 0.0;
	tikh_epsilon = 0.0;
	/* calculate Laplacian of sigma and epsilon */
	for (i=2;i<=NX-1;i=i+IDX){
	   for (j=2;j<=NY-1;j=j+IDY){

	       laplace_sigma = (*matTE).sigmar[j][i-1] + (*matTE).sigmar[j][i+1] + (*matTE).sigmar[j-1][i] + (*matTE).sigmar[j+1][i] + 4.0 * (*matTE).sigmar[j][i];
	       laplace_epsilon = (*matTE).epsilonr[j][i-1] + (*matTE).epsilonr[j][i+1] + (*matTE).epsilonr[j-1][i] + (*matTE).epsilonr[j+1][i] + 4.0 * (*matTE).epsilonr[j][i];

	       tikh_sigma += 0.5 * BETA_MAT1 * (*matTE).sigmar[j][i] * laplace_sigma;
	       tikh_epsilon += 0.5 * BETA_MAT2 * (*matTE).epsilonr[j][i] * laplace_epsilon;

	   }
	}

	/* estimate hyperparameters (*fwiTE).lambda_1 and (*fwiTE).lambda_2 for initial model */
	if(EST_HYPER==1){

	   /* sigma */
	   if(fabs(tikh_sigma)>0.0){
	       (*fwiTE).lambda_1 = L2 / tikh_sigma;
	   }else{
	       (*fwiTE).lambda_1 = 0.0;	       
	   }

	   /* epsilon */
	   if(fabs(tikh_epsilon)>0.0){
	       (*fwiTE).lambda_2 = L2 / tikh_epsilon;
	   }else{
	       (*fwiTE).lambda_2 = 0.0;	       
	   }

	   EST_HYPER=0;

	   if(MYID==0){
	      printf("Hyperparameters: lambda_1 = %e \t lambda_2 = %e \n",(*fwiTE).lambda_1,(*fwiTE).lambda_2);
	   }

	}

	/* update objective function */
	l2_tikhonov = L2 + (*fwiTE).lambda_1 * tikh_sigma + (*fwiTE).lambda_2 * tikh_epsilon;

	return l2_tikhonov;
}
