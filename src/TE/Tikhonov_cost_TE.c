/*------------------------------------------------------------------------
 * Contribution of Tikhonov regularization term to cost function
 * 
 * Daniel Koehn
 * Kiel, 26.08.2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

float Tikhonov_cost_TE(struct fwiTE *fwiTE, struct matTE *matTE, float L2, int iter){

	extern int NX, NY, IDX, IDY, MYID;
	extern float MAT1_NORM, MAT2_NORM, BETA_MAT1, BETA_MAT2, DH;
	extern float LAMBDA_1, LAMBDA_2;

	int i, j, h;
	float laplace_sigma, laplace_epsilon, tikh_sigma, tikh_epsilon;
	float l2_tikhonov;


	/* calculate normalized material parameters */
	scale_grad((*matTE).sigma,BETA_MAT1/MAT1_NORM,(*matTE).sigmar,NX,NY);
	scale_grad((*matTE).epsilon,BETA_MAT2/MAT2_NORM,(*matTE).epsilonr,NX,NY);

	tikh_sigma = 0.0;
	tikh_epsilon = 0.0;
	/* calculate Laplacian of sigma and epsilon */
	for (i=2;i<=NX-1;i=i+IDX){
	   for (j=2;j<=NY-1;j=j+IDY){

	       laplace_sigma = ((*matTE).sigmar[j][i-1] + (*matTE).sigmar[j][i+1] + (*matTE).sigmar[j-1][i] + (*matTE).sigmar[j+1][i] - 4.0 * (*matTE).sigmar[j][i]) / (DH*DH);
	       laplace_epsilon = ((*matTE).epsilonr[j][i-1] + (*matTE).epsilonr[j][i+1] + (*matTE).epsilonr[j-1][i] + (*matTE).epsilonr[j+1][i] - 4.0 * (*matTE).epsilonr[j][i])  / (DH*DH);

	       tikh_sigma += 0.5 * BETA_MAT1 * laplace_sigma * laplace_sigma;
	       tikh_epsilon += 0.5 * BETA_MAT2 * laplace_epsilon * laplace_epsilon;

	   }
	}

	/* update objective function */
	if((LAMBDA_1>0.0)||(LAMBDA_2>0.0)){
	   l2_tikhonov = L2 + LAMBDA_1 * tikh_sigma + LAMBDA_2 * tikh_epsilon;
	}else{
           l2_tikhonov = L2;
	}

	return l2_tikhonov;
}
