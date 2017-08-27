/*------------------------------------------------------------
 * Contribution of Tikhonov regularization term to gradient
 * 
 * Daniel Koehn
 * Kiel, 26.08.2017
 * -----------------------------------------------------------*/

#include "fd.h"

void Tikhonov_grad_TE(struct fwiTE *fwiTE, struct matTE *matTE, int iter){

	extern int NX, NY, IDX, IDY;
	extern float MAT1_NORM, MAT2_NORM, BETA_MAT1, BETA_MAT2;

	int i, j, h;
	float laplace_sigma, laplace_epsilon;

	/* calculate normalized material parameters */
	scale_grad((*matTE).sigma,1.0/MAT1_NORM,(*matTE).sigmar,NX,NY);
	scale_grad((*matTE).epsilon,1.0/MAT2_NORM,(*matTE).epsilonr,NX,NY);

	/* calculate Laplacian of sigma and epsilon and update gradient */
	for (i=2;i<=NX-1;i=i+IDX){
	   for (j=2;j<=NY-1;j=j+IDY){

	       laplace_sigma = (*matTE).sigmar[j][i-1] + (*matTE).sigmar[j][i+1] + (*matTE).sigmar[j-1][i] + (*matTE).sigmar[j+1][i] + 4.0 * (*matTE).sigmar[j][i];
	       laplace_epsilon = (*matTE).epsilonr[j][i-1] + (*matTE).epsilonr[j][i+1] + (*matTE).epsilonr[j-1][i] + (*matTE).epsilonr[j+1][i] + 4.0 * (*matTE).epsilonr[j][i];

	       (*fwiTE).grad_sigma[j][i] += (*fwiTE).lambda_1 * BETA_MAT1 * laplace_sigma;
	       (*fwiTE).grad_epsilon[j][i] += (*fwiTE).lambda_2 * BETA_MAT2 * laplace_epsilon;

	   }
	}

}
