/*------------------------------------------------------------------------
 *  Average shear modulus
 *
 *  D. Koehn
 *  Kiel, 15.08.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void av_mu(struct matSH *matSH, int i, int j){	

	(*matSH).bpm = 0.25 * ((*matSH).mu[j][i] + (*matSH).mu[j][i+1] + (*matSH).mu[j-1][i] + (*matSH).mu[j-1][i+1]);

	(*matSH).bmp = 0.25 * ((*matSH).mu[j][i] + (*matSH).mu[j][i-1] + (*matSH).mu[j+1][i] + (*matSH).mu[j+1][i-1]);

	(*matSH).bpp = 0.25 * ((*matSH).mu[j][i] + (*matSH).mu[j][i+1] + (*matSH).mu[j+1][i] + (*matSH).mu[j+1][i+1]);

	(*matSH).bmm = 0.25 * ((*matSH).mu[j][i] + (*matSH).mu[j][i-1] + (*matSH).mu[j-1][i] + (*matSH).mu[j-1][i-1]);

	(*matSH).b00 = (*matSH).mu[j][i];

	(*matSH).bp0 = 0.5 * ((*matSH).mu[j][i] + (*matSH).mu[j][i+1]);

	(*matSH).bm0 = 0.5 * ((*matSH).mu[j][i] + (*matSH).mu[j][i-1]);

	(*matSH).b0p = 0.5 * ((*matSH).mu[j][i] + (*matSH).mu[j+1][i]);

	(*matSH).b0m = 0.5 * ((*matSH).mu[j][i] + (*matSH).mu[j-1][i]);

	/* printf("bpm = %f, bmp = %f, bpp = %f, bmm = %f, b00 = %f \n",(*matSH).bpm,(*matSH).bmp,(*matSH).bpp,(*matSH).bmm,(*matSH).b00); */

	/*(*matSH).bpm = (*matSH).mu[j][i];

	(*matSH).bmp = (*matSH).mu[j][i];

	(*matSH).bpp = (*matSH).mu[j][i];

	(*matSH).bmm = (*matSH).mu[j][i];

	(*matSH).b00 = (*matSH).mu[j][i];

	(*matSH).bp0 = (*matSH).mu[j][i];

	(*matSH).bm0 = (*matSH).mu[j][i];

	(*matSH).b0p = (*matSH).mu[j][i];

	(*matSH).b0m = (*matSH).mu[j][i];*/

}




