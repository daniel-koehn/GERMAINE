/*------------------------------------------------------------------------
 *  Average density
 *
 *  D. Koehn
 *  Kiel, 15.08.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void av_rho(struct matAC *matAC, int i, int j){	

	/*(*matAC).bpm = 0.25 * ((*matAC).b[j][i] + (*matAC).b[j][i+1] + (*matAC).b[j-1][i] + (*matAC).b[j-1][i+1]);

	(*matAC).bmp = 0.25 * ((*matAC).b[j][i] + (*matAC).b[j][i-1] + (*matAC).b[j+1][i] + (*matAC).b[j+1][i-1]);

	(*matAC).bpp = 0.25 * ((*matAC).b[j][i] + (*matAC).b[j][i+1] + (*matAC).b[j+1][i] + (*matAC).b[j+1][i+1]);

	(*matAC).bmm = 0.25 * ((*matAC).b[j][i] + (*matAC).b[j][i-1] + (*matAC).b[j-1][i] + (*matAC).b[j-1][i-1]);

	(*matAC).b00 = (*matAC).b[j][i];

	(*matAC).bp0 = 0.5 * ((*matAC).b[j][i] + (*matAC).b[j][i+1]);

	(*matAC).bm0 = 0.5 * ((*matAC).b[j][i] + (*matAC).b[j][i-1]);

	(*matAC).b0p = 0.5 * ((*matAC).b[j][i] + (*matAC).b[j+1][i]);

	(*matAC).b0m = 0.5 * ((*matAC).b[j][i] + (*matAC).b[j-1][i]);*/

	/* printf("bpm = %f, bmp = %f, bpp = %f, bmm = %f, b00 = %f \n",(*matAC).bpm,(*matAC).bmp,(*matAC).bpp,(*matAC).bmm,(*matAC).b00); */

	(*matAC).bpm = (*matAC).b[j][i];

	(*matAC).bmp = (*matAC).b[j][i];

	(*matAC).bpp = (*matAC).b[j][i];

	(*matAC).bmm = (*matAC).b[j][i];

	(*matAC).b00 = (*matAC).b[j][i];

	(*matAC).bp0 = (*matAC).b[j][i];

	(*matAC).bm0 = (*matAC).b[j][i];

	(*matAC).b0p = (*matAC).b[j][i];

	(*matAC).b0m = (*matAC).b[j][i];

}




