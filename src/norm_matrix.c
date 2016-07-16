/* Euclidean norm of a matrix
 *
 * Daniel Koehn
 * Kiel, the 14.12.2015
 */

#include "fd.h"

float norm_matrix(float **A, int NX, int NY){

	float sum;
	int i, j;

        sum = 0.0;

	for (i=1;i<=NX;i++){
	   for (j=1;j<=NY;j++){
		  sum += A[j][i]*A[j][i];
	   }
	}

        sum = sqrt(sum);
   
	return sum;	
}
