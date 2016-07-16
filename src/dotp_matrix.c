/* Calculate dot product of two matrices
 *
 *
 * Daniel Koehn
 * Kiel, the 14.12.2015
 */

#include "fd.h"

float dotp_matrix(float ** A, float ** B, int NX, int NY){

	float sum;
	int i,j;

        sum = 0.0;        
	for (i=1;i<=NX;i++){
	   for (j=1;j<=NY;j++){
		  sum += A[j][i]*B[j][i];
	   }
	}       

	return sum;	
}


