/*
 * Summation of gradients over all MPI processes
 *
 * Daniel Koehn
 * Kiel, 08.01.2016
 */

#include "fd.h"

void exchange_grad_MPI(float ** grad){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j, h;
	float * grad_vec;
	
	grad_vec = (float *)malloc((NX+1)*(NY+1) * sizeof(float));
	
	/* store gradient in vector */
	h=1;
	for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){

               grad_vec[h] = grad[j][i];
	       h++;

            }
        }
	
	MPI_Barrier(MPI_COMM_WORLD);
	
        MPI_Bcast(grad_vec,NX*NY,MPI_FLOAT,0,MPI_COMM_WORLD);

        /* store gradient vector in matrix again */
	h=1;
	for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){

               grad[j][i] = grad_vec[h];
	       h++;

            }
        }

	free( grad_vec );
	
}



