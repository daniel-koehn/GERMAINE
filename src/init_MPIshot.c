/*
 * Initiate MPI parallelization by distributing shots over MPI processes
 *
 * Daniel Koehn
 * Kiel, 11.12.2015
 */

#include "fd.h"

void init_MPIshot(int nshots){

        /* global variables */
	extern int NP, MYID, NSHOT1, NSHOT2;

	/* local variables */
	int i, j;
	
	/* distribute shots over MPI processes */
	NSHOT1 = (nshots / NP) * MYID;

	if (nshots % NP > MYID){
	  NSHOT1 += MYID;
	  NSHOT2 = NSHOT1 + (nshots / NP) + 1;
	}else{
	  NSHOT1 += nshots % NP;
	  NSHOT2 = NSHOT1 + (nshots / NP);
	}

	NSHOT1++;
	NSHOT2++;

	/*printf("MYID = %d \t NSHOT1 = %d \t NSHOT2 = %d \n",MYID,NSHOT1,NSHOT2);*/

}



