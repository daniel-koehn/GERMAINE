/*
 * Initiate MPI parallelization by distributing shots over MPI processes
 *
 * Daniel Koehn
 * Kiel, 11.12.2015
 */

#include "fd.h"

void init_MPIshot(int nshots){

        /* global variables */
	extern int NPROCSHOT, COLOR, NSHOT1, NSHOT2;
	/* extern int MYID, NP, MYID_SHOT; */

	/* local variables */
	int i, j;
	
	/* distribute shots over MPI processes */
	NSHOT1 = (nshots / NPROCSHOT) * COLOR;

	if (nshots % NPROCSHOT > COLOR){
	  NSHOT1 += COLOR;
	  NSHOT2 = NSHOT1 + (nshots / NPROCSHOT) + 1;
	}else{
	  NSHOT1 += nshots % NPROCSHOT;
	  NSHOT2 = NSHOT1 + (nshots / NPROCSHOT);
	}

	NSHOT1++;
	NSHOT2++;

	/*printf("WORLD RANK/SIZE: %d/%d \t shot_comm COLOR/RANK/SIZE/NSHOT1/NSHOT2: %d/%d/%d/%d/%d\n", MYID, NP, COLOR, MYID_SHOT, NPROCSHOT, NSHOT1, NSHOT2);*/

}



