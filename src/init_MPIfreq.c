/*
 * Initiate MPI parallelization by distributing frequencies over MPI processes
 *
 * Daniel Koehn
 * Kiel, 25.06.2017
 */

#include "fd.h"

void init_MPIfreq(){

        /* global variables */
	extern int NPROCFREQ, MYID_SHOT, NFREQ1, NFREQ2, NF;
	extern int MYID, NP, MYID_SHOT, COLOR;

	/* local variables */
	int i, j;
	
	/* distribute frequencies over MPI processes */
	NFREQ1 = (NF / NPROCFREQ) * MYID_SHOT;

	if (NF % NPROCFREQ > MYID_SHOT){
	  NFREQ1 += MYID_SHOT;
	  NFREQ2 = NFREQ1 + (NF / NPROCFREQ) + 1;
	}else{
	  NFREQ1 += NF % NPROCFREQ;
	  NFREQ2 = NFREQ1 + (NF / NPROCFREQ);
	}

	NFREQ1++;
	NFREQ2++;

	/* printf("WORLD RANK/SIZE: %d/%d \t shot_comm COLOR/RANK/SIZE/NFREQ1/NFREQ2: %d/%d/%d/%d/%d\n", MYID, NP, COLOR, MYID_SHOT, NPROCFREQ, NFREQ1, NFREQ2); */

}



