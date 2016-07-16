/*------------------------------------------------------------------------
 *   Output of memory requirements                          
 *   
 *   D. Koehn
 *   Kiel, 18.12.2015
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info_mem(FILE *fp, int NLBFGS_vec, int ntr){

	/* global variables */
	extern char LOG_FILE[STRING_SIZE];
	extern int MYID, LOG, NX, NY, NLBFGS, INVMAT;

        /* local variables */
        int fac1;
        float memeikonal, memlbfgs, memfatt, memtotal, fac_float;

	/*allocate memory for dynamic, static and buffer arrays */
	fac1=NX*NY;
	fac_float=sizeof(float)*pow(2.0,-20.0);

        memeikonal = (3.0 * fac1 + ntr) * fac_float;
        memlbfgs = (2.0 * NLBFGS_vec * (1 + NLBFGS) + 3.0 * NLBFGS) * fac_float;
        memfatt = (5.0 * fac1 + 2.0 *ntr) * fac_float; 

        if(INVMAT==0){
	   memtotal = memeikonal;
        }

        if(INVMAT==1){
	   memtotal = memfatt + memlbfgs + memeikonal;
        }

	fprintf(fp,"\n Memory requirements \n");
	fprintf(fp," ------------------- \n\n");
	fprintf(fp," Size of FD grid: NX=%d \t NY=%d\n",NX,NY);
	fprintf(fp," Each process is now trying to allocate memory for:\n");
	fprintf(fp," Eikonal solver:%23.2f MB\n", memeikonal);
        if(INVMAT==1){        
           fprintf(fp," FATT:%33.2f MB\n", memfatt);
           fprintf(fp," l-BFGS:%31.2f MB\n", memlbfgs);
        }
	/* fprintf(fp," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
	fprintf(fp," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0)); */
	fprintf(fp," ------------------------------------------------ \n");
	fprintf(fp," Total memory required: \t %6.2f MB.\n\n", memtotal);
	
	fprintf(fp," ... memory allocation was successfull.\n\n");
 
	/* if (LOG==1) fprintf(fp," standard output. \n");
	    else    fprintf(fp," %s.%i .\n",LOG_FILE,MYID); */

}
