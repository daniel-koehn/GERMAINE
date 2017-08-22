/*------------------------------------------------------------------------
 *  Read FD seismograms from file
 *  
 *  D. Koehn
 *  Kiel, 07/08/2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void read_seis_AC(struct waveAC *waveAC, int nshots, int ntr, int nstage){

	/* global variables */
        extern char DATA_DIR[STRING_SIZE];        
        extern int MYID, INFO, NF;

	/* local variables */
	int i, j, k, index;
        char pickfile_char[STRING_SIZE];
	float tmp, tmp1;	

        if(MYID==0){
	    printf("\n Read FD seismograms for stage %d ...\n",nstage);
        }
	

        FILE *fp;
        sprintf(pickfile_char,"%s_p_stage_%d.bin",DATA_DIR,nstage);
        fp=fopen(pickfile_char,"rb");

	for(k=1;k<=NF;k++){ 	
	    for(j=1;j<=nshots;j++){
	        for(i=1;i<=ntr;i++){ 

		    index =  i + ntr * (j-1) + ntr * nshots * (k-1);

	            fread(&tmp, sizeof(float), 1, fp);
	            fread(&tmp1, sizeof(float), 1, fp);

		    (*waveAC).pobsr[index] = tmp;
		    (*waveAC).pobsi[index] = tmp1; 

		}
	    }
	}

	fclose(fp);      	

}
