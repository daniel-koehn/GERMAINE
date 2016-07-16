/*------------------------------------------------------------------------
 *  Write FD seismograms to file
 *  
 *  D. Koehn
 *  Kiel, 21/06/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void write_seis_AC(struct waveAC *waveAC, int ishot, int ntr, int nstage, int nfreq){

	/* global variables */
        extern char PICKS_FILE[STRING_SIZE];        
        extern int MYID, INFO, NONZERO;

	/* local variables */
	int i;
        char pickfile_char[STRING_SIZE];
	float tmp, tmp1;
	
        FILE *fp;

        if((MYID==0)&&(INFO==1)){
	    printf("\n Write FD seismograms for shot %d stage %d frequency %d ...\n",ishot,nstage,nfreq);
        }

        sprintf(pickfile_char,"%s_p_shot_%d_stage_%d_nfreq_%d.bin",PICKS_FILE,ishot,nstage, nfreq);

        fp=fopen(pickfile_char,"wb");

	for(i=1;i<=ntr;i++){
	    // fprintf(fp,"%e \t %e \n",(*waveAC).precr[i],(*waveAC).preci[i]); 
            tmp = (*waveAC).precr[i]; 
	    tmp1 = (*waveAC).preci[i];
 
	    fwrite(&tmp, sizeof(float), 1, fp);
	    fwrite(&tmp1, sizeof(float), 1, fp);
	}

        fclose(fp);
	
}
