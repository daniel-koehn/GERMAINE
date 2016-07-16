/*------------------------------------------------------------------------
 *   write local model to file              
 *   last update 16/02/02   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void writemod_vec(char modfile[STRING_SIZE], double * array, int format){


	/* extern variables */
	extern int MYID, NX, NY, IDX, IDY;
	extern FILE *FP;

	float amp; 
	int i, j, h;
	FILE *fpmod;
	char file[STRING_SIZE];

	fprintf(FP,"\n\n PE %d is writing model to \n",MYID);
	sprintf(file,"%s",modfile);
	fprintf(FP,"\t%s\n\n", file);
	fpmod=fopen(file,"wb");
        h=0;
	for (i=1;i<=NX;i+=IDX){
	    for (j=1;j<=NY;j+=IDY){
                amp = (float)(array[h]);
		fwrite(&amp, sizeof(float), 1, fpmod);
		h++;
	    }
     	}
				
	fclose(fpmod);


}


