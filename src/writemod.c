/*------------------------------------------------------------------------
 *   write local model to file              
 *   last update 16/02/02   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void writemod(char modfile[STRING_SIZE], float ** array, int format){


	/* extern variables */
	extern int MYID, NX, NY, IDX, IDY;
	extern FILE *FP;


	int i, j;
	FILE *fpmod;
	char file[STRING_SIZE];

	printf("\n\n PE %d is writing model to \n",MYID);
	sprintf(file,"%s",modfile);
	printf("\t%s\n\n", file);
	fpmod=fopen(file,"w");
	for (i=1;i<=NX;i+=IDX)
	for (j=1;j<=NY;j+=IDY)
		writedsk(fpmod,array[j][i],format);
				
	fclose(fpmod);


}


