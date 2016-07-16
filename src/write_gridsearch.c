/*------------------------------------------------------------------------
 *  Write grid search results                           
 *  
 *  D. Koehn
 *  Kiel, 19.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void write_gridsearch(float L2, float vp0, float grad0, int count){

	/* global variables */
        extern char GRIDSEARCH_FILE[STRING_SIZE];      

	/* local variables */
	int i;
	
        FILE *fp;

        if(count==1){
	   fp=fopen(GRIDSEARCH_FILE,"w");
	   fprintf(fp,"Objfcn \t         vp0 \t         grad0 \n");
        }else{
           fp=fopen(GRIDSEARCH_FILE,"a");
        }

        fprintf(fp,"%e \t %e \t %e \n",L2,vp0,grad0);

        fclose(fp);

}
