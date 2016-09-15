/*------------------------------------------------------------------------
 *  Assemble 1D linear gradient model for grid search
 *
 *  D. Koehn
 *  Kiel, 19.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void model(float  ** Vp){

        /* global variables */
	extern int NX, NY;
	extern float DH;
        extern char  MFILE[STRING_SIZE];

	/* local variables */
	int i, j;
        float vp0 = 1370.0, grad0 = 0.9, vptrue, vptmp;
        char filename[STRING_SIZE];
	
	FILE *fp_vp;

        /*printf("vp0 = %e \t grad0 = %e \n",vp0,grad0);*/
	
	sprintf(filename,"%s.vp",MFILE);
	fp_vp=fopen(filename,"r");
	if (fp_vp==NULL) err(" Could not open model file for Vp ! ");

	/* loop over global grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){

		     /* 1D linear gradient model */
                     Vp[j][i] = grad0 * ((j-1)*DH) + vp0;
		     
		     /* add water layer from true model */
                     fread(&vptrue, sizeof(float), 1, fp_vp);
                     if((i==1)&&(j==1)){vptmp = vptrue;}
		     if(vptrue==vptmp){
		        Vp[j][i] = vptmp;
		     }
                     
		}
	}

        fclose(fp_vp);
   
        /* output of model */
	sprintf(filename,"%s.rajzel.vp",MFILE);
	writemod(filename,Vp,3);
	MPI_Barrier(MPI_COMM_WORLD);
	
}
