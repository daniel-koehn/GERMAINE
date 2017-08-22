/*------------------------------------------------------------------------
 * Module for extraction of Preconditioned Conjugate Gradient Method (PCG)
 * for the TE-mode case
 * 
 * Daniel Koehn
 * Kiel, 21.08.2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

void extract_PCG_TE(float * PCG_old, struct fwiTE *fwiTE){

	extern int NX, NY, IDX, IDY, MYID;
	extern char JACOBIAN[STRING_SIZE];

	int i, j, h;
	float tmp;
	char jac[225];

	FILE *FP;
	
        h=1;
	/* extract sigma gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 (*fwiTE).Hgrad_sigma[j][i] = PCG_old[h];

                 h++;
	   }
	}

	/* extract epsilon gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 (*fwiTE).Hgrad_epsilon[j][i] = PCG_old[h];

                 h++;
	   }
	}

	if(MYID==0){

	   /* save sigma gradient */
	   sprintf(jac,"%s_sigma_p.old",JACOBIAN);
	   FP=fopen(jac,"wb");

           for (i=1;i<=NX;i=i+IDX){
              for (j=1;j<=NY;j=j+IDY){
		 tmp = (*fwiTE).Hgrad_sigma[j][i];
                 fwrite(&tmp,sizeof(float),1,FP);
              }
           }
	
	   fclose(FP);

	   /* save epsilon gradient */
	   sprintf(jac,"%s_epsilon_p.old",JACOBIAN);
	   FP=fopen(jac,"wb");

           for (i=1;i<=NX;i=i+IDX){
              for (j=1;j<=NY;j=j+IDY){
		 tmp = (*fwiTE).Hgrad_epsilon[j][i];
                 fwrite(&tmp,sizeof(float),1,FP);
              }
           }
	
	   fclose(FP);


	}


}
