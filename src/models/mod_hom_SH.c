/*------------------------------------------------------------------------
 *  Homogeneous SH model
 *
 *  D. Koehn
 *  Kiel, 05.10.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void model_SH(struct matSH *matSH){

        /* global variables */
	extern int NX, NY, NX0, NY0, NPML;
	extern int FREE_SURF;
	extern float DH;
        extern char  MFILE[STRING_SIZE];

	/* local variables */
	int i, j, yshift;
        float vs = 200.0, rho = 1500.0;
        char filename[STRING_SIZE];

	if(FREE_SURF==1){yshift=0;}
	if(FREE_SURF==0){yshift=NPML;}

	/* loop over global grid */
	for (i=1;i<=NX0;i++){
		for (j=1;j<=NY0;j++){

			(*matSH).vs[j+yshift][i+NPML] = vs;
			(*matSH).rho[j+yshift][i+NPML] = rho;
				
		}
	}

	/* extend model */

	/* left boundary */
	for (i=1;i<=NPML;i++){
		for (j=1;j<=NY;j++){
			(*matSH).vs[j][i] = (*matSH).vs[j][NPML+1];
			(*matSH).rho[j][i] = (*matSH).rho[j][NPML+1];				
		}
	}
	
	/* right boundary */
	for (i = NPML + NX0 + 1;i <= NX;i++){
		for (j=1;j<=NY;j++){
			(*matSH).vs[j][i] = (*matSH).vs[j][NPML + NX0];
			(*matSH).rho[j][i] = (*matSH).rho[j][NPML + NX0];
		}
	}

	/* top boundary */
	if(FREE_SURF==0){
	    for (i=1;i<=NX;i++){
	        for (j=1;j<=NPML;j++){
		    (*matSH).vs[j][i] = (*matSH).vs[NPML+1][i];
		    (*matSH).rho[j][i] = (*matSH).rho[NPML+1][i];				
		}
	    }
	}

	/* bottom boundary */
	for (i=1;i<=NX;i++){
		for (j = NY - NPML + 1;j <= NY;j++){
		    (*matSH).vs[j][i] = (*matSH).vs[NY-NPML][i];
		    (*matSH).rho[j][i] = (*matSH).rho[NY-NPML][i];				
		}
	}	

	/* each PE writes his model to disk */
	sprintf(filename,"%s.germaine.vs",MFILE);
	writemod(filename,(*matSH).vs,3);

	sprintf(filename,"%s.germaine.rho",MFILE);
	writemod(filename,(*matSH).rho,3);

}




