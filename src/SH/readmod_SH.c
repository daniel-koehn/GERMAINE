/*------------------------------------------------------------------------
 *  Read S-wave velocity and density model
 *
 *  D. Koehn
 *  Kiel, 26.05.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readmod_SH(struct matSH *matSH){

	extern int NX, NY, NX0, NY0, MYID, NPML, FREE_SURF;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	int i, j, yshift;
	float vs, rho;
	FILE *fp_vs, *fp_rho;
	char filename[STRING_SIZE];

	if(FREE_SURF==1){yshift=0;}
	if(FREE_SURF==0){yshift=NPML;}

	fprintf(FP,"\n... reading S-wave velocity model from file...\n");
           
	/* read S-wave velocity model */
	/* -------------------------- */
	fprintf(FP,"\t Vs:\n\t %s.vp\n\n",MFILE);
	sprintf(filename,"%s.vs",MFILE);
	fp_vs=fopen(filename,"r");
	if (fp_vs==NULL) err(" Could not open model file for Vs ! ");           
	
	fprintf(FP,"\n... reading density model from file...\n");
           
	/* read density model */
	/* ------------------ */
	fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	sprintf(filename,"%s.rho",MFILE);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) err(" Could not open model file for Density ! ");           

   
	/* loop over global grid */
	for (i=1;i<=NX0;i++){
		for (j=1;j<=NY0;j++){

			fread(&vs, sizeof(float), 1, fp_vs);
			(*matSH).vs[j+yshift][i+NPML] = vs;

			fread(&rho, sizeof(float), 1, fp_rho);
			(*matSH).rho[j+yshift][i+NPML] = rho;
				
		}
	}

	fclose(fp_vs);
	fclose(fp_rho);

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




