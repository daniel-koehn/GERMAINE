/*------------------------------------------------------------------------
 *  Read P-wave velocity model
 *
 *  D. Koehn
 *  Kiel, 7.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readmod(struct matAC *matAC){

	extern int NX, NY, NX0, NY0, MYID, NPML, FREE_SURF;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	int i, j, yshift;
	float vp, rho;
	FILE *fp_vp, *fp_rho;
	char filename[STRING_SIZE];

	if(FREE_SURF==1){yshift=0;}
	if(FREE_SURF==0){yshift=NPML;}

	/* read density and seismic velocities */
	/* ----------------------------------- */

	fprintf(FP,"\n... reading P-wave velocity model from file...\n");           

	fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
	sprintf(filename,"%s.vp",MFILE);
	fp_vp=fopen(filename,"r");
	if (fp_vp==NULL) err(" Could not open model file for Vp ! ");           

	fprintf(FP,"\n... reading density model from file...\n");           

	fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	sprintf(filename,"%s.rho",MFILE);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) err(" Could not open model file for density ! ");
	   
	/* loop over global grid */
	for (i=1;i<=NX0;i++){
		for (j=1;j<=NY0;j++){

			fread(&vp, sizeof(float), 1, fp_vp);
			(*matAC).vp[j+yshift][i+NPML] = vp;				

			fread(&rho, sizeof(float), 1, fp_rho);
			(*matAC).rho[j+yshift][i+NPML] = rho;				

		}
	}

	fclose(fp_vp);
	fclose(fp_rho);

	/* extend model */

	/* left boundary */
	for (i=1;i<=NPML;i++){
		for (j=1;j<=NY;j++){
			(*matAC).vp[j][i] = (*matAC).vp[j][NPML+1];
			(*matAC).rho[j][i] = (*matAC).rho[j][NPML+1];				
		}
	}
	
	/* right boundary */
	for (i = NPML + NX0 + 1;i <= NX;i++){
		for (j=1;j<=NY;j++){
			(*matAC).vp[j][i] = (*matAC).vp[j][NPML + NX0];				
			(*matAC).rho[j][i] = (*matAC).rho[j][NPML + NX0];				
		}
	}

	/* top boundary */
	if(FREE_SURF==0){
	    for (i=1;i<=NX;i++){
	        for (j=1;j<=NPML;j++){
		    (*matAC).vp[j][i] = (*matAC).vp[NPML+1][i];				
		    (*matAC).rho[j][i] = (*matAC).rho[NPML+1][i];				
		}
	    }
	}

	/* bottom boundary */
	for (i=1;i<=NX;i++){
		for (j = NY - NPML + 1;j <= NY;j++){
			(*matAC).vp[j][i] = (*matAC).vp[NY-NPML][i];				
			(*matAC).rho[j][i] = (*matAC).rho[NY-NPML][i];				
		}
	}	

	/* each PE writes his model to disk */
	sprintf(filename,"%s.germaine.vp",MFILE);
	writemod(filename,(*matAC).vp,3);

	sprintf(filename,"%s.germaine.rho",MFILE);
	writemod(filename,(*matAC).rho,3);


}




