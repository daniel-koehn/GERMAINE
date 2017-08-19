/*------------------------------------------------------------------------
 *  Read sigma and epsilon model
 *
 *  D. Koehn
 *  Kiel, 19.08.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readmod_TE(struct matTE *matTE){

	extern int NX, NY, NX0, NY0, MYID, NPML, FREE_SURF;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	int i, j, yshift;
	float sigma, epsilon;
	FILE *fp_sig, *fp_eps;
	char filename[STRING_SIZE];

	if(FREE_SURF==1){yshift=0;}
	if(FREE_SURF==0){yshift=NPML;}

	/* read conductivity and permittivity models */
	/* ----------------------------------------- */

	fprintf(FP,"\n... reading conductivity model from file...\n");           

	fprintf(FP,"\t Sigma:\n\t %s.sig\n\n",MFILE);
	sprintf(filename,"%s.sig",MFILE);
	fp_sig=fopen(filename,"r");
	if (fp_sig==NULL) err(" Could not open model file for Sigma ! ");           

	fprintf(FP,"\n... reading permittivity model from file...\n");           

	fprintf(FP,"\t Epsilon:\n\t %s.rho\n\n",MFILE);
	sprintf(filename,"%s.eps",MFILE);
	fp_eps=fopen(filename,"r");
	if (fp_eps==NULL) err(" Could not open model file for Epsilon ! ");
	   
	/* loop over global grid */
	for (i=1;i<=NX0;i++){
		for (j=1;j<=NY0;j++){

			fread(&sigma, sizeof(float), 1, fp_sig);
			(*matTE).sigma[j+yshift][i+NPML] = sigma;				

			fread(&epsilon, sizeof(float), 1, fp_eps);
			(*matTE).epsilon[j+yshift][i+NPML] = epsilon;				

		}
	}

	fclose(fp_sig);
	fclose(fp_eps);

	/* extend model */

	/* left boundary */
	for (i=1;i<=NPML;i++){
		for (j=1;j<=NY;j++){
			(*matTE).sigma[j][i] = (*matTE).sigma[j][NPML+1];
			(*matTE).epsilon[j][i] = (*matTE).epsilon[j][NPML+1];				
		}
	}
	
	/* right boundary */
	for (i = NPML + NX0 + 1;i <= NX;i++){
		for (j=1;j<=NY;j++){
			(*matTE).sigma[j][i] = (*matTE).sigma[j][NPML + NX0];				
			(*matTE).epsilon[j][i] = (*matTE).epsilon[j][NPML + NX0];				
		}
	}

	/* top boundary */
	if(FREE_SURF==0){
	    for (i=1;i<=NX;i++){
	        for (j=1;j<=NPML;j++){
		    (*matTE).sigma[j][i] = (*matTE).sigma[NPML+1][i];				
		    (*matTE).epsilon[j][i] = (*matTE).epsilon[NPML+1][i];				
		}
	    }
	}

	/* bottom boundary */
	for (i=1;i<=NX;i++){
		for (j = NY - NPML + 1;j <= NY;j++){
			(*matTE).sigma[j][i] = (*matTE).sigma[NY-NPML][i];				
			(*matTE).epsilon[j][i] = (*matTE).epsilon[NY-NPML][i];				
		}
	}	

	/* each PE writes his model to disk */
	sprintf(filename,"%s.germaine.sig",MFILE);
	writemod(filename,(*matTE).sigma,3);

	sprintf(filename,"%s.germaine.eps",MFILE);
	writemod(filename,(*matTE).epsilon,3);


}




