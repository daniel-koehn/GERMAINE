/*
 *  Taper gradient at source and receiver positions
 * 
 *  D. Koehn
 *  Kiel, 13.12.2015
 *
 */

#include "fd.h"

void taper_grad_shot1(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int ishot)
{

	/* gloabal variables */
        extern float DH;
	extern FILE *FP;
	
	/* local variables */
	int i, j, n;
	int i1, j1;
        float damp_src_rec=1e-5;

	char modfile[STRING_SIZE];
	
	extern int FILTSIZE;

        /* *********************************** */
        /* Taper source and receiver positions */
        /* *********************************** */
                
	/* set gradient near source positions to zero */
	for (n=1;n<=nshots;n++){
	
	     i1 = iround(srcpos[1][n]/DH);
	     j1 = iround(srcpos[2][n]/DH);
		
             for (i=i1-FILTSIZE;i<=i1+FILTSIZE;i++){
		  for (j=j1-FILTSIZE;j<=j1+FILTSIZE;j++){
		
		       waveconv[j][i] = damp_src_rec;
			
		       
		  }
	     }	       
	}
	
	
	/* set gradient near receiver positions to zero */
	for (n=1;n<=ntr;n++){
	
	     i1 = recpos[1][n];
	     j1 = recpos[2][n];
		
             for (i=i1-FILTSIZE;i<=i1+FILTSIZE;i++){
		  for (j=j1-FILTSIZE;j<=j1+FILTSIZE;j++){
		
		       waveconv[j][i] = damp_src_rec;
			
		       
		  }
	     }	       
	}	
	
	/*MPI_Barrier(MPI_COMM_WORLD);
        sprintf(modfile,"taper_coeff_%i.bin",ishot);
        writemod(modfile,taper_coeff,3); */
	
}



