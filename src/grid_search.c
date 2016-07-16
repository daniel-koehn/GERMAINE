/*------------------------------------------------------------------------
 *  1D linear gradient model by grid search
 *
 *  D. Koehn
 *  Kiel, 19.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void grid_search(float ** Vp, float ** S, float ** TT, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr){

	/* declaration of global variables */
        extern int SEISMO, NX, NY, MYID, NSHOT1, NSHOT2;
        extern float VP0_1, VP0_2, DVP0, GRAD0_1, GRAD0_2, DGRAD0;
    
        /* declaration of local variables */
        int ishot, ngrad0, nvp0, i, j, nobj, count; 
        float L2sum, L2, vp0, grad0;

	if (MYID==0){
	   printf("\n\n\n ------------------------------------------------------------------\n");
	   printf("\n\n\n                1D linear gradient grid search                     \n");
	   printf("\n\n\n ------------------------------------------------------------------\n");
	}

        /* suppress modelled FA traveltime output during FATT */
        SEISMO = 0;

        /* calculate number of objective function evaluations */
        nvp0 = ceil((VP0_2-VP0_1)/DVP0);
        ngrad0 = ceil((GRAD0_2-GRAD0_1)/DGRAD0);
        nobj = nvp0 * ngrad0;


        if (MYID==0) printf("nvp0 =  %d, ngrad0 = %d\n",nvp0,ngrad0);

        /* loop over parameter space */
        count = 1;
        vp0 = VP0_1;
        for (i=1;i<=nvp0;i++){

                grad0 = GRAD0_1;
		for (j=1;j<=ngrad0;j++){

                        if (MYID==0){
	                    printf("Model %d of %d \n",count,nobj);
	                }

                        /* assemble new model and calculate slowness */
                        model_gridsearch(Vp,vp0,grad0);
                        // calc_S(Vp,S);

			/* evaluate objective function */
			L2 = 0.0;

                        /*printf("MYID = %d \t NSHOT1 = %d \t NSHOT2 = %d \n",MYID,NSHOT1,NSHOT2);*/

			for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

			     /* solve forward problem with Eikonal solver */
			     // eikonal(S,TT,Tmod,srcpos,ishot,nshots,recpos,ntr);

			     /* calculate traveltime residuals at receiver positions */
		 	     // L2+=calc_FA_res(Tmod,Tobs,Tres,ntr,ishot);
				

			} /* end of loop over shots */   


			/* assemble objective function from all MPI processes */
			L2sum = 0.0;

                        /*printf("Before MPI_Allreduce MYID = %d \t L2 = %e \n",MYID,L2);*/
                        
                        MPI_Barrier(MPI_COMM_WORLD);

			MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
                        L2 = 0.5 * L2sum;

                        /*printf("After MPI_Allreduce MYID = %d \t L2 = %e \n",MYID,L2);*/

                        /* output of objective function value and model parameters */
                        if(MYID==0){write_gridsearch(L2,vp0,grad0,count);}
                        
                        grad0+=DGRAD0;  /* update gradient */
                        count++;
                }
                
                vp0+=DVP0;  /* update P-wave velocity @ surface*/
        }	
                	    
}
