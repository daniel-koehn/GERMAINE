/*------------------------------------------------------------------------
 *  Solve TE-mode forward problem with FDFD for all shots
 *
 *  D. Koehn
 *  Kiel, 19.08.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void forward_shot_TE(struct waveAC *waveAC, struct PML_AC *PML_AC, struct matTE *matTE, float ** srcpos, int nshots, int ** recpos, int ntr, int nstage, int nfreq){

	/* declaration of global variables */
        extern int NSHOT1, NSHOT2, NONZERO, NX, NY, NXNY, NF;
        extern int SNAP, SEISMO, MYID, INFO, N_STREAMER, READ_REC;
        extern float DH;
        extern char SNAP_FILE[STRING_SIZE];
	extern FILE * FP;
    
        /* declaration of local variables */
        int ishot, status, nxsrc, nysrc, i;
    	double *null = (double *) NULL ;

	int     *Ap, *Ai;
	double  *Ax, *Az, *xr, *xi; 
	double  time1, time2;
        char filename[STRING_SIZE];
        void *Symbolic, *Numeric;

	/* Allocate memory for compressed sparse column form and solution vector */
	Ap = malloc(sizeof(int)*(NXNY+1));
	Ai = malloc(sizeof(int)*NONZERO);
	Ax = malloc(sizeof(double)*NONZERO);
	Az = malloc(sizeof(double)*NONZERO);
	xr = malloc(sizeof(double)*NONZERO);
	xi = malloc(sizeof(double)*NONZERO);

        /* assemble acoustic impedance matrix */
        init_A_TE_9p_pml(PML_AC,matTE,waveAC);

        /* convert triplet to compressed sparse column format */
        status = umfpack_zi_triplet_to_col(NXNY,NXNY,NONZERO,(*waveAC).irow,(*waveAC).icol,(*waveAC).Ar,(*waveAC).Ai,Ap,Ai,Ax,Az,NULL);

	/* Here is something buggy (*waveAC).Ar != Ax and (*waveAC).Ai != Az.
           Therefore, set  Ax = (*waveAC).Ar and Az = (*waveAC).Ai */
	for (i=0;i<NONZERO;i++){
	     Ax[i] = (*waveAC).Ar[i];
	     Az[i] = (*waveAC).Ai[i];
        }

	if((MYID==0)&&(INFO==1)){
	  printf("\n==================================== \n");
	  printf("\n *****  LU factorization **********  \n");
	  printf("\n==================================== \n\n");
	  time1=MPI_Wtime(); 
	}

        /* symbolic factorization */
	status = umfpack_zi_symbolic(NXNY, NXNY, Ap, Ai, Ax, Az, &Symbolic, null, null);

        /* sparse LU decomposition */
	status = umfpack_zi_numeric(Ap, Ai, Ax, Az, Symbolic, &Numeric, null, null);
        umfpack_zi_free_symbolic (&Symbolic);

	if((MYID==0)&&(INFO==1)){
	  time2=MPI_Wtime();
	  printf("\n Finished after %4.2f s \n",time2-time1);
	}

	if((MYID==0)&&(INFO==1)){
	  printf("\n============================================================================================================= \n");
	  printf("\n *****  Solve TE-mode forward problem by FDFD for shot %d - %d (f = %4.2f MHz) on MPI process no. %d **********  \n",NSHOT1,NSHOT2-1,(*waveAC).freq/1e6, MYID);
	  printf("\n============================================================================================================= \n\n");				
	  time1=MPI_Wtime(); 
	}

        /* loop over all shots */
	for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

		/* read receiver positions from receiver files for each shot */
		if(READ_REC==1){
		    acq.recpos=receiver(FP, &ntr, 1);	                         
      		}

                /* define source vector RHS */
                RHS_source_AC(waveAC,srcpos,ishot);

                /* solve forward problem by forward and back substitution */
    		status = umfpack_zi_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xi, (*waveAC).RHSr, (*waveAC).RHSi, Numeric, null, null);

		/* convert vector xr/xi to pr/pi */
		vec2mat((*waveAC).pr,(*waveAC).pi,xr,xi);

		/* write real part of pressure wavefield to file */
		if(SNAP==1){
		   sprintf(filename,"%s_shot_%d.p",SNAP_FILE,ishot);
		   writemod_true(filename,(*waveAC).pr,3);
		}

		/* calculate FD seismogram files */
		if(SEISMO==1){
		   calc_seis_AC(waveAC,acq.recpos,ntr,ishot,nshots,nfreq);
		}

		if(READ_REC==1){
		   free_imatrix(acq.recpos,1,3,1,ntr);
		   ntr=0;
 		}

	}


	if((MYID==0)&&(INFO==1)){
	  time2=MPI_Wtime();
	  printf("\n Finished after %4.2f s \n",time2-time1);
	}
         	    
	/* free memory */
    	free(Ap); free(Ai); free(Ax); free(Az); free(xr); free(xi); 

	umfpack_zi_free_numeric (&Numeric);	

}
