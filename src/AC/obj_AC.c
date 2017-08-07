/*------------------------------------------------------------------------
 *  Calculate objective function
 *
 *  D. Koehn
 *  Kiel, 29.06.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float obj_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage){

	/* declaration of global variables */
        extern int SEISMO, NX, NY, NSHOT1, NSHOT2, INFO, NONZERO, NXNY;
        extern int SEISMO, MYID, INFO, NF, N_STREAMER, READ_REC;
        extern int NFREQ1, NFREQ2;
        extern float DH;
	extern char SNAP_FILE[STRING_SIZE];
	extern FILE *FP;
    
        /* declaration of local variables */
        int ishot, i, nfreq; 
        float L2sum, L2;
        int status, nxsrc, nysrc;
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

        /* suppress modelled seismogram output during FWI */
        SEISMO = 0;

	/* set objective function to zero */
        L2 = 0.0;

	/* loop over frequencies at each stage */
        for(nfreq=NFREQ1;nfreq<NFREQ2;nfreq++){

		/* set frequency on local MPI process */
		(*waveAC).freq = (*waveAC).stage_freq[nfreq];

		/* define PML damping profiles */
		pml_pro(PML_AC,waveAC);

		init_mat_AC(waveAC,matAC);

		/* assemble acoustic impedance matrix */
		init_A_AC_9p_pml(PML_AC,matAC,waveAC);	

		/* convert triplet to compressed sparse column format */
		status = umfpack_zi_triplet_to_col(NXNY,NXNY,NONZERO,(*waveAC).irow,(*waveAC).icol,(*waveAC).Ar,(*waveAC).Ai,Ap,Ai,Ax,Az,NULL);

		/* Here is something buggy (*waveAC).Ar != Ax and (*waveAC).Ai != Az.
		   Therefore, set  Ax = (*waveAC).Ar and Az = (*waveAC).Ai */
		for (i=0;i<NONZERO;i++){
		     Ax[i] = (*waveAC).Ar[i];
		     Az[i] = (*waveAC).Ai[i];
		}

		/* symbolic factorization */
		status = umfpack_zi_symbolic(NXNY, NXNY, Ap, Ai, Ax, Az, &Symbolic, null, null);

		/* sparse LU decomposition */
		status = umfpack_zi_numeric(Ap, Ai, Ax, Az, Symbolic, &Numeric, null, null);
		umfpack_zi_free_symbolic (&Symbolic);

		/* loop over shots */ 
		for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

		     /* solve forward problem by FDFD */
                     /* ----------------------------- */

		     /* read receiver positions from receiver files for each shot */
		     if(READ_REC==1){

			acq.recpos=receiver(FP, &ntr, 1);

		        (*fwiAC).presr = vector(1,ntr);
		        (*fwiAC).presi = vector(1,ntr);			                         

	      	     }	

		     /* define source vector RHS */
                     RHS_source_AC(waveAC,srcpos,ishot);

		     /* solve forward problem by forward and back substitution */
	    	     status = umfpack_zi_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xi, (*waveAC).RHSr, (*waveAC).RHSi, Numeric, null, null);

		     /* convert vector xr/xi to pr/pi */
		     vec2mat((*waveAC).pr,(*waveAC).pi,xr,xi);

		     /*sprintf(filename,"%s_for_shot_%d.bin",SNAP_FILE,ishot);
		     writemod(filename,(*waveAC).pr,3);*/

		     /* calculate seismograms at receiver positions */		     
		     calc_seis_AC(waveAC,acq.recpos,ntr,ishot,nshots,nfreq);

		     /* store forward wavefield */
		     store_mat((*waveAC).pr, (*fwiAC).forwardr, NX, NY);
		     store_mat((*waveAC).pi, (*fwiAC).forwardi, NX, NY);

		     /* calculate FD residuals at receiver positions */
		     /* -------------------------------------------- */
		     L2 += calc_res_AC(fwiAC,waveAC,ntr,ishot,nstage,nfreq);		     

		     /* de-allocate memory */
		     if(READ_REC==1){
		        free_imatrix(acq.recpos,1,3,1,ntr);
		        free_vector((*fwiAC).presr,1,ntr);
		        free_vector((*fwiAC).presi,1,ntr);
		        ntr=0;
		     }

		} /* end of loop over shots (forward) */

	umfpack_zi_free_numeric (&Numeric);
 
	} /* end of loop over frequencies */

        /* assemble objective function from all MPI processes */

	/* printf("L2 before MPI_Allreduce = %e  on MYID = %d \n", L2, MYID); */

	L2sum = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        L2 = 0.5 * L2sum;	

	/* printf("L2 after MPI_Allreduce = %e  on MYID = %d \n", L2, MYID); */

	/* free memory */
    	free(Ap); free(Ai); free(Ax); free(Az); free(xr); free(xi); 

return L2;
                	    
}
