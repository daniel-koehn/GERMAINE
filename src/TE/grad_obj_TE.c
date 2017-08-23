/*------------------------------------------------------------------------
 *  Calculate gradient and objective function value for TE-mode problem
 *
 *  D. Koehn
 *  Kiel, 20.08.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float grad_obj_TE(struct fwiTE *fwiTE, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matTE *matTE, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage){

	/* declaration of global variables */
        extern int SEISMO, NX, NY, NSHOT1, NSHOT2, NF, INFO, NONZERO, NXNY;
        extern int SEISMO, MYID, INFO, NF, N_STREAMER, READ_REC;
	extern int SPATFILTER, SWS_TAPER_GRAD_HOR, SWS_TAPER_FILE;
	extern int NFREQ1, NFREQ2;
        extern float DH, FC_low, FC_high, MAT1_NORM, MAT2_NORM;
	extern float SCALE_GRAD1, SCALE_GRAD2;
	extern char SNAP_FILE[STRING_SIZE];
	extern FILE *FP;
    
        /* declaration of local variables */
        int ishot, i, nfreq; 
        float L2sum, L2, ** grad_shot_sigma, ** grad_shot_epsilon;
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

        grad_shot_sigma =  matrix(1,NY,1,NX);
        grad_shot_epsilon =  matrix(1,NY,1,NX);

        /* suppress modelled seismogram output during FWI */
        SEISMO = 0;

	/* set gradient matrix and objective function to zero before next iteration */
	init_grad((*fwiTE).grad_sigma);
	init_grad((*fwiTE).grad_epsilon);
        L2 = 0.0;

	/* loop over frequencies at each stage */
	for(nfreq=NFREQ1;nfreq<NFREQ2;nfreq++){

		/* set frequency on local MPI process */
		(*waveAC).freq = (*waveAC).stage_freq[nfreq];

		/* define PML damping profiles */
		pml_pro(PML_AC,waveAC);

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

		if(MYID==0){
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

		if(MYID==0){
		  time2=MPI_Wtime();
		  printf("\n Finished after %4.2f s \n",time2-time1);
		}

	        if(MYID==0){
		   printf("\n============================================================================================================== \n");
	  	   printf("\n   *** Evaluate gradient and objective function for shot no. %d - %d on MPI-process %d (Freq. %d of %d) ***    \n", NSHOT1, NSHOT2-1, MYID, nfreq, NF);
		   printf("\n============================================================================================================== \n");
	        }

		/* loop over shots */ 
		for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

		     /* set gradient for each shot to zero before next iteration */
		     init_grad(grad_shot_sigma);
  		     init_grad(grad_shot_epsilon);

		     /* solve forward problem by FDFD */
                     /* ----------------------------- */

		     /* read receiver positions from receiver files for each shot */
		     if(READ_REC==1){

			acq.recpos=receiver(FP, &ntr, 1);

		        (*waveAC).presr = vector(1,ntr);
		        (*waveAC).presi = vector(1,ntr);				                         

	      	     }

		     /* define source vector RHS */
                     RHS_source_AC(waveAC,srcpos,ishot);

		     /* solve forward problem by forward and back substitution */
	    	     status = umfpack_zi_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xi, (*waveAC).RHSr, (*waveAC).RHSi, Numeric, null, null);

		     /* convert vector xr/xi to pr/pi */
		     vec2mat((*waveAC).pr,(*waveAC).pi,xr,xi);

		     /*sprintf(filename,"%s_for_shot_%d.bin",SNAP_FILE,ishot);
		     writemod(filename,(*waveTE).pr,3);*/

		     /* calculate seismograms at receiver positions */
		     calc_seis_AC(waveAC,acq.recpos,ntr,ishot,nshots,nfreq);

		     /* store forward wavefield */
		     store_mat((*waveAC).pr, (*fwiTE).forwardr, NX, NY);
		     store_mat((*waveAC).pi, (*fwiTE).forwardi, NX, NY);

		     /* calculate FD residuals at receiver positions */
		     /* -------------------------------------------- */
		     L2 += calc_res_AC(waveAC,ntr,ishot,nstage,nfreq);

		     /* solve adjoint problem by FDFD */
		     /* ----------------------------- */

		     /* define source vector RHS for adjoint wavefield */
		     RHS_source_AC_adj(waveAC,acq.recpos,ntr);

 		     /* solve adjoint problem by forward and back substitution */
	    	     status = umfpack_zi_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xi, (*waveAC).RHSr, (*waveAC).RHSi, Numeric, null, null);

		     /* convert vector xr/xi to pr/pi */
		     vec2mat((*waveAC).pr,(*waveAC).pi,xr,xi);

		     /*sprintf(filename,"%s_for_shot_%d.bin",SNAP_FILE,ishot);
		     writemod(filename,(*waveAC).pr,3);*/

		     /* assemble gradient for each shot */
		     ass_grad_TE(fwiTE, waveAC, matTE, grad_shot_sigma, grad_shot_epsilon, srcpos, nshots, acq.recpos, ntr, ishot);

		     /* de-allocate memory */
		     if(READ_REC==1){
		        free_imatrix(acq.recpos,1,3,1,ntr);
		        free_vector((*waveAC).presr,1,ntr);
		        free_vector((*waveAC).presi,1,ntr);
		        ntr=0;
		     }

		} /* end of loop over shots (forward and adjoint) */

	umfpack_zi_free_numeric (&Numeric);
 
	} /* end of loop over frequencies */

        /* assemble objective function from all MPI processes */
	/* printf("L2 before MPI_Allreduce = %e  on MYID = %d \n", L2, MYID); */

	/* Assemble gradients from all MPI processes  */
	sum_grad_MPI((*fwiTE).grad_sigma);
	cp_grad_frame((*fwiTE).grad_sigma);

	sum_grad_MPI((*fwiTE).grad_epsilon);
	cp_grad_frame((*fwiTE).grad_epsilon);

	L2sum = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        L2 = 0.5 * L2sum;

	/* calculate gradients for normalized material parameters */	
	scale_grad((*fwiTE).grad_sigma,MAT1_NORM,(*fwiTE).grad_sigma,NX,NY);
	scale_grad((*fwiTE).grad_epsilon,MAT2_NORM,(*fwiTE).grad_epsilon,NX,NY);

	/* calculate L2norm of gradients */
	SCALE_GRAD1 = norm_matrix((*fwiTE).grad_sigma,NX,NY); 
	SCALE_GRAD2 = norm_matrix((*fwiTE).grad_epsilon,NX,NY);

	/* printf("L2 after MPI_Allreduce = %e  on MYID = %d \n", L2, MYID); */

        /* deallocate memory */
        free_matrix(grad_shot_sigma,1,NY,1,NX);
        free_matrix(grad_shot_epsilon,1,NY,1,NX);

	/* free memory */
    	free(Ap); free(Ai); free(Ax); free(Az); free(xr); free(xi);

return L2;
                	    
}
