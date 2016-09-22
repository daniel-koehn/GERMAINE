/*--------------------------------------------------------------------------------------------------------------
 *  Calculate diagonal elements of approximate Hessian according to Pratt et al. (1998), Operto et al. (2006)
 *
 *  D. Koehn
 *  Kiel, 18.09.2016
 *  ------------------------------------------------------------------------------------------------------------*/

#include "fd.h"

void hessian_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage){

	/* declaration of global variables */
        extern int SEISMO, NX, NY, NSHOT1, NSHOT2, NF, INFO, NONZERO, NXNY;
        extern int SEISMO, MYID, INFO, NF, LAPLACE, N_STREAMER, READ_REC;
	extern int SPATFILTER, SWS_TAPER_GRAD_HOR, SWS_TAPER_FILE;
        extern float DH, FC_low, FC_high;
	extern char SNAP_FILE[STRING_SIZE];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
    
        /* declaration of local variables */
        int ishot, i, j, nfreq, tr; 
        int status, nxsrc, nysrc;
    	double *null = (double *) NULL ;
	int     *Ap, *Ai;
	double  *Ax, *Az, *xr, *xi; 
	double  time1, time2;
	complex float J;
	float tmp;
        char filename[STRING_SIZE];
        void *Symbolic, *Numeric;
	FILE *FP1;

	/* Allocate memory for compressed sparse column form and solution vector */
	Ap = malloc(sizeof(int)*(NXNY+1));
	Ai = malloc(sizeof(int)*NONZERO);
	Ax = malloc(sizeof(double)*NONZERO);
	Az = malloc(sizeof(double)*NONZERO);
	xr = malloc(sizeof(double)*NONZERO);
	xi = malloc(sizeof(double)*NONZERO);	

        /* suppress modelled seismogram output during Hessian construction */
        SEISMO = 0;

	/* set Hessian to zero */
	init_grad((*fwiAC).hess);

	/* estimate frequency sample interval and set first frequency */
	(*waveAC).dfreq = (FC_high-FC_low) / NF;
        (*waveAC).freq = FC_low;

	/* loop over frequencies at each stage */
	for(nfreq=1;nfreq<=NF;nfreq++){

		/* set squared angular frequency*/
		(*waveAC).omega2 = pow(2.0*M_PI*(*waveAC).freq,2.0);

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
		   printf("\n========================================== \n");
	  	   printf("\n   *** Evaluate approximate Hessian ***    \n");
		   printf("\n========================================== \n");
	        }

		/* loop over shots */ 
		for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

		     /* solve forward problem by FDFD */
                     /* ----------------------------- */

		     /* read receiver positions from receiver files for each shot */
		     if(READ_REC==1){

			acq.recpos=receiver(FP, &ntr, 1);

		        /*(*fwiAC).presr = vector(1,ntr);
		        (*fwiAC).presi = vector(1,ntr);*/
	
			/* Allocate memory for FD seismograms */
			/*alloc_seis_AC(waveAC,ntr);*/

			if(N_STREAMER>0){
			  free_imatrix(acq.recpos,1,3,1,ntr);
			}			                         

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
		     calc_seis_AC(waveAC,acq.recpos,ntr);

		     /* store forward wavefield */
		     store_mat((*waveAC).pr, (*fwiAC).forwardr, NX, NY);
		     store_mat((*waveAC).pi, (*fwiAC).forwardi, NX, NY);

		     /* calculate Greens function at receiver positions */
		     /* ----------------------------------------------- */
			
		     for (tr=1;tr<=ntr;tr++){
		
		         /* define source vector RHS for adjoint wavefield */
			 RHS_source_AC_hess(waveAC,fwiAC,acq.recpos,tr);

 		         /* solve adjoint problem by forward and back substitution */
	    	         status = umfpack_zi_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xi, (*waveAC).RHSr, (*waveAC).RHSi, Numeric, null, null);

		         /* convert vector xr/xi to pr/pi */
		         vec2mat((*waveAC).pr,(*waveAC).pi,xr,xi);		        

			 /* assemble approximate Hessian */
			 for (i=1;i<=NX;i++){
			     for (j=1;j<=NY;j++){

			       J = - 2.0 * (*waveAC).omega2 * ((*fwiAC).forwardr[j][i] + (*fwiAC).forwardi[j][i] * I)  
				         * (1.0 / ((*matAC).vp[j][i] * (*matAC).vp[j][i] * (*matAC).vp[j][i])) * ((*waveAC).pr[j][i] + (*waveAC).pi[j][i] * I);

			       (*fwiAC).hess[j][i] += creal( J * conjf(J));

			     }
			 }

		         /* de-allocate memory */
		         if(READ_REC==1){
		             free_imatrix(acq.recpos,1,3,1,ntr);
			     free_vector((*waveAC).precr,1,ntr);
		             free_vector((*waveAC).preci,1,ntr);
		             free_vector((*fwiAC).presr,1,ntr);
		             free_vector((*fwiAC).presi,1,ntr);
		             ntr=0;
		         }

		     }

		} /* end of loop over shots (forward and adjoint) */  

	(*waveAC).freq += (*waveAC).dfreq; 

	umfpack_zi_free_numeric (&Numeric);
 
	} /* end of loop over frequencies */

	/* Assemble Hessian from all MPI processes  */
	sum_grad_MPI((*fwiAC).hess);
	//cp_grad_frame((*fwiAC).hess);

	/* output of Hessian */
	sprintf(filename,"%s_hess",JACOBIAN);
	FP1=fopen(filename,"wb");
	
	for (i=1;i<=NX;i++){   
           for (j=1;j<=NY;j++){
		 tmp = (*fwiAC).hess[j][i];
                 fwrite(&tmp,sizeof(float),1,FP1);
	   }
        }
        
	fclose(FP1);

	/* free memory */
    	free(Ap); free(Ai); free(Ax); free(Az); free(xr); free(xi);
                	    
}
