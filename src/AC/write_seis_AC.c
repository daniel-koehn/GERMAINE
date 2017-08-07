/*------------------------------------------------------------------------
 *  Write FD seismograms to file
 *  
 *  D. Koehn
 *  Kiel, 21/06/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void write_seis_AC(struct waveAC *waveAC, int nshots, int ntr, int nstage){

	/* global variables */
        extern char PICKS_FILE[STRING_SIZE];        
        extern int MYID, INFO, NONZERO, NF;

	/* local variables */
	int i, j, k, index;
        char pickfile_char[STRING_SIZE];
	float tmp, tmp1;
	float * precr_vec_local, * precr_vec, * preci_vec_local, * preci_vec;		

        if((MYID==0)&&(INFO==1)){
	    printf("\n Write FD seismograms for stage %d ...\n",nstage);
        }

	/* gather data from all MPI processes */
	precr_vec = vector(1,ntr*NF*nshots);
	precr_vec_local = vector(1,ntr*NF*nshots);
	preci_vec = vector(1,ntr*NF*nshots);
	preci_vec_local = vector(1,ntr*NF*nshots);


	for(i=1;i<=(ntr*nshots*NF);i++){

	    precr_vec_local[i] = (*waveAC).precr[i];
	    preci_vec_local[i] = (*waveAC).preci[i];

	    precr_vec[i] = 0.0;
	    preci_vec[i] = 0.0;

	}

	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Allreduce(precr_vec_local,precr_vec,ntr*nshots*NF,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(preci_vec_local,preci_vec,ntr*nshots*NF,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

        if(MYID==0){	

           FILE *fp;
           sprintf(pickfile_char,"%s_p_stage_%d.bin",PICKS_FILE,nstage);
           fp=fopen(pickfile_char,"wb");

	   for(k=1;k<=NF;k++){ 	
	       for(j=1;j<=nshots;j++){
	           for(i=1;i<=ntr;i++){ 

		       index =  i + ntr * (j-1) + ntr * nshots * (k-1);

                       tmp = precr_vec[index]; 
	               tmp1 = preci_vec[index];
 
	               fwrite(&tmp, sizeof(float), 1, fp);
	               fwrite(&tmp1, sizeof(float), 1, fp);

		   }
	       }
	   }

	   fclose(fp);
	}        

	free_vector(precr_vec,1,ntr*NF*nshots);
	free_vector(precr_vec_local,1,ntr*NF*nshots);	

	free_vector(preci_vec,1,ntr*NF*nshots);
	free_vector(preci_vec_local,1,ntr*NF*nshots);	

}
