/*---------------------------------------------------------------------------------------
 *  GERMAINE: 2D frequency domain finite-difference modelling an FWI for acoustic media 
 *
 *  Authors:
 *  -----------  
 * 
 *  D. Koehn    (FDFD code + updates)
 *  D. De Nil   (FDFD code + updates)
 *  
 *  In case of questions contact the author:
 *	Dr. Daniel Koehn, Kiel University, Institute of Geoscience,
 *	Otto-Hahn-Platz 1, D-24098 Kiel, Germany, ph: +49 431 880 3878,
 *	mailto:dkoehn@geophysik.uni-kiel.de,
 *	Homepage: http://www.geophysik.uni-kiel.de/~dkoehn
 *
 *  GERMAINE is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, version 2.0 of the License only. 
 *  
 *  GERMAINE is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *  
 *  You should have received a copy of the GNU General Public License 
 *  along with GERMAINE (see file LICENSE.md) 
 *
 *  If you show modelling/inversion results in a paper or presentation please 
 *  give a reference to the following papers: 
 *  
 *  Thank you for your co-operation, 
 *  Daniel Koehn
 * 
 *  ---------------------------------------------------------------------------------------*/

#include "fd.h"           /* general include file */
#include "globvar.h"      /* definition of global variables  */

int main(int argc, char **argv){

/* variables in main */
double time1, time8;
char ext[10], *fileinp, *fileinp1;	

/* Initialize MPI environment */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

setvbuf(stdout, NULL, _IONBF, 0);

if (MYID == 0){
   time1=MPI_Wtime(); 
   clock();
}
		
/* print program name, version etc to stdout*/
if (MYID == 0) info(stdout);

/* read parameters from parameter-file (stdin) */
fileinp=argv[1];
fileinp1=argv[2];
FP=fopen(fileinp,"r");
if(FP==NULL) {
	if (MYID == 0){
		printf("\n==================================================================\n");
		printf(" Cannot open GERMAINE input file %s \n",fileinp);
		printf("\n==================================================================\n\n");
		err(" --- ");
	}
}


/* read input file *.inp */
read_par(FP);

if (MYID == 0) note(stdout);

/* define frequency domain or Laplace domain modelling/FWI */
LAPLACE = 0;

/* solve acoustic forward problem by FDFD */
/* -------------------------------------- */
if(INVMAT==0){
  forward_AC(fileinp1);
}

/* 2D acoustic FDFD Full Waveform Inversion */
/* ---------------------------------------- */
if(INVMAT==1){
   fwi_FD_AC(fileinp1);
}
          
/* 1D linear gradient model estimation by grid search */
/* -------------------------------------------------- */
/*if(INVMAT==2){
   grid_search(Vp,S,TT,Tmod,Tobs,Tres,srcpos,nshots,recpos,ntr);
}*/

MPI_Barrier(MPI_COMM_WORLD);

if(MYID==0){
	printf("\n **Info from main (written by PE %d): \n",MYID);
	printf(" CPU time of program per PE: %li seconds.\n",clock()/CLOCKS_PER_SEC);
	time8=MPI_Wtime();
	printf(" Total real time of program: %4.2f seconds.\n",time8-time1);		
}

// fclose(FP1);

MPI_Finalize();
return 0;	

}/*main*/
