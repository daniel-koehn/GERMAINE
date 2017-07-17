/*------------------------------------------------------------------------
 *  Initiate impedance matrix for 2D elastic SH problem with PMLs 
 *  using a 9p-mixed grid according to 
 * 
 *  Z. Chen, D. Cheng, W. Feng, H. Yang, 2013, An optimal 9-point finite 
 *  difference scheme for the Helmholtz equation with PML, Int. J. Numer. 
 *  Anal. Model., 10, 389-410.
 *
 *  D. Koehn
 *  Kiel, 26.05.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void init_A_SH_9p_pml(struct PML_AC *PML_AC, struct matSH *matSH, struct waveAC *waveAC){

	extern int NX, NY, NONZERO;
        extern float DH, S;	

	/* local variables */
	int i, j, k;
        complex float tmp, Omega2;
        complex float tmpA, tmpB, tmpC;
        float b, d, e, idh2;
	//SuiteSparse_long count;
	int count, ishift;

        /* define FD parameters */
	b = 0.7926;
	d = 0.3768;
        e = -0.0064;

        idh2 = 1.0/(DH*DH);

	/* assemble impedance matrix */        

        k=0;        /* index of main diagonal  */	
        count=0;    /* count non-zero elements */

	/* set squared complex angular frequency*/
	Omega2 = cpowf(((2.0*M_PI*(*waveAC).freq) + (I * S)),2.0);

        for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		    /* NW gridpoint */
		    if((i > 1) && (j > 1)){

                       tmpA = (*PML_AC).Axr[j-1][i-1] + (*PML_AC).Axi[j-1][i-1] * I;
		       tmpB = (*PML_AC).Byr[j-1][i-1] + (*PML_AC).Byi[j-1][i-1] * I;
		       tmpC = (*PML_AC).Cr[j-1][i-1] + (*PML_AC).Ci[j-1][i-1] * I;

		       tmp = ((1-b)*idh2/2.0) * (tmpA + tmpB) + (e/4.0) * tmpC * Omega2 * (*matSH).rho[j-1][i-1];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - 1 - NX;
                       count++;

		    }   

		    /* N gridpoint */
		    if(j > 1){

                       tmpA = (*PML_AC).Ar[j-1][i] + (*PML_AC).Ai[j-1][i] * I;
		       tmpB = (*PML_AC).Byr[j-1][i] + (*PML_AC).Byi[j-1][i] * I;
		       tmpC = (*PML_AC).Cr[j-1][i] + (*PML_AC).Ci[j-1][i] * I;

		       tmp = -((1-b)*idh2) * tmpA + (b*idh2) * tmpB + (d/4.0) * tmpC *  Omega2 * (*matSH).rho[j-1][i];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - NX;
                       count++;

		    }  

		    /* NE gridpoint */
		    if((i < NX) && (j > 1)){

                       tmpA = (*PML_AC).Axr[j-1][i] + (*PML_AC).Axi[j-1][i] * I;
		       tmpB = (*PML_AC).Byr[j-1][i+1] + (*PML_AC).Byi[j-1][i+1] * I;
		       tmpC = (*PML_AC).Cr[j-1][i+1] + (*PML_AC).Ci[j-1][i+1] * I;

		       tmp = ((1-b)*idh2/2.0) * (tmpA + tmpB) + (e/4.0) * tmpC * Omega2 * (*matSH).rho[j-1][i+1];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1 - NX;
                       count++;

		    }  

		    /* W gridpoint */
		    if(i > 1){

                       tmpA = (*PML_AC).Axr[j][i-1] + (*PML_AC).Axi[j][i-1] * I;
		       tmpB = (*PML_AC).Br[j][i-1] + (*PML_AC).Bi[j][i-1] * I;
		       tmpC = (*PML_AC).Cr[j][i-1] + (*PML_AC).Ci[j][i-1] * I;

		       tmp = (b*idh2) * tmpA - ((1-b)*idh2) * tmpB + (d/4.0) * tmpC * Omega2 * (*matSH).rho[j][i-1];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k-1;
                       count++;

		    }

		    /* central gridpoint */
                    tmpA = (*PML_AC).Ar[j][i] + (*PML_AC).Ai[j][i] * I;
		    tmpB = (*PML_AC).Br[j][i] + (*PML_AC).Bi[j][i] * I;
		    tmpC = (*PML_AC).Cr[j][i] + (*PML_AC).Ci[j][i] * I;

		    tmp = (1-d-e) * tmpC * Omega2 * (*matSH).rho[j][i] - (2.0*b*idh2) * (tmpA+tmpB);

                    (*waveAC).Ar[count] = creal(tmp); 
                    (*waveAC).Ai[count] = cimag(tmp);
                    (*waveAC).irow[count] = k;
                    (*waveAC).icol[count] = k;
                    count++;

		    /* E gridpoint */
		    if(i < NX){

                       tmpA = (*PML_AC).Axr[j][i] + (*PML_AC).Axi[j][i] * I;
		       tmpB = (*PML_AC).Br[j][i+1] + (*PML_AC).Bi[j][i+1] * I;
		       tmpC = (*PML_AC).Cr[j][i+1] + (*PML_AC).Ci[j][i+1] * I;

		       tmp = (b*idh2) * tmpA - ((1-b)*idh2) * tmpB + (d/4.0) * tmpC * Omega2 * (*matSH).rho[j][i+1];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1;
                       count++;

		    }

		    /* SW gridpoint */
		    if( (j < NY) && (i > 1) ){

                       tmpA = (*PML_AC).Axr[j+1][i-1] + (*PML_AC).Axi[j+1][i-1] * I;
		       tmpB = (*PML_AC).Byr[j][i-1] + (*PML_AC).Byi[j][i-1] * I;
		       tmpC = (*PML_AC).Cr[j+1][i-1] + (*PML_AC).Ci[j+1][i-1] * I;

		       tmp = ((1-b)*idh2/2.0) * (tmpA + tmpB) + (e/4.0) * tmpC * Omega2 * (*matSH).rho[j+1][i-1];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k - 1 + NX;
                       count++;

		    }

		    /* S gridpoint */
		    if(j < NY){

                       tmpA = (*PML_AC).Ar[j+1][i] + (*PML_AC).Ai[j+1][i] * I;
		       tmpB = (*PML_AC).Byr[j][i] + (*PML_AC).Byi[j][i] * I;
		       tmpC = (*PML_AC).Cr[j+1][i] + (*PML_AC).Ci[j+1][i] * I;

		       tmp = -((1-b)*idh2) * tmpA + (b*idh2) * tmpB + (d/4.0) * tmpC * Omega2 * (*matSH).rho[j+1][i];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + NX;
                       count++;

		    }     
            
		    /* SE gridpoint */
		    if((i < NX) && (j < NY)) {

                       tmpA = (*PML_AC).Axr[j+1][i] + (*PML_AC).Axi[j+1][i] * I;
		       tmpB = (*PML_AC).Byr[j][i+1] + (*PML_AC).Byi[j][i+1] * I;
		       tmpC = (*PML_AC).Cr[j+1][i+1] + (*PML_AC).Ci[j+1][i+1] * I;

		       tmp = ((1-b)*idh2/2.0) * (tmpA + tmpB) + (e/4.0) * tmpC * Omega2 * (*matSH).rho[j+1][i+1];

                       (*waveAC).Ar[count] = creal(tmp); 
                       (*waveAC).Ai[count] = cimag(tmp);
                       (*waveAC).irow[count] = k;
                       (*waveAC).icol[count] = k + 1 + NX;
                       count++;

		    }

                    k++;                       

		}
	}

}




