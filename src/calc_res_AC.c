/*------------------------------------------------------------------------
 *  Calculate FD residuals and objective function                           
 *  
 *  D. Koehn
 *  Kiel, 23/06/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float calc_res_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, int ntr, int ishot, int nstage, int nfreq){

        /* global variables */
	extern int LAPLACE, STF_INV, INVMAT;
	extern char DATA_DIR[STRING_SIZE];

	/* local variables */
	int i;
        char pickfile_char[STRING_SIZE];
        float l2, pobsr, pobsi;
        complex float res, wiennom, wiendenom, wien;

        FILE *fp;

        /* open FA pick file */
      	sprintf(pickfile_char,"%s_p_shot_%d_stage_%d_nfreq_%d.bin",DATA_DIR,ishot,nstage,nfreq);
	fp=fopen(pickfile_char,"rb");
	if (fp == NULL) {
		err(" field data file could not be opened !");
	}

	l2 = 0.0;

        /* initiate wiener deconvolution */
        wiennom = 0.0 + 0.0 * I;
        wiendenom = 0.0 + 0.0 * I;

	for(i=1;i<=ntr;i++){

            /* read FD seismic data */
            //fscanf(fp,"%e%e",&pobsr,&pobsi);
	    fread(&pobsr, sizeof(float), 1, fp);  
	    fread(&pobsi, sizeof(float), 1, fp);  

	    /*printf("resr = %e \n",creal(res));
	    printf("resi = %e \n",cimag(res));*/
		
	    // estimate nominator and denominator of Wiener deconvolution
	    if(STF_INV==1){
	        if((pobsr>1e-20)&&(pobsi>1e-20)){
		    wiennom += conj((*waveAC).precr[i] + (*waveAC).preci[i] * I) * (pobsr + pobsi * I);
		    wiendenom += conj(pobsr + pobsi * I) * (pobsr + pobsi * I);
	        } 
 	    }

	}

        fclose(fp);

        /* estimate STF */
        if(STF_INV==1){
            wien = wiennom/wiendenom;	
	}

	if(STF_INV==0){
	    wien = 1.0 + 0.0 * I;
        }

        (*fwiAC).stfr = creal(wien);
        (*fwiAC).stfi = cimag(wien);

	// printf("stfr = %e \t stfi = %e \n",(*fwiAC).stfr,(*fwiAC).stfi);

	sprintf(pickfile_char,"%s_p_shot_%d_stage_%d_nfreq_%d.bin",DATA_DIR,ishot,nstage,nfreq);
	fp=fopen(pickfile_char,"rb");

	for(i=1;i<=ntr;i++){

            /* read FD seismic data */
            //fscanf(fp,"%e%e",&pobsr,&pobsi);
	    fread(&pobsr, sizeof(float), 1, fp);  
	    fread(&pobsi, sizeof(float), 1, fp);  

	    /*printf("resr = %e \n",creal(res));
	    printf("resi = %e \n",cimag(res));*/

	    if(LAPLACE==0){

		/* calculate complex data residuals ... */
    		/* ... for FDFD-FWI */
		if(INVMAT==1){res = ((pobsr + pobsi * I) - (wien*((*waveAC).precr[i] + (*waveAC).preci[i] * I)))/cabsf(wien);}

    		/* ... for FDFD-RTM */
		if(INVMAT==2){res = pobsr + pobsi * I;}

	    	(*fwiAC).presr[i] = creal(res);
	    	(*fwiAC).presi[i] = cimag(res);

		if((pobsr<1e-20)&&(pobsi<1e-20)){
		   (*fwiAC).presr[i] = 0.0;
	    	   (*fwiAC).presi[i] = 0.0;
		}

	    }

	    if(LAPLACE==1){
	    	// (*fwiAC).presr[i] = log((*waveAC).precr[i]/pobsr)/(*waveAC).precr[i];
	    	(*fwiAC).presr[i] = pobsr - (*waveAC).precr[i];
	    	(*fwiAC).presi[i] = 0.0;
	    }


            /* calculate objective function */
	    if(LAPLACE==0){l2 += creal(res * conj(res));}
	    // if(LAPLACE==1){l2 += pow(log((*waveAC).precr[i]/pobsr),2.0);}
	    if(LAPLACE==1){l2 += pow(pobsr - (*waveAC).precr[i],2.0);}
 
	}

        fclose(fp);
	/*printf("l2 = %e \n",l2);*/

        return l2;
}
