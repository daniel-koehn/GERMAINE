/*------------------------------------------------------------------------
 *   Read time-domain source wavelet from SU file 
 *   and apply discrete Fourier transform for given 
 *   frequency 
 *    
 *   Daniel Koehn                              
 *   Kiel, 15.07.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "segy.h"

void read_stf_dft(struct waveAC *waveAC, float *amp){

	/* declaration of global variables */
        extern char WAVELET_NAME[STRING_SIZE];
	extern int PHYSICS;

	/* declaration of local variables */
	int i, j, h, ns, ndt;
	segy tr;
	float dump, *wavelet, dt, t;
        FILE *fpdata;		

	/* Read source wavelet from SU file */
	/* -------------------------------- */
	
	fpdata = fopen(WAVELET_NAME,"r");
        if (fpdata==NULL) err(" Time domain source wavelet SU file could not be opened !");

	/* read header */   			
	fread(&tr,240,1,fpdata);

	/* read number of time samples ns and sample interval dt */
	ns = tr.ns;

	ndt = tr.dt;

	if(PHYSICS<=3){dt = ndt / 1e6;  /* convert units [mus] to [s] */}
	if(PHYSICS==4){dt = ndt / 1e9;  /* convert units [ns] to [s] */}

	/* read source wavelet */
	fread(&tr.data,4,ns,fpdata);

	/* allocate memory for source wavelet */
	wavelet = vector(1,ns);						
			  
	/* read source wavelet */
	h=ns;
	for(j=1;j<=ns;j++){
	    dump=tr.data[j];
	    wavelet[h]=dump;
	    h--;
	}

	fclose(fpdata);

	/* Calculate DFT of time-domain source wavelet */
	/* ------------------------------------------- */

	amp[1] = 0.0;
	amp[2] = 0.0;

	for(j=1;j<=ns;j++){
   
                t = j * dt;
            
                amp[1] += wavelet[j] * cos(2.0 * t * (*waveAC).freq * M_PI) * dt;
                amp[2] += wavelet[j] * sin(2.0 * t * (*waveAC).freq * M_PI) * dt;
 
	}

	/* free memory for source wavelet */
	free_vector(wavelet,1,ns);  
			
}
