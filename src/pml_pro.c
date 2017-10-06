/*----------------------------------------------------------------------------------------------------------------
 *  Define damping profiles within PML boundary condition according to 
 * 
 *  Hustedt, B., Operto, S. and Virieux, J. (2004) Mixed-grid and staggered-grid finite difference
 *  methods for frequency domain acoustic wave modelling. Geophysical Journal International, 157:1269–1296.
 *
 *  Operto, S. and Virieux, J. (2006) Practical aspects of frequency-domain finite-difference modelling of 
 *  seismic wave propagation. Lecture Notes summer school "Seismic Imaging of Complex Structures from Multicomponent 
 *  global offset data by Full Waveform Inversion"
 *  
 *  Operto, S., Virieux, J., Ribodetti, A. and Anderson, J.A. (2009) Finite-difference frequency-domain modeling 
 *  of viscoacoustic wave propagation in 2D tilted transversely isotropic (TTI) media . Geophysics 74(5):T75-T95.
 *  
 *  D. Koehn
 *  Kiel, 15.08.2017
 *  -------------------------------------------------------------------------------------------------------------- */

#include "fd.h"

void pml_pro(struct PML_AC *PML_AC, struct waveAC *waveAC){

	extern int NX, NY, NPML, FREE_SURF;
        extern float DH, A0_PML, S;
	extern char SNAP_FILE[STRING_SIZE];	

	/* local variables */
	int i;
        complex float ci, OMEGA_PML, ax, ay;
        float lPML, x, y, xh, yh, pih, betad, betadh, betam;
	float xmax, xh1;
	float alpha, alphah, eps, epsh;
	char filename[STRING_SIZE];
	FILE *fpmod;

	pih = M_PI / 2.0;

        /* define complex i */
        ci = 0.0 + 1.0 * I;

	/* define thickness of PML [m] */
	lPML = NPML * DH;

	/* define PML damping parameters */
	betam = 0.0;
        OMEGA_PML = 2.0 * M_PI * ((*waveAC).freq - I * S);	

	/* initialize PML profiles */
	for (i=0;i<=NX+1;i++){;
	    (*PML_AC).dampxr[i] = 1.0;
	    (*PML_AC).dampxi[i] = 0.0;
	    (*PML_AC).dampxhr[i] = 1.0;
	    (*PML_AC).dampxhi[i] = 0.0;
	}

	for (i=0;i<=NY+1;i++){;
	    (*PML_AC).dampyr[i] = 1.0;
	    (*PML_AC).dampyi[i] = 0.0;
	    (*PML_AC).dampyhr[i] = 1.0;
	    (*PML_AC).dampyhi[i] = 0.0;
	}

        /* calculate PML damping profiles in x- and y-direction */
	/* ---------------------------------------------------- */
	for (i=1;i<=NPML;i++){  

	    x = i * DH;
	    xh = x + 0.5 * DH;
     	    
	    /* define PML damping profiles in x- and y-direction */

	    /* "classical" PML profiles (Berenger 1994, Chew & Weedon 1994) */
	    /* ------------------------------------------------------------ */

	    /* compute eps and alpha */
            eps = A0_PML * (1.0 - cos((lPML-x) * pih / lPML));
	    	
	    /* compute epsh and alphah */
            epsh = A0_PML * (1.0 - cos((lPML-xh) * pih / lPML));
	    
            /*(*PML_AC).dampxr[i] = creal( 1.0 / (1.0  + ci * eps  / OMEGA_PML));            
            (*PML_AC).dampxi[i] = cimag( 1.0 / (1.0  + ci * eps  / OMEGA_PML));            

            (*PML_AC).dampxhr[i] = creal( 1.0 / (1.0 + ci * epsh / OMEGA_PML));
            (*PML_AC).dampxhi[i] = cimag( 1.0 / (1.0 + ci * epsh / OMEGA_PML));*/

            (*PML_AC).dampxr[i] = creal( 1.0 / (1.0  + eps  / (ci * OMEGA_PML)));            
            (*PML_AC).dampxi[i] = cimag( 1.0 / (1.0  + eps  / (ci * OMEGA_PML)));            

            (*PML_AC).dampxhr[i] = creal( 1.0 / (1.0 + epsh / (ci * OMEGA_PML)));
            (*PML_AC).dampxhi[i] = cimag( 1.0 / (1.0 + epsh / (ci * OMEGA_PML)));

	    /* C-PML profiles (Kuzuoglu & Mittra; 1996; Roden & Gedney; 2000; Komatitsch & Martin; 2007) */
	    /* ----------------------------------------------------------------------------------------- */

	    /*betad = x * betam / lPML;
	    betadh = xh * betam / lPML;

	    alpha = 1.0;
	    alphah = 1.0;

            (*PML_AC).dampxr[i] = creal( 1.0 / (alpha  + ci * eps  / (OMEGA_PML + ci * betad)));            
            (*PML_AC).dampxi[i] = cimag( 1.0 / (alpha  + ci * eps  / (OMEGA_PML + ci * betad)));            

            (*PML_AC).dampxhr[i] = creal( 1.0 / (alphah + ci * epsh / (OMEGA_PML + ci * betadh)));
            (*PML_AC).dampxhi[i] = cimag( 1.0 / (alphah + ci * epsh / (OMEGA_PML + ci * betadh)));*/

	    (*PML_AC).dampyr[i] = (*PML_AC).dampxr[i];
	    (*PML_AC).dampyi[i] = (*PML_AC).dampxi[i];
	    (*PML_AC).dampyhr[i] = (*PML_AC).dampxhr[i];
	    (*PML_AC).dampyhi[i] = (*PML_AC).dampxhi[i];


            /* mirror PMLs at opposite site of computational domain */
	    (*PML_AC).dampxr[NX-i+1] = (*PML_AC).dampxr[i];
	    (*PML_AC).dampxi[NX-i+1] = (*PML_AC).dampxi[i];
	    (*PML_AC).dampxhr[NX-i+1] = (*PML_AC).dampxhr[i];
	    (*PML_AC).dampxhi[NX-i+1] = (*PML_AC).dampxhi[i];

	    (*PML_AC).dampyr[NY-i+1] = (*PML_AC).dampxr[i];
	    (*PML_AC).dampyi[NY-i+1] = (*PML_AC).dampxi[i];
	    (*PML_AC).dampyhr[NY-i+1] = (*PML_AC).dampxhr[i];
	    (*PML_AC).dampyhi[NY-i+1] = (*PML_AC).dampxhi[i];

	    /* deactivate PML for y <= NPML if FREE_SURF == 1*/
	    if(FREE_SURF==1){
	        (*PML_AC).dampyr[i] = 1.0;
	        (*PML_AC).dampyi[i] = 0.0;
	        (*PML_AC).dampyhr[i] = 1.0;
	        (*PML_AC).dampyhi[i] = 0.0;
	    }

	}

	(*PML_AC).dampxr[0] = (*PML_AC).dampxr[1];
	(*PML_AC).dampxi[0] = (*PML_AC).dampxi[1];
	(*PML_AC).dampxhr[0] = (*PML_AC).dampxhr[1];
	(*PML_AC).dampxhi[0] = (*PML_AC).dampxhi[1];

	(*PML_AC).dampxr[NX+1] = (*PML_AC).dampxr[NX];
	(*PML_AC).dampxi[NX+1] = (*PML_AC).dampxi[NX];
	(*PML_AC).dampxhr[NX+1] = (*PML_AC).dampxhr[NX];
	(*PML_AC).dampxhi[NX+1] = (*PML_AC).dampxhi[NX];

	(*PML_AC).dampyr[0] = (*PML_AC).dampyr[1];
	(*PML_AC).dampyi[0] = (*PML_AC).dampyi[1];
	(*PML_AC).dampyhr[0] = (*PML_AC).dampyhr[1];
	(*PML_AC).dampyhi[0] = (*PML_AC).dampyhi[1];

	(*PML_AC).dampyr[NY+1] = (*PML_AC).dampyr[NY];
	(*PML_AC).dampyi[NY+1] = (*PML_AC).dampyi[NY];
	(*PML_AC).dampyhr[NY+1] = (*PML_AC).dampyhr[NY];
	(*PML_AC).dampyhi[NY+1] = (*PML_AC).dampyhi[NY];        


	/*fpmod=fopen("pml.dat","w");
	for (i=0;i<=NY+1;i++){fprintf(fpmod,"%e \n",(*PML_AC).dampyhr[i]);}						
	fclose(fpmod);*/

       
}




