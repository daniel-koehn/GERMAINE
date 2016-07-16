/*
 * taper gradient with a gaussian frame to damp inversion artefacts near the sources and receivers
 sws == 1 vertical taper (for tomography geometry)
 sws == 2 horizontal taper (for reflection geometry)
 sws == 3 local radial taper at the source and receiver positions
 sws == 4  
 */

#include "fd.h"

void taper_grad(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int sws)
{

	/* extern variables */

        extern float DH;
	extern int NX, NY, NPML;
	extern int MYID;
	extern FILE *FP;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2, VTON, SRTON;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1;

	extern int GRADT1, GRADT2, GRADT3, GRADT4;
	float amp, a, *window, grad_tap;
	char modfile[STRING_SIZE];
	
	extern float SRTRADIUS, EXP_TAPER_GRAD_HOR;
	extern int SRTSHAPE, FILTSIZE;
        float x, y, taper;
        float **taper_coeff;
	float R, frh, f0, r;

	FILE *fp_taper;

        taper_coeff= matrix(1,NY,1,NX);


    /* =============== */
    /* Vertical taper  */
    /* =============== */
	
	if(sws==1){
	
	/* define taper geometry */
	taperlength=GRADT2-GRADT1;
	taperlength2=GRADT4-GRADT3;
	ifw = GRADT4-GRADT1;
	
	/*printf("%d \t %d \t %d \t %d \n",GRADT1, GRADT2, GRADT3, GRADT4);
	printf("%d \t %d \t %d \n",taperlength, taperlength2,ifw);*/
	
	if (MYID==0)
	{
		fprintf(FP,"\n **Message from taper_grid (printed by PE %d):\n",MYID);
		fprintf(FP," Coefficients for gradient taper are now calculated.\n");
	}
	
	window=vector(1,ifw);
	 
	 /* Gaussian window */  
        a=3;    /* damping coefficient */
        for (i=1;i<=ifw;i++){
	window[i] = 1.0;
	
	    if(i<=taperlength){
	       window[i] = exp(-(1.0/2.0)*(a*(i-taperlength)/(taperlength/2.0))*(a*(i-taperlength)/(taperlength/2.0)));
	    }
	    
	    if(i>=(ifw-taperlength2)){
	       window[i] = exp(-(1.0/2.0)*(a*(i-(ifw-taperlength2))/(taperlength2/2.0))*(a*(i-(ifw-taperlength2))/(taperlength2/2.0)));
	    }
	    
	}
	
	/* loop over global grid */
	for (j=1;j<=NY;j++){
	
	   h=1;
	   for (i=1;i<=NX;i++){

                        grad_tap=0.0;
			
			if((i>GRADT1)&&(i<GRADT4)){
			   grad_tap=window[h];
			   h++;
			}
			   			

			taper_coeff[j][i]=grad_tap;

			
		}
	}	
	
	/* apply taper on local gradient */
	for (j=1;j<=NY;j++){
	   for (i=1;i<=NX;i++){
           waveconv[j][i]*=taper_coeff[j][i];
           }
	}
	

	sprintf(modfile,"taper_coeff_vert.bin");
	writemod(modfile,taper_coeff,3); 


	free_vector(window,1,ifw);

	}

        
	/* ======================== */
	/* Horizontal Taper         */
	/* ======================== */
	
	if(sws==2){
        /* define taper geometry */
        taperlength=GRADT2-GRADT1;
        taperlength2=GRADT4-GRADT3;
        ifw = GRADT2-GRADT1+1;

        if(MYID==0){
           printf("%d \t %d \t %d \t %d \n",GRADT1, GRADT2, GRADT3, GRADT4);
           printf("%d \t %d \t %d \n",taperlength, taperlength2,ifw);
        }

        if (MYID==0)
        {
                fprintf(FP,"\n **Message from taper_grid (printed by PE %d):\n",MYID);
                fprintf(FP," Coefficients for gradient taper are now calculated.\n");
        }

        window=vector(1,ifw);
         
        /* Gaussian window */  
        a=3;    /* damping coefficient */
        for (i=1;i<=ifw;i++){
        window[i] = 1.0;

            if(i<=taperlength){
               window[i] = exp(-(1.0/2.0)*(a*(i-taperlength)/(taperlength/2.0))*(a*(i-taperlength)/(taperlength/2.0)));
            }
          
        }

        /* loop over global grid */
	h=1;
        for (j=1;j<=NY;j++){
           for (i=1;i<=NX;i++){

                        grad_tap=0.0;

                        if((j>=GRADT1+NPML)&&(j<=GRADT2+NPML)){
                           grad_tap=window[h];
                        }

                        if(j>GRADT2+NPML){
                          
                          /*grad_tap=((float)(j*DH))*((float)(j*DH))*((float)(j*DH));*/
			  /*grad_tap=(float)(j*DH);*/ 
			  
			  /*grad_tap=((float)((j*j)/(GRADT2*GRADT2)));*/
			  /*grad_tap=1.0;*/

                          grad_tap=pow((float)(j*DH),EXP_TAPER_GRAD_HOR);

                        }
                        
                        /*grad_tap=((float)(j*DH));*/
			

                        taper_coeff[j][i]=grad_tap;
                  
                }
		
	   if(j>=GRADT1+NPML){h++;}
        }

        /* apply taper on local gradient */
        for (j=1;j<=NY;j++){
           for (i=1;i<=NX;i++){
           waveconv[j][i]*=taper_coeff[j][i];
           }
        }

        sprintf(modfile,"taper_coeff_hor.bin");
        writemod(modfile,taper_coeff,3);

        free_vector(window,1,ifw);

        } /* end of sws==2 */

        /* =================================== */
        /* taper source and receiver positions */
	/* =================================== */
	
	if(sws==3) {

        R = SRTRADIUS;
	frh = 0.1;
	f0 = 1e-3;
	a = - (log(frh-f0)-log(1.0-f0))/log(2.0);

        /*****************************/
        /* Taper at source positions */
        /*****************************/
	for (n=1;n<=nshots;n++) {
	     for (iy=1;iy<=NY;iy++){
                  for (ix=1;ix<=NX;ix++){
	
	               /* calculate global coordinates */
		       x = ix * DH;
		       y = iy * DH; 
		       r = sqrt(pow((x-srcpos[1][n]),2.0)+pow((y-srcpos[2][n]),2.0));
		       taper = 1.0;
		  
		       if(r<=R){
		          taper = f0 + (1.0 - f0) * pow((0.5 - 0.5 *cos(PI*r/R)),a);		 
		       }
		  
		       taper_coeff[iy][ix] = taper;
		       waveconv[iy][ix] *= taper;
		       
                  }
	     }
        }
	
        /*****************************/
        /* Taper at receiver positions */
        /*****************************/
	for (n=1;n<=ntr;n++) {
	     for (iy=1;iy<=NY;iy++){
                  for (ix=1;ix<=NX;ix++){
	
		       x = ix * DH;
		       y = iy * DH; 
		       r = sqrt(pow((x-(recpos[1][n]*DH)),2.0)+pow((y-(recpos[2][n]*DH)),2.0));
		       taper = 1.0;
		  
		       if(r<=R){
		         taper = f0 + (1.0 - f0) * pow((0.5 - 0.5 *cos(PI*r/R)),a);
		       }
		  
		           taper_coeff[iy][ix] = taper;
		           waveconv[iy][ix] *= taper;
		      
                  }
	     }
	}
	
	
	/*
        sprintf(modfile,"taper_coeff_%i.bin",ishot);
        writemod(modfile,taper_coeff,3); */		

	}

    /* ======================== */
    /* Read Taper from file     */  
    /* ======================== */
        
    if((sws>=4)&&(sws<=6)){
        
          if (MYID==0)
          {
                  fprintf(FP,"\n **Message from taper_grid (printed by PE %d):\n",MYID);
                  fprintf(FP," Coefficients for gradient taper are now calculated.\n");
          }
        
        fp_taper=fopen("taper.bin","r");

        /* loop over global grid */
        for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){
           
                fread(&grad_tap, sizeof(float), 1, fp_taper);                                                                                       
                taper_coeff[j][i]=grad_tap;
                             
            }
        }
 
        for (j=1;j<=NY;j++){   
           for (i=1;i<=NX;i++){     
              waveconv[j][i]*=taper_coeff[j][i];
           }                            
        }   
                                         
        fclose(fp_taper);
        
        sprintf(modfile,"taper_coeff_file.bin");
        writemod(modfile,taper_coeff,3);
                              
        }

        free_matrix(taper_coeff,1,NY,1,NX);
 
}



