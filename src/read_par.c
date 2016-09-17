/**/
/*------------------------------------------------------------------------
 *   Read FD-Parameters from Stdin                           
 *   last update 05/04/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

/* reading input-parameter from input-file or stdin for
viscoelastic finite-differnce modelling with fdveps */

#include "fd.h"

void read_par(FILE *fp_in){

/* declaration of extern variables */
extern int   NX, NY, SNAP, SNAP_FORMAT;
extern float DH;
extern int SEISMO, READMOD, READ_REC;
extern int LOG, TAPER, TAPERLENGTH;
extern float REFREC[4];
extern char  MFILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
extern char JACOBIAN[STRING_SIZE],DATA_DIR[STRING_SIZE];
extern int  MYID, IDX, IDY, NPML; 
extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT;
extern int GRAD_METHOD, NFREQ, STF_INV;
extern int MODEL_FILTER, FILT_SIZE;

extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
extern int SWS_TAPER_FILE;
extern float SRTRADIUS, EXP_TAPER_GRAD_HOR, A0_PML;
extern int SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
extern int MIN_ITER;
extern char INV_MODELFILE[STRING_SIZE];
extern int nfstart, nf;
extern int nfstart_jac, nf_jac;
extern float VPUPPERLIM, VPLOWERLIM;
extern char PARA[STRING_SIZE];

extern float FC_START, FC_END, FC_INCR, EPS_HESS;

extern int LNORM;

extern int LINESEARCH,STEPMAX;
extern float EPS_SCALE, SCALEFAC, C1, C2;

extern float PRO;

extern int TIMEWIN, NORMALIZE;
extern float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
extern char PICKS_FILE[STRING_SIZE];

extern char MISFIT_LOG_FILE[STRING_SIZE]; 


extern int N_STREAMER;
extern float REC_INCR_X, REC_INCR_Y;

extern int NLBFGS, PCG_BETA;

extern float VP0_1, VP0_2, DVP0, GRAD0_1, GRAD0_2, DGRAD0;
extern char GRIDSEARCH_FILE[STRING_SIZE];

/* definition of local variables */
char s[74];
int  c=0, lineno=0, l;

   while ((c=getc(fp_in)) != EOF){
      if ((c=='\n') && (getc(fp_in)!='#')){     
	 lineno++;
	 /* printf(" reading line %d \n",lineno);*/
	 switch (lineno){
	 case 1 : 
            fscanf(fp_in,"%s =%i",s,&INVMAT);
            break;
	 case 2 :
	    fscanf(fp_in,"%s =%i",s,&NX);
	    break;
	 case 3 :
	    fscanf(fp_in,"%s =%i",s,&NY);
	    break;
	 case 4 :
	    fscanf(fp_in,"%s =%f",s,&DH);
	    break;
	 case 5 :
	    fscanf(fp_in,"%s =%s",s,SOURCE_FILE);
	    break;
	 case 6 :
	    fscanf(fp_in,"%s =%i",s,&READMOD);
	    break;
	 case 7 :
	    fscanf(fp_in,"%s =%s",s,MFILE);
	    break;
	 case 8 :
	    fscanf(fp_in,"%s =%i",s,&NPML);
	    break;
	 case 9 :
	    fscanf(fp_in,"%s =%f",s,&A0_PML);
	    break;
	 case 10 :
	    fscanf(fp_in,"%s =%i",s,&SNAP);
	    break;
	 case 11 :
	    fscanf(fp_in,"%s =%i",s,&IDX);
	    break;
	 case 12 :
	    fscanf(fp_in,"%s =%i",s,&IDY);
	    break;
	 case 13 :
	    fscanf(fp_in,"%s =%i",s,&SNAP_FORMAT);
	    break;
	 case 14 :
	    fscanf(fp_in,"%s =%s",s,SNAP_FILE);
	    break;
	 case 15 :
	    fscanf(fp_in,"%s =%i",s,&SEISMO);
	    break;
	 case 16 :
	    fscanf(fp_in,"%s =%s",s,REC_FILE);
	    break;
	 case 17 :
	    fscanf(fp_in,"%s =%i",s,&READ_REC);
	    break;
	 case 18 :
	    fscanf(fp_in,"%s =%f ,%f",s,&REFREC[1],&REFREC[2]);
	    break;
	 case 19 :
	    fscanf(fp_in,"%s =%i",s,&N_STREAMER);
	    break;
	 case 20 :
	    fscanf(fp_in,"%s =%f",s,&REC_INCR_X);
	    break;
	 case 21 :
	    fscanf(fp_in,"%s =%f",s,&REC_INCR_Y);
	    break;
	 case 22 :
	   fscanf(fp_in,"%s =%s",s,PICKS_FILE);                         
            break;
	 case 23 :
	    fscanf(fp_in,"%s =%s",s,LOG_FILE);
	    break;     			
	 case 24 :
	    fscanf(fp_in,"%s =%i",s,&LOG);
	    break;
	 case 25 :
	    fscanf(fp_in,"%s =%i",s,&ITERMAX);
	    break;
         case 26 :
            fscanf(fp_in,"%s =%i",s,&STF_INV);
            break;
	 case 27 :
	    fscanf(fp_in,"%s =%s",s,JACOBIAN);
	    break;   
	 case 28 :
	    fscanf(fp_in,"%s =%s",s,DATA_DIR);
	    break;    
	 case 29 :
	    fscanf(fp_in,"%s =%i",s,&TAPER);
	    break;
	 case 30 :
	    fscanf(fp_in,"%s =%i",s,&TAPERLENGTH);
	    break; 
	 case 31 :
	    fscanf(fp_in,"%s =%f",s,&EPS_HESS);
	    break;        	
	 case 32 :
	    fscanf(fp_in,"%s =%i, %i, %i, %i",s,&GRADT1,&GRADT2,&GRADT3,&GRADT4);
	    break; 
	 case 33 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_GRAD_VERT);
            break;            	        
	 case 34 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_GRAD_HOR);
            break; 
         case 35 :
            fscanf(fp_in,"%s =%f",s,&EXP_TAPER_GRAD_HOR);
            break;
	 case 36 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_GRAD_SOURCES);
            break; 
	 case 37 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_CIRCULAR_PER_SHOT);
            break;    
	 case 38 :
            fscanf(fp_in,"%s =%i",s,&SRTSHAPE);
            break;   
	 case 39 :
            fscanf(fp_in,"%s =%f",s,&SRTRADIUS);
            break; 
	 case 40 :
            fscanf(fp_in,"%s =%i",s,&FILTSIZE);
            break;
         case 41 :
            fscanf(fp_in,"%s =%i",s,&SWS_TAPER_FILE);
            break;                          
	 case 42 :
            fscanf(fp_in,"%s =%s",s,INV_MODELFILE);
            break; 
	 case 43 :
            fscanf(fp_in,"%s =%i",s,&nfstart);
            break;  
	 case 44 :
            fscanf(fp_in,"%s =%i",s,&nf);
            break;   
	 case 45 :
            fscanf(fp_in,"%s =%i",s,&nfstart_jac);
            break;
	 case 46 :
            fscanf(fp_in,"%s =%i",s,&nf_jac);
            break;
	 case 47 :
            fscanf(fp_in,"%s =%f",s,&VPUPPERLIM);
            break;
	 case 48 :
            fscanf(fp_in,"%s =%f",s,&VPLOWERLIM);
            break;
	 case 49 :         
            fscanf(fp_in,"%s =%i",s,&GRAD_METHOD);                         
            break;
	 case 50 :         
            fscanf(fp_in,"%s =%i",s,&PCG_BETA);                         
            break;
	 case 51 :         
            fscanf(fp_in,"%s =%i",s,&NLBFGS);                         
            break;   
	 case 52 :         
            fscanf(fp_in,"%s =%i",s,&MODEL_FILTER);                         
            break;  
	 case 53 :         
            fscanf(fp_in,"%s =%i",s,&FILT_SIZE);                         
            break;
         case 54 :    
            fscanf(fp_in,"%s =%i",s,&LINESEARCH);
            break;
	 case 55 :
	   fscanf(fp_in,"%s =%f",s,&EPS_SCALE);                         
            break;
	 case 56 :
	    fscanf(fp_in,"%s =%f, %f",s,&C1,&C2);
	    break; 
	 case 57 :
	   fscanf(fp_in,"%s =%i",s,&STEPMAX);                         
            break;
	 case 58 :
	   fscanf(fp_in,"%s =%f",s,&SCALEFAC);                         
            break;
	 case 59 :
	   fscanf(fp_in,"%s =%s",s,&MISFIT_LOG_FILE);                         
            break; 
	 case 60 :
	   fscanf(fp_in,"%s =%i",s,&MIN_ITER);                         
            break; 
	 case 61 :
	    fscanf(fp_in,"%s =%f, %f, %f",s,&VP0_1,&VP0_2,&DVP0);
	    break;
         case 62 :
	    fscanf(fp_in,"%s =%f, %f, %f",s,&GRAD0_1,&GRAD0_2,&DGRAD0);
	    break;
	 case 63 :
	    fscanf(fp_in,"%s =%s",s,GRIDSEARCH_FILE);
	    break; 
	 default:
	    break;
	 }
	 }
      }

fclose(fp_in);

}
