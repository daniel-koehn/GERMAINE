/*------------------------------------------------------------------------
 *   globvar.h - global variables for RAJZEL   
 *   
 *   Daniel Koehn
 *   Kiel, 12.12.2015
 *  ----------------------------------------------------------------------*/

/* definition of global variables used in the finite difference programs*/
/* generally, for the names of the global variables
   uppercase letters are used */

float XS, YS, DH;
float REFREC[4]={0.0, 0.0, 0.0, 0.0};
int   SEISMO, NSRC=1, READMOD;
int   NX, NY, NT, QUELLART, QUELLTYP, SNAP, SNAP_FORMAT, LOG, NTRG, INFO;
int   NXNY, DC, NXG, NYG, IDX, IDY, NX0, NY0, READ_REC;
char  SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE];
char  MFILE[STRING_SIZE], REC_FILE[STRING_SIZE];
char  LOG_FILE[STRING_SIZE];
FILE *FP;

/* Mpi-variables */
int   NP, NPSP, NPROC, MYID, IENDX, IENDY;
int   NSHOT1, NSHOT2;
int   POS[3], INDEX[5];     
const int TAG1=1,TAG2=2, TAG3=3, TAG4=4, TAG5=5,TAG6=6;

/* TDFWI Code DENISE_elastic*/
char JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE];
int TAPER, TAPERLENGTH;
int GRADT1,GRADT2,GRADT3,GRADT4;
int ITERMAX, REC1, REC2, INVMAT;
int HESSIAN, GRAD_METHOD, NLBFGS;
float FC_HESS_START, FC_HESS_INC;
int MODEL_FILTER, FILT_SIZE;
float EPSILON, MUN, EPSILON_u, EPSILON_rho;
 
int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
int SWS_TAPER_FILE;
float SRTRADIUS, EXP_TAPER_GRAD_HOR;

int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
int INV_RHO_ITER, INV_VS_ITER, INV_VP_ITER;
int MIN_ITER, GRAD_FORM;

char INV_MODELFILE[STRING_SIZE];
int nfstart, nf;

int nfstart_jac, nf_jac;

float VPUPPERLIM, VPLOWERLIM;

int INV_STF, N_STF, N_STF_START;
float EPS_STF, OFFSETC_STF;
char PARA[STRING_SIZE];

int TIME_FILT, ORDER, SHOTINC;
float FC_START, FC_END, FC_INCR;

int LNORM;

int LINESEARCH, STEPMAX;
float EPS_SCALE, SCALEFAC, C1, C2;

float PRO;

int TIMEWIN, NORMALIZE;
float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
char PICKS_FILE[STRING_SIZE];
char MISFIT_LOG_FILE[STRING_SIZE]; 

int GRAD_FILTER, FILT_SIZE_GRAD, FILT_SIZE_GRAD1 ;

/* parameters for energy weighted preconditioning */
int EPRECOND;

/* parameters for FD Hessian calculation */
int NFREQ;
float FC_HESS_START, FC_HESS_INC, FC;

/* parameters for wavenumber domain  damping */
float WD_DAMP, WD_DAMP1;

/* parameters for offset muting */
float OFFSETC;
int OFFSET_MUTE;

/* parameters for envelope calculation */
int ENV;

/* parameters for towed streamer */
int N_STREAMER;
float REC_INCR_X, REC_INCR_Y;

/* grid search parameters */
float VP0_1, VP0_2, DVP0, GRAD0_1, GRAD0_2, DGRAD0;
char GRIDSEARCH_FILE[STRING_SIZE];

/* FDFD PML parameters */
int NPML;
float A0_PML, OMEGA_PML;

/* Free surface */
int FREE_SURF;

/* parameters for SuiteSparse */
int NONZERO;

/* frequency range for each stage */
float FC_low, FC_high;
int NF;

/* parameters for Hessian approximation according to Shin et al. (2001) */
float MAX_HESS, EPS_HESS;

/* parameters for Laplace modelling/FWI */
int LAPLACE;

/* source time function inversion via Wiener deconvolution */
int STF_INV;
