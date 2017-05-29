/*------------------------------------------------------------------------
 *  fd.h - include file for GERMAINE          
 *
 *  Daniel Koehn
 *  Kiel, 19.08.2016
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <complex.h>
//#include <suitesparse/umfpack.h>
#include <umfpack.h> 

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 96
#define STRING_SIZE 74
#define STRING_SIZE2 256
#define REQUEST_COUNT 4

/* declaration of data-structures */

/* acoustic wavefield variables */
struct waveAC{
   float freq, dfreq, omega2;
   float **pr, **pi; 
   float *precr, *preci;
   double *RHSr, *RHSi, *Ar, *Ai;
   int *irow, *icol;
} waveAC; 

/* acoustic FWI variables */
struct fwiAC{
   float  ** lam, ** vp_old, ** grad, ** Hgrad, ** gradm, **hess;
   float ** forwardr, ** forwardi;
   float * presr, * presi;
   float stfr, stfi;
} fwiAC;

/* acoustic material parameters */
struct matAC{
   float  ** vp, ** ivp2, **k2;	
} matAC; 

/* Acquisition geometry */
struct acq{
   int ** recpos;
   float ** srcpos;
} acq;

/* acoustic PML variables*/
struct PML_AC{
   float ** Ar, ** Ai, ** Br, ** Bi, ** Cr, ** Ci;
   float ** Axr, ** Axi, ** Byr, ** Byi;
   float * d_x, * d_y, * b_x, * b_y, * a_x, * a_y; 
} PML_AC;

/* elastic SH-wavefield variables */

/* elastic SH-FWI variables */
struct fwiSH{
   float  ** mu, ** vs_old, ** grad, ** Hgrad, ** gradm, **hess;
   float ** forwardr, ** forwardi;
   float * presr, * presi;
   float stfr, stfi;
} fwiSH;

/* elastic SH material parameters */
struct matSH{
   float  **rho, ** vs, ** ivs2, **k2;	
} matSH; 


/* declaration of acoustic functions */

void alloc_fwiAC(struct fwiAC *fwiAC, int ntr);

void alloc_matAC(struct matAC *matAC);

void alloc_seis_AC(struct waveAC *waveAC, int ntr);

void alloc_waveAC(struct waveAC *waveAC, struct PML_AC *PML_AC);

void apply_hess_AC(float ** grad, float ** hess);

void ass_grad_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct matAC *matAC, float ** grad_shot, float **srcpos,  int nshots, int **recpos, int ntr, int ishot);

float calc_res_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, int ntr, int ishot, int nstage, int nfreq);

void calc_seis_AC(struct waveAC *waveAC, int ** recpos, int ntr);

void extract_PCG_AC(float * PCG_old, float ** waveconv);

void free_seis_AC(struct waveAC *waveAC, int ntr);

void forward_AC(char *fileinp1);

void forward_shot_AC(struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int nstage, int nfreq);

void fwi_FD_AC(char *fileinp1);

float grad_obj_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage);

void hessian_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage);

void hessian_shin_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage);

void init_A_AC_9p_pml(struct PML_AC *PML_AC, struct matAC *matAC, struct waveAC *waveAC);

void init_mat_AC(struct waveAC *waveAC, struct matAC *matAC);

float obj_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage);

float parabolicls_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage, float alpha, float L2);

void RHS_source_AC(struct waveAC *waveAC, float ** srcpos, int ishot);

void RHS_source_AC_adj(struct waveAC *waveAC, struct fwiAC *fwiAC, int ** recpos, int ntr);

void RTM_FD_AC(char *fileinp1);

void RTM_AC_out(float ** Vp);

void store_PCG_AC(float * PCG_old, float ** waveconv);

void write_seis_AC(struct waveAC *waveAC, int ishot, int ntr, int nstage, int nfreq);

/* declaration of elastic SH functions */

void alloc_matSH(struct matSH *matSH);

void forward_SH(char *fileinp1);

void forward_shot_SH(struct waveAC *waveAC, struct PML_AC *PML_AC, struct matSH *matSH, float ** srcpos, int nshots, int ** recpos, int ntr, int nstage, int nfreq);

void init_A_SH_9p_pml(struct PML_AC *PML_AC, struct matSH *matSH, struct waveAC *waveAC);

void pml_pro_SH(struct PML_AC *PML_AC, struct waveAC *waveAC, struct matSH *matSH);

void readmod_SH(struct matSH *matSH);

/* declaration of general functions */

float calc_mat_change(float  **  waveconv, float **  pi, float **  pinp1, int iter, float eps_scale, int itest, int nfstart);

void calc_mat_change_wolfe(float  **  Hgrad, float **  vp, float **  vp_old, float eps_scale, int itest);

void calc_nonzero();

float calc_opt_step(float *  L2t, float * epst, int sws);

void calc_S(float **  Vp, float **  S);

void check_descent(float ** waveconv, float ** gradp, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, int iter);

void cp_grad_frame(float ** A);

void cp_vec(struct waveAC *waveAC, int * Ti, int * Tj, double * Tx, double * Tz);

void descent(float ** grad, float ** gradm);

float dotp(float * vec1, float *vec2, int n1, int n2, int sw);

float dotp_matrix(float ** A, float ** B, int NX, int NY);

void exchange_grad_MPI(float ** grad);

void gauss_filt(float ** waveconv);

void grid_search(float ** Vp, float ** S, float ** TT, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr);

void info_mem(FILE *fp, int NLBFGS_vec, int ntr);

void info(FILE *fp);

void init_grad(float ** A);

void init_MPIshot(int nshots);

void LBFGS(float ** Hgrad, float ** grad, float ** gradm, int iter, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, float * alpha_LBFGS, 
float **Vp, float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec);
           
double LU_decomp(double  **A, double *x, double *b,int n);

float minimum_m(float **mat, int nx, int ny);
float maximum_m(float **mat, int nx, int ny);

void model(float  ** Vp);

void model_gridsearch(float  ** Vp, float vp0, float grad0);

void model_out(float ** Vp, int iter);		  

float norm(float ** waveconv, int iter, int sws);

float norm1(float ** TT, float ** TTold);

float norm_matrix(float **A, int NX, int NY);

void note(FILE *fp);

void PCG(float * PCG_new, float * PCG_old, float * PCG_dir, int PCG_class);

void pml_pro(struct PML_AC *PML_AC, struct waveAC *waveAC);

void precond(float ** grad);

float readdsk(FILE *fp_in, int format);

void read_par(FILE *fp_in);

void read_par_inv(FILE *fp,int nstage,int stagemax);

void readmod(struct matAC *matAC);

int **receiver(FILE *fp, int *ntr, int ishot);

void rot_LBFGS_vec(float * y_LBFGS, float * s_LBFGS, int NLBFGS, int NLBFGS_vec);

float **sources(int *nsrc);

void solvelin(float  **AA, float *bb, float *x, int e, int method);

void store_mat(float ** A, float ** B, int n, int m);

void sum_grad_MPI(float ** grad);

void taper_grad(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int sws);

void taper_grad_hor(float ** grad);

void taper_grad_file(float ** grad);

void taper_grad_shot(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int ishot);

void taper_grad_shot1(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int ishot);

void tripd(float *d, float *e, float *b, int n);

void vec2mat(float **pr, float **pi, double *xr, double *xi);

float wolfels(float ** waveconv, float ** gradp, float ** Vp, float ** S, float ** TT, float ** lam, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr,
              int iter, float alpha, float L2);

void write_gridsearch(float L2, float vp0, float grad0, int count);

void write_par(FILE *fp);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float ** array, int format);

void writemod_true(char modfile[STRING_SIZE], float ** array, int format);

void writemod_vec(char modfile[STRING_SIZE], double * array, int format);

void zero_LBFGS(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float * r_LBFGS, 
                 float * alpha_LBFGS, float * beta_LBFGS, float * rho_LBFGS);

void zero_LBFGS1(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS);

/* utility functions */
void err(char err_text[]);
void warning(char warn_text[]);
double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
double *dvector(int nl, int nh);
float **fmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);

float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh); 
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, 
int ndh);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);
