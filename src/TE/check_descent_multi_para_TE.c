/*------------------------------------------------------------------------------------
 *  Check if search direction is a descent direction
 *  
 *  
 *  D. Koehn
 *  Kiel, 23.08.2017
 *  ----------------------------------------------------------------------------------*/

#include "fd.h"

void check_descent_multi_para_TE(struct fwiTE *fwiTE, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float * r_LBFGS, 
                                 float * alpha_LBFGS, float * beta_LBFGS, float * rho_LBFGS){

	/* declaration of global variables */
        extern int NLBFGS, GRAD_METHOD, NX, NY, MYID, LBFGS_RESET;

        /* declaration of local variables */
        float nom, denom, p;        

        /* check if search direction is a descent direction */
	nom = 0.0;
        denom = 0.0;

        nom = dotp_matrix((*fwiTE).Hgrad_sigma,(*fwiTE).grad_sigma,NX,NY);
        nom += dotp_matrix((*fwiTE).Hgrad_epsilon,(*fwiTE).grad_epsilon,NX,NY);

        denom = dotp_matrix((*fwiTE).grad_sigma,(*fwiTE).grad_sigma,NX,NY);
        denom += dotp_matrix((*fwiTE).grad_epsilon,(*fwiTE).grad_epsilon,NX,NY);

        p = -nom/denom;

        if(MYID==0){
           printf("p = %e \n",p);
        }

        if (p<0.0){
        
            if(GRAD_METHOD==2){

               if(MYID==0){
                  printf("Loss of descent, reset l-bfgs history \n");
               }

	       LBFGS_RESET = 1;
               //zero_LBFGS1(NLBFGS, NLBFGS_vec, y_LBFGS, s_LBFGS);
	       zero_LBFGS(NLBFGS, NLBFGS_vec, y_LBFGS, s_LBFGS, q_LBFGS, r_LBFGS, alpha_LBFGS, beta_LBFGS, rho_LBFGS);

            }

        }
                	    
}
