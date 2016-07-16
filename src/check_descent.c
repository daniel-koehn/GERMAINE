/*------------------------------------------------------------------------------------
 *  Check if search direction is a descent direction
 *  
 *  
 *  D. Koehn
 *  Kiel, 14.12.2015
 *  ----------------------------------------------------------------------------------*/

#include "fd.h"

void check_descent(float ** Hgrad, float ** grad, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, int iter){

	/* declaration of global variables */
        extern int NLBFGS, GRAD_METHOD, NX, NY, MYID;

        /* declaration of local variables */
        float nom, denom, p;        

        /* check if search direction is a descent direction */
        nom = dotp_matrix(Hgrad,grad,NX,NY);
        denom = dotp_matrix(grad,grad,NX,NY);

        p = -nom/denom;

        if(MYID==0){
           printf("p = %e \n",p);
        }

        if (p<0.0){
        
            if(GRAD_METHOD==2){

               if(MYID==0){
                  printf("Loss of descent, reset l-bfgs history \n");
               }

               zero_LBFGS1(NLBFGS, NLBFGS_vec, y_LBFGS, s_LBFGS);

            }

        }
                	    
}
