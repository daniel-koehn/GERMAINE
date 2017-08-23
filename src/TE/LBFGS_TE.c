/*----------------------------------------------------------------------------------------------------------
 *
 * Quasi-Newton limited memory - Broyden-Fletcher-Goldfarb-Shanno (l-BFGS) method (Nocedal & Wright, 2006)
 * 
 * 
 * Daniel Koehn
 * Kiel, 23.08.2017
 * ---------------------------------------------------------------------------------------------------------*/

#include "fd.h"

void LBFGS_TE(struct fwiTE *fwiTE, struct matTE *matTE, int iter, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, float * alpha_LBFGS, 
float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec){

        /* global variables */
	extern int NX, NY, IDX, IDY, MYID, LBFGS_RESET;
	extern float MAT1_NORM, MAT2_NORM;
	extern char JACOBIAN[STRING_SIZE];
	
        /* local variables */
	char jac[225], jac1[225];
	int i, j, k, h, h1, h2;
	float gradplastiter, beta;
	float gamma_LBFGS, sum_nom, sum_denom;
        float LBFGSTMP, LBFGSTMP1, LBFGSTMP2, LBFGSTMP3, modellastiter, norm_fac, norm_fac_u, norm_fac_rho;
        float beta_LBFGS_1, tmp;
        int ki, itershift, iter1;
	FILE *FP3, *FP4, *FP6, *FP5, *FP7;
	
        itershift = 1;

if(MYID==0){ /* Apply l-BFGS only on MPI process 0 */

/* normalize material parameters */
scale_grad((*matTE).sigma,1.0/MAT1_NORM,(*matTE).sigmar,NX,NY);
scale_grad((*matTE).epsilon,1.0/MAT2_NORM,(*matTE).epsilonr,NX,NY);

/* calculate H^-1 * gradm, using the L-BFGS method, if iter > 1          */
/* --------------------------------------------------------------------- */

/* if(iter>1){ */
if(LBFGS_RESET==0){

   /* load old models and gradients - Vp and store them in the LBFGS vectors */
   /* ----------------------------------------------------------------------- */

   /*iter1 = iter-itershift;*/ /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */
   h = NLBFGS_vec*(LBFGS_pointer-1) + 1; /* locate current initial position in LBFGS-vector */

   sprintf(jac,"%s_p_sigma_grad.old",JACOBIAN);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_sigma.old",JACOBIAN);
   FP7=fopen(jac1,"rb");
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS[h] = (*fwiTE).grad_sigma[j][i]-gradplastiter;

    	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS[h] = (*matTE).sigmar[j][i]-modellastiter;  
          
          h++;
          
       }
     }
     
     fclose(FP6);
     fclose(FP7);

   sprintf(jac,"%s_p_epsilon_grad.old",JACOBIAN);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_epsilon.old",JACOBIAN);
   FP7=fopen(jac1,"rb");
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS[h] = (*fwiTE).grad_epsilon[j][i]-gradplastiter;

    	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS[h] = (*matTE).epsilonr[j][i]-modellastiter;  
          
          h++;
          
       }
     }
     
     fclose(FP6);
     fclose(FP7);     

     /* calculate improved first guess Hessian gamma_LBFGS */
     h1 = NLBFGS_vec*(LBFGS_pointer-1) + 1;
     h2 = NLBFGS_vec*LBFGS_pointer; 
     
     sum_nom = dotp(y_LBFGS,s_LBFGS,h1,h2,0);
     sum_denom = dotp(y_LBFGS,y_LBFGS,h1,h2,0);
     gamma_LBFGS = sum_nom/sum_denom;

     printf("MYID = %d \n",MYID);     
     printf("gamma_LBFGS = %e \n",gamma_LBFGS);
         
     /* update variable rho for all LBFGS-iterations and all parameter classes*/
     for(k=1;k<=NLBFGS;k++){
          
        h1 = NLBFGS_vec*(k-1) + 1;
        h2 = NLBFGS_vec*k;
        sum_nom = dotp(y_LBFGS,s_LBFGS,h1,h2,0); 
	
	if(fabs(sum_nom)>0.0){
	  rho_LBFGS[k] = 1.0/sum_nom;
	}
	else{
	  rho_LBFGS[k] = 0.0;
	} 
	  
	if(MYID==0){                                                
	   printf("rho_LBFGS = %e of k = %d \n",rho_LBFGS[k],k);
        }
	                                                       
     }
     
     /* save q_LBFGS for all material parameters */    
        
      h=1;                                                                                                                                                     

     for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
          
	     q_LBFGS[h] = (*fwiTE).gradm_sigma[j][i];
	     h++;	   
	      
         }
     }

     for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
          
	     q_LBFGS[h] = (*fwiTE).gradm_epsilon[j][i];
	     h++;	   
	      
         }
     }

     /* update alpha_LBFGS and q_LBFGS */
     for(k=NLBFGS;k>=1;k--){
		
       h1 = NLBFGS_vec*(k-1) + 1;
       h2 = NLBFGS_vec*k;
       sum_nom = dotp(s_LBFGS,q_LBFGS,h1,h2,1);
       alpha_LBFGS[k] = rho_LBFGS[k] * sum_nom;
       
       /* update q for all material parameters */
       h = NLBFGS_vec*(k-1) + 1;
       for (i=1;i<=NLBFGS_vec;i++){
           q_LBFGS[i] = q_LBFGS[i] - alpha_LBFGS[k] * y_LBFGS[h];
           h++;
       }
     }
	 
       /* Multiply gradient with approximated Hessian */
       for (i=1;i<=NLBFGS_vec;i++){
           r_LBFGS[i] = gamma_LBFGS * q_LBFGS[i];
       }

     /* calculate H^-1 * gradm[j][i] */
     for(k=1;k<=NLBFGS;k++){
        
        h1 = NLBFGS_vec*(k-1) + 1;
        h2 = NLBFGS_vec*k;
        /* calculate beta_LBFGS*/   
        sum_nom = dotp(y_LBFGS,r_LBFGS,h1,h2,1);
        beta_LBFGS_1 = rho_LBFGS[k] * sum_nom;

        h = NLBFGS_vec*(k-1) + 1;
        for (i=1;i<=NLBFGS_vec;i++){
	   r_LBFGS[i] = r_LBFGS[i] + s_LBFGS[h]*(alpha_LBFGS[k]-beta_LBFGS_1);
	   h++;
        }
         
     }

     /* update gradients */
     h=1;
                                                                                 
     /* sigma */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            (*fwiTE).Hgrad_sigma[j][i] = r_LBFGS[h];
	    h++;
		  
        }
     }

     /* epsilon */
     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            (*fwiTE).Hgrad_epsilon[j][i] = r_LBFGS[h];
	    h++;
		  
        }
     }


}

/* if(iter==1){*/
if(LBFGS_RESET==1){

     for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
               
            (*fwiTE).Hgrad_sigma[j][i] = (*fwiTE).gradm_sigma[j][i];
            (*fwiTE).Hgrad_epsilon[j][i] = (*fwiTE).gradm_epsilon[j][i];
		  
        }
     }

     LBFGS_RESET=0;
}

cp_grad_frame((*fwiTE).Hgrad_sigma);
cp_grad_frame((*fwiTE).Hgrad_epsilon);

/* save old models Sigma */
/* --------------------- */

/* save old model */
sprintf(jac,"%s_p_sigma.old",JACOBIAN);
FP3=fopen(jac,"wb");

for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
        tmp = (*matTE).sigmar[j][i];
        fwrite(&tmp,sizeof(float),1,FP3);
    }
}
	
fclose(FP3);

/* save old gradient */
sprintf(jac,"%s_p_sigma_grad.old",JACOBIAN);
FP3=fopen(jac,"wb");

for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
        tmp = (*fwiTE).grad_sigma[j][i];
        fwrite(&tmp,sizeof(float),1,FP3);
    }
}
	
fclose(FP3);
	
/* save H^-1 * g approximation*/
sprintf(jac,"%s_c_sigma_grad.old",JACOBIAN);
FP3=fopen(jac,"wb");
	
for (i=1;i<=NX;i=i+IDX){   
    for (j=1;j<=NY;j=j+IDY){
        tmp = (*fwiTE).Hgrad_sigma[j][i];
        fwrite(&tmp,sizeof(float),1,FP3);
    }
}
        
fclose(FP3);	

/* save old models Epsilon */
/* ----------------------- */

/* save old model */
sprintf(jac,"%s_p_epsilon.old",JACOBIAN);
FP3=fopen(jac,"wb");

for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
        tmp = (*matTE).epsilonr[j][i];
        fwrite(&tmp,sizeof(float),1,FP3);
    }
}
	
fclose(FP3);

/* save old gradient */
sprintf(jac,"%s_p_epsilon_grad.old",JACOBIAN);
FP3=fopen(jac,"wb");

for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
        tmp = (*fwiTE).grad_epsilon[j][i];
        fwrite(&tmp,sizeof(float),1,FP3);
    }
}
	
fclose(FP3);
	
/* save H^-1 * g approximation*/
sprintf(jac,"%s_c_epsilon_grad.old",JACOBIAN);
FP3=fopen(jac,"wb");
	
for (i=1;i<=NX;i=i+IDX){   
    for (j=1;j<=NY;j=j+IDY){
        tmp = (*fwiTE).Hgrad_epsilon[j][i];
        fwrite(&tmp,sizeof(float),1,FP3);
    }
}
        
fclose(FP3);	


} /* end MYID==0*/

/* distribute gradient to all MPI processes */
MPI_Barrier(MPI_COMM_WORLD);
exchange_grad_MPI((*fwiTE).Hgrad_sigma);
exchange_grad_MPI((*fwiTE).Hgrad_epsilon);
	
}
