/*------------------------------------------------------------------------
 * Non-linear (Preconditioned) Conjugate Gradient Method (NL-PCG)
 *
 * 
 * Daniel Koehn
 * Kiel, 14.12.2015
 * ----------------------------------------------------------------------*/

#include "fd.h"

void PCG(float ** Hgrad, float ** grad, int iter){

        /* global variables */
	extern int NX, NY, GRAD_METHOD, MYID;
	extern char JACOBIAN[STRING_SIZE];
	
        /* local variables */
	char jac[225], jac2[225];
	int i, j;
	float betaz, betan, gradlastiter, gradclastiter, betar, beta;
	FILE *fp, *fp1;
	
	/* calculate conjugate gradient direction, if iter > 1 (after Mora 1987) */
	/* --------------------------------------------------------------------- */
	if(GRAD_METHOD!=3){

	   if(iter>1){
	   
	      sprintf(jac,"%s_p.old",JACOBIAN);
	      fp=fopen(jac,"rb");

	      sprintf(jac2,"%s_c.old",JACOBIAN);
	      fp1=fopen(jac2,"rb");
	   
	      /* apply scalar product to obtain the coefficient beta */
	      betaz = 0.0;
	      betan = 0.0;
	      for (i=1;i<=NX;i++){
		   for (j=1;j<=NY;j++){
	   	  
		        fread(&gradlastiter,sizeof(float),1,fp);
		  
			/* Polak and Ribiere */
			betaz += grad[j][i] * (grad[j][i] - gradlastiter);
			betan += gradlastiter * gradlastiter;
		  
			/* Fletcher and Reeves */
			/*betaz += (grad[j][i]) * (gradg[j][i]);
			betan += (gradlastiter) * (gradglastiter);*/
		  
		   }
	      }
	     
	      /*printf("TEST: vor exchange (MYID=%d): beta = betaz/betan = %e/%e = %e\n",MYID,betaz,betan,betaz/betan);*/

	      /*beta = 0.0f;
	      if(betan !=0.0f) beta = betaz/betan;*/

	      beta = 0.0;
	      beta = betaz/betan;          

	      /* direction reset */
	      if(beta<0.0){beta = 0.0;}
	     
	      // printf("\n\nTEST: after exchange (MYID=%d): beta = %e / %e = %e\n",MYID,betaz,betan,beta);
	     
	      /*fseek(FP6,0,SEEK_SET);*/
	     
	      for (i=1;i<=NX;i++){
		   for (j=1;j<=NY;j++){	
	   
			fread(&gradclastiter,sizeof(float),1,fp1);
		        Hgrad[j][i] = grad[j][i] + gradclastiter * beta;
		  
		   }
	      }
	     
	      fclose(fp);
	   
	      if(iter>=2){fclose(fp1);}

	   } /* end of iter > 1 */

	   if (iter==1){
	       for (i=1;i<=NX;i++){
		    for (j=1;j<=NY;j++){      
		         Hgrad[j][i] = grad[j][i];
		    } 	
		  
	       }
	   }

	   /* output of conjugate gradient direction */
	   sprintf(jac2,"%s_c.old",JACOBIAN);
	   fp=fopen(jac2,"wb");

	   for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		     fwrite(&Hgrad[j][i],sizeof(float),1,fp);
		}
	   }

	   fclose(fp);

	} /* end of if GRAD_METHOD!=3*/

	/* output of preconditioned gradient */
	sprintf(jac,"%s_p.old",JACOBIAN);
	fp=fopen(jac,"wb");

	/* output of the preconditioned gradient */
	for (i=1;i<=NX;i++){
	     for (j=1;j<=NY;j++){
		  fwrite(&grad[j][i],sizeof(float),1,fp);
	     }
	}

	fclose(fp);

}
