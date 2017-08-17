/*------------------------------------------------------------------------------------
 *  Calculate step length which satisfies the Wolfe condition 
 *  from T. van Leeuwen: https://github.com/tleeuwen/Penalty-Method/
 *  (adapted from http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3)
 *  
 *  D. Koehn
 *  Kiel, 14.12.2015
 *  ----------------------------------------------------------------------------------*/

#include "fd.h"

float wolfels_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct PML_AC *PML_AC, struct matAC *matAC, float ** srcpos, int nshots, int ** recpos, int ntr, int iter, int nstage, float alpha, float L2){

	/* declaration of global variables */
        extern int NX, NY, GRAD_METHOD, STEPMAX, MYID;
        extern float EPS_SCALE, C1, C2, SCALEFAC;

        /* declaration of local variables */
        float normg, mu, nu, ft, ** Snp1;
        float ** gt, ** Vpnp1, g0s0, gts0, maxgrad, maxvp;     
	int lsiter, done;

        gt =  matrix(1,NY,1,NX);

        /* copy  -> gt */
	store_mat((*fwiAC).grad,gt,NX,NY);


        /* calculate initial step length */
        if(iter==1){

           /* normg = norm_matrix(grad,NX,NY);
           alpha = 1.0/normg;*/

           maxgrad = maximum_m((*fwiAC).Hgrad,NX,NY);
           maxvp = maximum_m((*matAC).vp,NX,NY);

           alpha = EPS_SCALE * maxvp/maxgrad;

        }else{

           alpha = 1.0;

        }

        lsiter = 0;
	done = 0;
	mu = 0.0;
        nu = 1e30;
        alpha = alpha/SCALEFAC;

        if(MYID==0){
           printf(" Estimate step length by Wolfe line search \n");
           printf(" ----------------------------------------- \n");
        }

	while(done!=1){

              if(nu < 1e30){
                 alpha = (nu + mu)/SCALEFAC;
              }else{
                 alpha = SCALEFAC*alpha;
              }
    
              if(lsiter < STEPMAX){

                 if(MYID==0){
                    printf("gradient evaluation %d ... \n",lsiter);
                    printf("alpha = %e ... \n",alpha);
                 }

          	 /* copy vp -> vp_old */
	  	 store_mat((*matAC).vp,(*fwiAC).vp_old,NX,NY);

          	 /* test Vp-update */
	  	 calc_mat_change_wolfe((*fwiAC).Hgrad,(*matAC).vp,(*fwiAC).vp_old,alpha,1);

		 ft = grad_obj_AC(fwiAC,waveAC,PML_AC,matAC,srcpos,nshots,recpos,ntr,iter,nstage);
 
	  	 /* copy vp_old -> vp */
	  	 store_mat((*fwiAC).vp_old,(*matAC).vp,NX,NY);

                 lsiter++;

              }else{

                 alpha = 0.0;
                 break;

              }
           
              g0s0=dotp_matrix(gt,(*fwiAC).Hgrad,NX,NY);
              gts0=dotp_matrix((*fwiAC).grad,(*fwiAC).Hgrad,NX,NY);

              if(MYID==0){
                 printf("Wolfe Conditions \n");
                 printf("---------------- \n");
                 printf("ft = %e \t <= L2  = %e \n",ft,L2);
                 printf("ft = %e \t <= L2 + C1*alpha*g0s0 = %e \n",ft,L2 + C1*alpha*g0s0);
                 printf("gts0 = %e \t >= C2*g0s0 = %e \n",gts0,C2*g0s0);
              }
 
              /* if (ft > (L2 + C1*alpha*g0s0)){
                  nu = alpha;
              }else if(gts0 < (C2*g0s0)){
                  mu = alpha;
              }else{
                  done = 1;
              }*/

              if(ft < L2){
		done = 1;
	      }

              if(ft > L2){
                alpha/=2.0;
	      }

        }

        if(MYID==0){
           printf("final alpha = %e \n",alpha);
        }

        free_matrix(gt,1,NY,1,NX);

return alpha;
                	    
}
