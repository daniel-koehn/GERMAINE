/*------------------------------------------------------------------------
 *   Output of model after each iteration and end of stage
 *   
 *   Daniel Koehn
 *   last update 22.08.2017
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void model_out_TE(struct matTE *matTE, int nstage, int stage_switch){

	/* global variables */
	extern int MYID;
	extern char INV_MODELFILE[STRING_SIZE];
	
	/* local variables */
        char modfile[STRING_SIZE];
	
	if(stage_switch==0){sprintf(modfile,"%s_sig.bin",INV_MODELFILE);}
	if(stage_switch==1){sprintf(modfile,"%s_sig_stage_%d.bin",INV_MODELFILE,nstage);}
	if(MYID==0){writemod_true(modfile,((*matTE).sigma),3);}

	if(stage_switch==0){sprintf(modfile,"%s_eps.bin",INV_MODELFILE);}
	if(stage_switch==1){sprintf(modfile,"%s_eps_stage_%d.bin",INV_MODELFILE,nstage);}
	if(MYID==0){writemod_true(modfile,((*matTE).epsilon),3);}

}

