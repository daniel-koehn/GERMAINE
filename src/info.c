/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *   last update 13.12.2015
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ***********************************************************\n");
	fprintf(fp," GERMAINE: Parallel 2D acoustic frequency domain            \n");
        fprintf(fp," finite-difference and full waveform inversion code         \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," FDFD solver and FWI code                                   \n");
	fprintf(fp," developed by D. Koehn and D. De Nil                        \n");
	fprintf(fp," Institute of Geosciences, Kiel University, Germany         \n\n");
	fprintf(fp," See COPYING file for copying and redistribution conditions.\n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");

}
