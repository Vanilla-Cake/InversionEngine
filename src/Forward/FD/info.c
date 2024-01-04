/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *   last update 16/02/02
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," *******************************************************************************\n");
	fprintf(fp," This is a parallel 2-D acoustic / (an)isotropic elastic finite-difference modeling code \n");
	fprintf(fp," *******************************************************************************\n");
	fprintf(fp,"\n");

}
