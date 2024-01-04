/*------------------------------------------------------------------------
 *  Check if parameters MODE and PHYSICS are reasonable
 *
 *  D. Koehn
 *  Kiel, 02.02.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void check_mode_phys(){

	/* declaration of extern variables */
        extern int MODE, PHYSICS, MYID;
	
        if (MYID==0){

		printf("\n **Message from check_model_phys (printed by PE %d):\n",MYID);
		printf("\n");
		printf(" -----------------------  Seismic modeling  ----------------------\n");		
		printf("\n\n");

		switch (PHYSICS){
		case 1 :
			printf(" PHYSICS=%d: Solve 2D isotropic elastic PSV problem.\n",PHYSICS);
			break;
		case 2 :
			printf(" PHYSICS=%d: Solve 2D acoustic problem.\n",PHYSICS);
			break;
		case 3 :
			printf(" PHYSICS=%d: Solve 2D elastic PSV VTI problem.\n",PHYSICS);
			break;
		case 4 :
			printf(" PHYSICS=%d: Solve 2D elastic PSV TTI problem.\n",PHYSICS);
			break;
		case 5 :
			printf(" PHYSICS=%d: Solve 2D isotropic elastic SH problem.\n",PHYSICS);
			break;

		default:
			err(" PHYSICS unkown ");
		}

		printf("\n\n");

        }
}
