/*------------------------------------------------------------------------
 *  2D acoustic/ (an)isotropic elastic time domain modeling Code 
 *
 * 
 *  ---------------------------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */

int main(int argc, char **argv){
char * fileinp;
FILE *fpsrc;

/* Initialize MPI environment */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

setvbuf(stdout, NULL, _IONBF, 0);
		
/* print program name, version etc to stdout*/
if (MYID == 0) info(stdout);

/* read parameters from parameter-file (stdin) */
fileinp=argv[1];
FILEINP1=argv[2];
FP=fopen(fileinp,"r");
if(FP==NULL) {
	if (MYID == 0){
		printf("\n==================================================================\n");
		printf(" Cannot open input file %s \n",fileinp);
		printf("\n==================================================================\n\n");
                err(" --- ");
	}
}

/* read input file *.inp */
read_par(FP);
 
/* Init shot parallelization*/
COLOR = MYID / (NPROCX * NPROCY);
MPI_Comm_split(MPI_COMM_WORLD, COLOR, MYID, &SHOT_COMM);
MPI_Comm_rank(SHOT_COMM, &MYID_SHOT);

/* Init subdomain communication*/
MPI_Comm_split(MPI_COMM_WORLD, MYID_SHOT, MYID, &DOMAIN_COMM);

NCOLORS = NP / (NPROCX * NPROCY);

/*printf("MYID: %d \t COLOR: %d \n", MYID, COLOR);
printf("NP: %d \t NCOLORS: %d \n", NP, NCOLORS);
printf("NX: %d \t NY: %d \n", NPROCX, NPROCY);*/

MPI_Barrier(MPI_COMM_WORLD);

count_src();

/*printf("Number of shots %d \n", NSHOTS);*/

/* check if parameters for PHYSICS and MODE are correct */
check_mode_phys();

/* ---------------------------------------------------- */
/* Forward 2D PSV-problem */
/* ---------------------------------------------------- */
if(PHYSICS==1){
  FD_PSV();
}

/* ---------------------------------------------------- */
/* Forward 2D AC-problem */
/* ---------------------------------------------------- */
if(PHYSICS==2){
  FD_AC();
}

/* ---------------------------------------------------- */
/* Forward 2D PSV VTI-problem */
/* ---------------------------------------------------- */
if(PHYSICS==3){
  FD_VTI();
}

/* ---------------------------------------------------- */
/* Forward 2D PSV TTI-problem */
/* ---------------------------------------------------- */
if(PHYSICS==4){
  FD_TTI();
}


MPI_Comm_free(&SHOT_COMM);

MPI_Finalize();
return 0;	

}/*main*/
