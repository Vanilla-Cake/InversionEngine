/*------------------------------------------------------------------------
 *  Read acoustic model properties (vp,density) from files  
 *
 *  D. Koehn
 *  Kiel, 10.06.2017
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod_AC(float  **  rho, float **  pi){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern int WRITEMOD;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, piv, vp;
	int i, j, ii, jj;
	FILE *fp_vp, *fp_rho;
	char filename[STRING_SIZE];


	fprintf(FP,"\n...reading model information from model-files...\n");
           
	/* read density and seismic velocities */
	/* ----------------------------------- */

	fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
	sprintf(filename,"%s.vp",MFILE);
	fp_vp=fopen(filename,"r");
	if (fp_vp==NULL) err(" Could not open model file for Vp ! ");

	fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	sprintf(filename,"%s.rho",MFILE);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) err(" Could not open model file for densities ! ");	   
	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			fread(&rhov, sizeof(float), 1, fp_rho);
				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
                                rho[jj][ii]=rhov;
                                pi[jj][ii]=vp;
				
				}
			}
		}
	

	fclose(fp_vp);
	fclose(fp_rho);
		
	/* each PE writes his model to disk */
	if(WRITEMOD){
	sprintf(filename,"%s.out.vp",MFILE);
	writemod(filename,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);	
	
	sprintf(filename,"%s.out.rho",MFILE);
	writemod(filename,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);
	}

	/* clean up temporary files */
	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(filename,"%s.out.vp.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

	sprintf(filename,"%s.out.rho.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);
}




