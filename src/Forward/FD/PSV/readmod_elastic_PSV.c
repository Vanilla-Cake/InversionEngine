/*------------------------------------------------------------------------
 *  Read elastic model properties (vp,vs,density) from files  
 *
 *  D. Koehn
 *  Kiel, 24.07.2016
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod_elastic_PSV(float  **  rho, float **  pi, float **  u, float **  eta_x, float **  eta_z, float **  lambda){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern int WRITEMOD;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;
	extern int TOPO;
		
	/* local variables */
	float rhov, vp, vs;
	float eta_xv,eta_zv,lambda_v;
	int i, j, ii, jj;
	FILE *fp_vs, *fp_vp, *fp_rho;
	FILE *fp_etax, *fp_etaz, *fp_lambda;
	char filename[STRING_SIZE];


	   fprintf(FP,"\n...reading model information from modell-files...\n");
           
	   /* read density and seismic velocities */
	   /* ----------------------------------- */

	   fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
	   sprintf(filename,"%s.vp",MFILE);
	   fp_vp=fopen(filename,"r");
	   if (fp_vp==NULL) err(" Could not open model file for Vp ! ");


	   fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
	   sprintf(filename,"%s.vs",MFILE);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open model file for Vs ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open model file for densities ! ");
	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&rhov, sizeof(float), 1, fp_rho);
				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
					u[jj][ii]=vs;
                    rho[jj][ii]=rhov;
                    pi[jj][ii]=vp;

					if(!TOPO){
						lambda[jj][ii] = rho[jj][ii] * ( pi[jj][ii] * pi[jj][ii] - 2.* u[jj][ii] * u[jj][ii] );
						eta_x[jj][ii] = rho[jj][ii] * pi[jj][ii] * pi[jj][ii];
						eta_z[jj][ii] = eta_x[jj][ii];
					}


				}
			}
		}

	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);
	
	
	if(TOPO)
	{
		// Jian Cao: Parameter-modified method for surface topography

		fprintf(FP,"\n...reading MODIFIED model parameter for implementing surface topography ...\n");

		fprintf(FP,"\t ETAX:\n\t %s.etax\n\n",MFILE);
		sprintf(filename,"%s.etax",MFILE);
		fp_etax=fopen(filename,"r");
		if (fp_etax==NULL) err(" Could not open model file for ETAX ! ");


		fprintf(FP,"\t ETAZ:\n\t %s.etaz\n\n",MFILE);
		sprintf(filename,"%s.etaz",MFILE);
		fp_etaz=fopen(filename,"r");
		if (fp_etaz==NULL) err(" Could not open model file for ETAZ ! ");

		fprintf(FP,"\t LAMBDA:\n\t %s.lambda\n\n",MFILE);
		sprintf(filename,"%s.lambda",MFILE);
		fp_lambda=fopen(filename,"r");
		if (fp_lambda==NULL) err(" Could not open model file for LAMBDA ! ");

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&eta_xv, sizeof(float), 1, fp_etax);
			fread(&eta_zv, sizeof(float), 1, fp_etaz);
			fread(&lambda_v, sizeof(float), 1, fp_lambda);
				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
					lambda[jj][ii]=lambda_v;
                    eta_x[jj][ii]=eta_xv;
                    eta_z[jj][ii]=eta_zv;

				}
			}
		}

		fclose(fp_etax);
		fclose(fp_etaz);
		fclose(fp_lambda);

	}
	
	/* each PE writes his model to disk */
        if(WRITEMOD){

	   sprintf(filename,"%s.out.vp",MFILE);
	   writemod(filename,pi,3);
 	   MPI_Barrier(MPI_COMM_WORLD);

	   if (MYID==0) mergemod(filename,3);
	
	   sprintf(filename,"%s.out.vs",MFILE);
           writemod(filename,u,3);
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

	sprintf(filename,"%s.out.vs.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

	sprintf(filename,"%s.out.rho.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);


}




