/*---------------------------------------------------------------------------------------
 *  Read TTI model properties (c11, c13, c33, c44, theta, density) from external files  
 *
 *  D. Koehn
 *  Kiel, 04.02.2017
 *  ------------------------------------------------------------------------------------- */


#include "fd.h"

void readmod_elastic_TTI(float  **  rho, float **  c11, float **  c13, float **  c33, float **  c44, float ** theta){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern int WRITEMOD;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, vp, vs, epsi, delt, thetav;
	int i, j, ii, jj;
	FILE *fp_vs, *fp_vp, *fp_epsi, *fp_delt, *fp_theta, *fp_rho;
	char filename[STRING_SIZE];

	fprintf(FP,"\n...reading model information from external model-files...\n");
           
	   
	/* read density and TTI elastic tensor parameters */
        /* ---------------------------------------------- */

	fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
	sprintf(filename,"%s.vp",MFILE);
	fp_vp=fopen(filename,"r");
	if (fp_vp==NULL) err(" Could not open model file for Vp ! ");

	fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
	sprintf(filename,"%s.vs",MFILE);
	fp_vs=fopen(filename,"r");
	if (fp_vs==NULL) err(" Could not open model file for Vs ! ");

	fprintf(FP,"\t Epsilon:\n\t %s.epsi\n\n",MFILE);
	sprintf(filename,"%s.epsi",MFILE);
	fp_epsi=fopen(filename,"r");
	if (fp_epsi==NULL) err(" Could not open model file for Epsilon ! ");

	fprintf(FP,"\t Delta:\n\t %s.delt\n\n",MFILE);
	sprintf(filename,"%s.delt",MFILE);
	fp_delt=fopen(filename,"r");
	if (fp_delt==NULL) err(" Could not open model file for Delta ! ");

	fprintf(FP,"\t theta:\n\t %s.theta\n\n",MFILE);
	sprintf(filename,"%s.theta",MFILE);
	fp_theta=fopen(filename,"r");
	if (fp_theta==NULL) err(" Could not open model file for theta ! ");

	fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	sprintf(filename,"%s.rho",MFILE);
	fp_rho=fopen(filename,"r");
	if (fp_rho==NULL) err(" Could not open model file for densities ! ");
	   

	/* loop over global grid */
	for (i=1;i<=NXG;i++){
	    for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&epsi, sizeof(float), 1, fp_epsi);
			fread(&delt, sizeof(float), 1, fp_delt);
	        fread(&thetav, sizeof(float), 1, fp_theta);
			fread(&rhov, sizeof(float), 1, fp_rho);
				
		/* only the PE which belongs to the current global gridpoint 
		is saving model parameters in his local arrays */
		if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){

		    ii=i-POS[1]*NX;
		    jj=j-POS[2]*NY;
    		    
			c33[jj][ii]=rhov*vp*vp;
			c44[jj][ii]=rhov*vs*vs;		
			c11[jj][ii]=(2.*epsi+1.)*c33[jj][ii];
			c13[jj][ii]=sqrt((2.*c33[jj][ii]*c33[jj][ii]*delt+(c33[jj][ii]-c44[jj][ii])*(c11[jj][ii]+c33[jj][ii]-2.*c44[jj][ii]))/2.)-c44[jj][ii];

		    theta[jj][ii]=thetav * PI / 180.0;
            rho[jj][ii]=rhov;                                                                    
				
		}
	    }
	}
	
	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_epsi);
	fclose(fp_delt);
	fclose(fp_theta);
	fclose(fp_rho);

	
    /* each PE writes his model to disk and PE 0 merges model files */
		if(WRITEMOD){
	sprintf(filename,"%s.out.c11",MFILE);
	writemod(filename,c11,3);
	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.out.c13",MFILE);
        writemod(filename,c13,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                           
	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.out.c33",MFILE);
	writemod(filename,c33,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

        sprintf(filename,"%s.out.c44",MFILE);
	writemod(filename,c44,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);

	sprintf(filename,"%s.out.theta",MFILE);
	writemod(filename,theta,3);
	MPI_Barrier(MPI_COMM_WORLD);
	                        
	if (MYID==0) mergemod(filename,3);
		}

	/* clean up temporary files */
	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(filename,"%s.out.c11.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

	sprintf(filename,"%s.out.c13.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

	sprintf(filename,"%s.out.c33.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

	sprintf(filename,"%s.out.c44.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

	sprintf(filename,"%s.out.theta.%i.%i",MFILE,POS[1],POS[2]);
	remove(filename);

}




