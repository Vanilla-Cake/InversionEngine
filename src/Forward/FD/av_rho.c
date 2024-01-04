/* $Id: av_rho.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/

#include "fd.h"

void av_rho(float **rho, float **rip, float **rjp){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern int TOPO;
	extern FILE *FP;
	extern char  MFILE[STRING_SIZE];


	/* local variables */
	float rjp_v, rip_v;
	int i, j, ii, jj;
	FILE *fp_rhoi, *fp_rhoj;
	char filename[STRING_SIZE];
		

	if(TOPO){

		// Jian Cao: Parameter-modified method for surface topography

		fprintf(FP,"\n...reading MODIFIED density parameter for implementing surface topography ...\n");

		fprintf(FP,"\t DENSITY_I:\n\t %s.rhoi\n\n",MFILE);
		sprintf(filename,"%s.rhoi",MFILE);
		fp_rhoi=fopen(filename,"r");
		if (fp_rhoi==NULL) err(" Could not open model file for DENSITY_I ! ");


		fprintf(FP,"\t DENSITY_J:\n\t %s.rhoj\n\n",MFILE);
		sprintf(filename,"%s.rhoj",MFILE);
		fp_rhoj=fopen(filename,"r");
		if (fp_rhoj==NULL) err(" Could not open model file for DENSITY_J ! ");

		/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&rip_v, sizeof(float), 1, fp_rhoi);
			fread(&rjp_v, sizeof(float), 1, fp_rhoj);
				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
                    rip[jj][ii]=rip_v;
                    rjp[jj][ii]=rjp_v;

				}
			}
		}

		fclose(fp_rhoi);
		fclose(fp_rhoj);


	}
	else
	{
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){
			
				rjp[j][i] = 1.0/(0.5*(rho[j][i]+rho[j+1][i])); 	
				rip[j][i] = 1.0/(0.5*(rho[j][i]+rho[j][i+1])); 	

				if((rho[j][i]<1e-4)&&(rho[j+1][i]<1e-4)){
				rjp[j][i] = 0.0;
				}

				if((rho[j][i]<1e-4)&&(rho[j][i+1]<1e-4)){
				rip[j][i] = 0.0;
				}	
		
		
			}
		}

	}
	



}
