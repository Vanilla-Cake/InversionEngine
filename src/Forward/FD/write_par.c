/*------------------------------------------------------------------------
 *  Write Modeling parameters                           
 *
 *  D. Koehn
 *  Kiel, 02.08.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp){

	/* declaration of extern variables */
	extern int   NX, NY, NT, QUELLART, FDORDER, MAXRELERROR;
	extern int  SNAP, SNAP_FORMAT, SNAP_SHOT, SRCREC, TAPER;
	extern float DH, TIME, DT, TS, DAMPING, FPML, npower, k_max_PML;
	extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, FW;
	extern int  READREC;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
	extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE];
	extern char SIGNAL_FILE[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
	extern char  MFILE[STRING_SIZE];
	extern int NP, NPROCX, NPROCY, MYID;
	extern int N_STREAMER, TOPO, PHYSICS;
    extern float REC_INCR_X, REC_INCR_Y;

	extern float FC_SPIKE_1, FC_SPIKE_2;
	
	
	/* definition of local variables */
	int l;
	

	fprintf(fp,"\n **Message from write_par (printed by PE %d):\n",MYID);
	fprintf(fp,"\n");
	fprintf(fp,"------------------------- Processors ------------------------\n");
	fprintf(fp," Number of PEs in x-direction (NPROCX): %d\n",NPROCX);
	fprintf(fp," Number of PEs in vertical direction (NPROCY): %d\n",NPROCY);
	fprintf(fp," Total number of PEs in use: %d\n",NP);
	fprintf(fp,"\n");
	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", NX);
	fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", NY);
	fprintf(fp," Grid-spacing (DH): %e meter\n", DH);
	fprintf(fp," Time of wave propagation (T): %e seconds\n",TIME);
	fprintf(fp," Timestep (DT): %e seconds\n", DT);
	fprintf(fp," Number of timesteps: %i \n",NT);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- FD ORDER -----------------------------\n");
	fprintf(fp," FDORDER = %d\n",FDORDER);
	fprintf(fp," MAXRELERROR = %d\n",MAXRELERROR);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");

	if (SRCREC){
		fprintf(fp," reading source positions, time delay, centre frequency \n");
		fprintf(fp," and initial amplitude from ASCII-file \n");
		fprintf(fp,"\t%s\n\n",SOURCE_FILE);}

	fprintf(fp," wavelet of source:");

	switch (QUELLART){
	case 1 :
		fprintf(fp," Ricker\n");
		break;
	case 2 :
		fprintf(fp," Fuchs-Mueller\n");
		break;
	case 3 :
		fprintf(fp," reading from \n\t %s\n",SIGNAL_FILE);
		break;
	case 4 :
		fprintf(fp," sinus raised to the power of 3.0 \n");
		break;
	case 5 :
		fprintf(fp," 1st derivative of Gaussian \n");
		break;
	case 6 :
	        fprintf(fp," Klauder\n");
	        fprintf(fp," FC_SPIKE_1 = %f, FC_SPIKE_2 = %f, TS = %f \n",FC_SPIKE_1,FC_SPIKE_2,TS);
	        break;			                                	
	default :
		err(" Sorry, incorrect specification of source wavelet ! ");
	}

        fprintf(fp,"\n\n");

	/*fprintf(fp," Type of source:");
	switch (QUELLTYP){
	case 1 :
		fprintf(fp," explosive source \n");
		break;
	case 2 :
		fprintf(fp," point source with directive force in x-direction\n");
		break;
	case 3 :
		fprintf(fp," point source with directive force in (vertical) y-direction\n");
		break;
	case 4 :
		fprintf(fp," rotated point source with directive force in x- and y-direction\n");
		break;
        case 5 :
                fprintf(fp," source defined by moment tensor\n");
                break;	                                 	
	default :
		err(" Sorry, wrong source type specification ! ");
	}
	
	fprintf(fp,"\n");*/

	if (SEISMO){
		fprintf(fp," ------------------------- RECEIVER  --------------------------\n");
		if (READREC){
                        if(READREC==1){
			    fprintf(fp," reading receiver positions from single file \n");
			    fprintf(fp,"\t%s.dat\n\n",REC_FILE);
			}
                        if(READREC==2){
			    fprintf(fp," reading receiver positions from multiple files \n");
			    fprintf(fp,"\t%s_shot_X.dat\n\n",REC_FILE);
			}
			fprintf(fp," reference_point_for_receiver_coordinate_system:\n");
			fprintf(fp," x=%f \ty=%f\t z=%f\n",REFREC[1], REFREC[2], REFREC[3]);
		}

                if (N_STREAMER){
                        fprintf(fp," ------------------------- Towed streamer  --------------------------\n");
                        fprintf(fp," Assuming the first N_STREAMER = %d receivers belong to a streamer \n",N_STREAMER);
                        fprintf(fp," Shifting the streamer by ... \n");
                        fprintf(fp," REC_INCR_X=%f in x-direction \n",REC_INCR_X);
                        fprintf(fp," REC_INCR_Y=%f in y-direction \n",REC_INCR_Y);
                }
	}

	fprintf(fp," ------------------------- FREE SURFACE ------------------------\n");
	if (FREE_SURF) fprintf(fp," free surface at the top of the model ! \n");
	else fprintf(fp," no free surface at the top of the model ! \n");
	fprintf(fp,"\n");

	fprintf(fp," ------------------------- CPML ---------------------\n");
	if (FW>0.0){
		fprintf(fp," width of absorbing frame is %i gridpoints.\n",FW);
		fprintf(fp," CPML damping applied. \n");
		fprintf(fp," Damping velocity in the PML frame in m/s: %f .\n",DAMPING);
		fprintf(fp," Frequency within the PML frame in Hz: %f \n",FPML);
		fprintf(fp," npower: %f \n",npower);
		fprintf(fp," k_max: %f \n",k_max_PML); 
	}
	else fprintf(fp," absorbing frame not installed ! \n");

	fprintf(fp," ------------------------- MODEL-FILES -------------------------\n");
	fprintf(fp," names of model-files: \n");
	fprintf(fp,"\t shear wave velocities:\n\t %s.vs\n",MFILE);
	// fprintf(fp,"\t tau for shear waves:\n\t %s.ts\n",MFILE);
	fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
	fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
	// fprintf(fp,"\t tau for P-waves:\n\t %s.tp\n",MFILE);

	if (SNAP){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of");
		switch(SNAP){
		case 1:
			fprintf(fp," x- and y-component");
			fprintf(fp," of particle velocity.\n");
			break;
		case 2:
			fprintf(fp," pressure field.\n");
			break;
		case 3:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			break;
		case 4:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			fprintf(fp," x- and y-component of particle velocity.\n");
			break;
		default:
			err(" sorry, incorrect value for SNAP ! \n");
		}

		fprintf(fp," \t write snapshots for shot SNAP_SHOT= %d \n", SNAP_SHOT);
		fprintf(fp," \t first (TSNAP1)= %8.5f s\n", TSNAP1);
		fprintf(fp," \t last (TSNAP2)=%8.5f s\n",TSNAP2);
		fprintf(fp," \t increment (TSNAPINC) =%8.5f s\n\n",TSNAPINC);
		fprintf(fp," \t first_and_last_horizontal(x)_gridpoint = %i, %i \n",1,NX);
		fprintf(fp," \t first_and_last_vertical_gridpoint = %i, %i \n",1,NY);
		fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);
		switch (SNAP_FORMAT){
		case 1 :
			err(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
			break;
		default:
			err(" Don't know the format for the Snapshot-data ! \n");
		}
	
		fprintf(fp,"\n\n");
	}
	if (SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
		if ((SEISMO==1) || (SEISMO==4)){
			fprintf(fp," seismograms of ");
			fprintf(fp," x- and y-component");
			fprintf(fp," of particle velocity.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s\n\t%s\n",SEIS_FILE_VX,SEIS_FILE_VY);
		}
		if ((SEISMO==2) || (SEISMO==4)){
			fprintf(fp," seismograms of pressure field (hydrophones).\n");
			fprintf(fp," output-file: \n ");
			fprintf(fp,"\t%s\n",SEIS_FILE_P);
		}
		if ((SEISMO==3) || (SEISMO==4)){
			fprintf(fp," seismograms of curl and div.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s\n\t%s\n",SEIS_FILE_CURL,SEIS_FILE_DIV);
			
		}		
	
		switch (SEIS_FORMAT){
		case 1 :
			fprintf(fp," The data is written in IEEE SU-format . \n");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary IEEE (4 byte per float)");
			break;
		default:
			err(" Sorry. I don't know the format for the seismic data ! \n");
		}
		fprintf(fp," samplingrate of seismic data: %f s\n",NDT*DT);
		fprintf(fp," Number of samples per trace: %i \n", iround(NT/NDT));
		fprintf(fp," ----------------------------------------------------------\n");
		fprintf(fp,"\n");
		fprintf(fp,"\n");
	}

	if(TOPO) {
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  Modeling with surface topography  ----------------------\n");
		if(PHYSICS != 1 ) err(" Sorry. It is only avaliable to isotropic elastic media. ! \n");
		fprintf(fp,"\n"); 
		if(!FREE_SURF) err(" FREE_SURF should be 1, because free-surface boundary conditions will be applied. ! \n");
		
	}
                                                   
	fprintf(fp,"\n");
	fprintf(fp," **************************************************************\n");
	fprintf(fp,"\n");
}
