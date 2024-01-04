/*  --------------------------------------------------------------------------
 *   Solving the (visco)-elastic 2D VTI-forward problem by finite-differences 
 *   for a single shot 
 *
 *   mode = 0 - forward modelling only, STF estimation or FWI gradient calculation
 *   mode = 1 - backpropagation of data residuals
 *   mode = 2 - evaluation objective function for step length estimation  
 * 
 *   
 *   D. Koehn
 *   Kiel, 02.02.2017
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void VTI(struct wavePSV *wavePSV, struct wavePSV_PML *wavePSV_PML, struct matVTI *matVTI, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, MPI_Request * req_send, MPI_Request * req_rec){

    /* global variables */
	extern float DT, DH, TSNAP1, TSNAP2, TSNAPINC;
	extern int MYID, MYID_SHOT, FDORDER, FW, SNAP_SHOT;
    extern int NX, NY, QUELLTYP;
    // extern int FREE_SURF,POS[3];
	extern int NDT, SEISMO;
    extern int SNAP, NT, NCOLORS;
	extern FILE *FP;
	extern int FRE_SNAP;
	extern float freq0;
	extern char SNAP_FILE[STRING_SIZE];
	extern MPI_Comm SHOT_COMM;

    /* local variables */
	int i,j,nt,lsamp,lsnap,nsnap, nd, hin1, imat, imat1, imat2, infoout;
    float tmp, tmp1, muss, lamss;
	float tt;
	char filename[STRING_SIZE], outfiles[STRING_SIZE];

    nd = FDORDER/2 + 1;

	tt = 0.0;
	/*MPI_Barrier(MPI_COMM_WORLD);*/
        
        if(MYID_SHOT==0){

		fprintf(FP,"\n==================================================================================\n");
		fprintf(FP, "\n *****  Starting simulation (forward model) for shot %d of %d  ********** \n", ishot, nshots);
		fprintf(FP,"\n==================================================================================\n");
		
		}

			    
	/* initialize PSV wavefields with zero */
		zero_denise_elast_PSV(-nd+1,NY+nd,-nd+1,NX+nd,(*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
                            (*wavePSV).ux,(*wavePSV).uy,(*wavePSV).uxy,(*wavePSV).pvxp1, (*wavePSV).pvyp1,(*wavePSV_PML).psi_sxx_x,
                            (*wavePSV_PML).psi_sxy_x,(*wavePSV_PML).psi_vxx,(*wavePSV_PML).psi_vyx,(*wavePSV_PML).psi_syy_y,(*wavePSV_PML).psi_sxy_y,
                            (*wavePSV_PML).psi_vyy,(*wavePSV_PML).psi_vxy,(*wavePSV_PML).psi_vxxs);	                                                        
	     
		if(FRE_SNAP)
		{
			for (i=1;i<=NX;i++){
				for (j=1;j<=NY;j++){
						(*wavePSV).vx_real[j][i] = 0.0;
						(*wavePSV).vx_imag[j][i] = 0.0;
						(*wavePSV).vy_real[j][i] = 0.0;
						(*wavePSV).vy_imag[j][i] = 0.0;
					}
			}
		}
	/*----------------------  loop over timesteps (forward model) ------------------*/

	lsnap=iround(TSNAP1/DT);  
	lsamp=NDT;
	nsnap=0;

	for (nt=1;nt<=NT;nt++){     
		        
		/* Check if simulation is still stable */
		/*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
		if (isnan((*wavePSV).pvy[NY/2][NX/2])) {
		   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,(*wavePSV).pvy[NY/2][NX/2]);
		   err(" Simulation is unstable !");}

	   infoout = !(nt%10000);

	   if (MYID_SHOT==0){
	      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
	      /*time3=MPI_Wtime();*/
	   }

        /* update of particle velocities */
		update_v_PML_PSV(1, NX, 1, NY, nt, (*wavePSV).pvx, (*wavePSV).pvxp1, (*wavePSV).pvxm1, (*wavePSV).pvy, (*wavePSV).pvyp1, (*wavePSV).pvym1, (*wavePSV).uttx, (*wavePSV).utty, (*wavePSV).psxx, (*wavePSV).psyy,       
						(*wavePSV).psxy, (*matVTI).prip, (*matVTI).prjp, (*acq).srcpos_loc,(*acq).signals,(*acq).signals,nsrc_loc,(*wavePSV_PML).absorb_coeff,hc,infoout, 0, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, 
						(*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, 
						(*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_sxx_x, (*wavePSV_PML).psi_syy_y, (*wavePSV_PML).psi_sxy_y, (*wavePSV_PML).psi_sxy_x);
	                 
		                                           
		/* exchange of particle velocities between PEs */
		exchange_v_PSV((*wavePSV).pvx,(*wavePSV).pvy, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);
		                                                       
	    update_s_elastic_PML_VTI(1, NX, 1, NY, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).ux, (*wavePSV).uy, (*wavePSV).uxy, (*wavePSV).uyx, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matVTI).c11, (*matVTI).c13, 
                                     (*matVTI).c33, (*matVTI).c44h, (*wavePSV_PML).absorb_coeff, hc, infoout, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).K_x_half, (*wavePSV_PML).a_x_half, 
                                     (*wavePSV_PML).b_x_half, (*wavePSV_PML).K_y, (*wavePSV_PML).a_y, (*wavePSV_PML).b_y, (*wavePSV_PML).K_y_half, (*wavePSV_PML).a_y_half, (*wavePSV_PML).b_y_half, (*wavePSV_PML).psi_vxx,  
                                     (*wavePSV_PML).psi_vyy, (*wavePSV_PML).psi_vxy, (*wavePSV_PML).psi_vyx, 0);  


	   /* explosive source */
	   if (QUELLTYP==1){
	   
	       psource(nt,(*wavePSV).psxx,(*wavePSV).psyy,(*acq).srcpos_loc,(*acq).signals,nsrc_loc,0);
	       
           }

	   /* moment tensor source */
	   if (QUELLTYP==5) 	
	   msource(nt,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,(*acq).srcpos_loc,(*acq).signals,nsrc_loc,0);

	//    if ((FREE_SURF) && (POS[2]==0)){
 
	//    	//    surface_elastic_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*matPSV).ppi, (*matPSV).pu, (*matPSV).prho, hc, 
	// 	// 			       (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs); 
	//    }

	   /* stress exchange between PEs */
	    exchange_s_PSV((*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy, 
	      (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, 
	      (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top,
	      req_send, req_rec);

		/* store amplitudes at receivers in section-arrays */
		if (SEISMO)
		{
			seismo_ssg_VTI(nt, ntr, (*acq).recpos_loc, (*seisPSV).sectionvx, (*seisPSV).sectionvy, 
				(*seisPSV).sectionp, (*seisPSV).sectioncurl, (*seisPSV).sectiondiv, 
				(*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, hc);
			/*lsamp+=NDT;*/
		}

	   /* WRITE SNAPSHOTS TO DISK */
	   if ((SNAP) && (ishot == SNAP_SHOT) && (nt==lsnap) && (nt<=iround(TSNAP2/DT))){

	      snap_VTI(FP,nt,++nsnap,(*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,hc);
	      lsnap=lsnap+iround(TSNAPINC/DT);

	   }

		/* WRITE FREQUENCY DOMAIN SNAPSHOTS TO DISK*/
	   if(FRE_SNAP)
	   {
		   tt = tt + DT;
		   for (i=1;i<=NX;i++){
		   		for (j=1;j<=NY;j++){
					   (*wavePSV).vx_real[j][i] = (*wavePSV).vx_real[j][i] + (*wavePSV).pvx[j][i]*cos(2.*PI*freq0*tt);
					   (*wavePSV).vx_imag[j][i] = (*wavePSV).vx_imag[j][i] + (*wavePSV).pvx[j][i]*sin(2.*PI*freq0*tt);
					   (*wavePSV).vy_real[j][i] = (*wavePSV).vy_real[j][i] + (*wavePSV).pvy[j][i]*cos(2.*PI*freq0*tt);
					   (*wavePSV).vy_imag[j][i] = (*wavePSV).vy_imag[j][i] + (*wavePSV).pvy[j][i]*sin(2.*PI*freq0*tt);
				   }
		   }

	   }

	   }/*--------------------  End  of loop over timesteps ----------*/		

		if(FRE_SNAP)
		{
			sprintf(filename,"%s.real.vx.shot%i",SNAP_FILE,ishot);
			writemod(filename,(*wavePSV).vx_real,3);
			MPI_Barrier(SHOT_COMM);
			if (MYID_SHOT==0) {
				mergemod(filename,3);	
				sprintf(outfiles,"rm %s*.*",filename);
				printf("%s\n",outfiles);
				system(outfiles);
			}

			sprintf(filename,"%s.imag.vx.shot%i",SNAP_FILE,ishot);
			writemod(filename,(*wavePSV).vx_imag,3);
			MPI_Barrier(SHOT_COMM);
			if (MYID_SHOT==0) {
				mergemod(filename,3);	
				sprintf(outfiles,"rm %s*.*",filename);
				printf("%s\n",outfiles);
				system(outfiles);
			}

			sprintf(filename,"%s.real.vy.shot%i",SNAP_FILE,ishot);
			writemod(filename,(*wavePSV).vy_real,3);
			MPI_Barrier(SHOT_COMM);
			if (MYID_SHOT==0) {
				mergemod(filename,3);	
				sprintf(outfiles,"rm %s*.*",filename);
				printf("%s\n",outfiles);
				system(outfiles);
			}

			sprintf(filename,"%s.imag.vy.shot%i",SNAP_FILE,ishot);
			writemod(filename,(*wavePSV).vy_imag,3);
			MPI_Barrier(SHOT_COMM);
			if (MYID_SHOT==0) {
				mergemod(filename,3);	
				sprintf(outfiles,"rm %s*.*",filename);
				printf("%s\n",outfiles);
				system(outfiles);
			}
			
		}
}
