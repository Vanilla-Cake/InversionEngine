/*  --------------------------------------------------------------------------
 *   Solving the (visco)-elastic 2D acoustic forward problem by finite-differences 
 *   for a single shot 
 *
 *   mode = 0 - forward modelling only, STF estimation or FWI gradient calculation
 *   mode = 1 - backpropagation of data residuals
 *   mode = 2 - evaluation objective function for step length estimation  
 * 
 *   
 *   D. Koehn
 *   Kiel, 10.06.2017
 *
 *  --------------------------------------------------------------------------*/

#include "fd.h"

void ac(struct waveAC *waveAC, struct waveAC_PML *waveAC_PML, struct matAC *matAC, struct mpiPSV *mpiPSV, 
         struct seisPSV *seisPSV, struct acq *acq, float *hc, int ishot, int nshots, int nsrc_loc, 
         int ns, int ntr, float **Ws, float **Wr, MPI_Request * req_send, MPI_Request * req_rec){

        /* global variables */
	extern float DT, DH, TSNAP1, TSNAP2, TSNAPINC;
	extern int MYID, MYID_SHOT, FDORDER, FW, SNAP_SHOT;
    extern int NX, NY, FREE_SURF, QUELLTYP;
	extern int POS[3], NDT, SEISMO;
    extern int SNAP, NT, NCOLORS;
	extern int FRE_SNAP;
	extern float freq0;
	extern FILE *FP;
	extern char SNAP_FILE[STRING_SIZE];
	extern MPI_Comm SHOT_COMM;

        /* local variables */
	int i,j,nt,lsamp,lsnap,nsnap, nd, hin1, imat, imat1, imat2, infoout;
    float tmp, tmp1, muss, lamss, tt;
	char filename[STRING_SIZE], outfiles[STRING_SIZE];

    nd = FDORDER/2 + 1;

	tt = 0.0;


	/*MPI_Barrier(MPI_COMM_WORLD);*/
        
	if (MYID_SHOT == 0)
	{
		fprintf(FP,"\n==================================================================================\n");
		fprintf(FP, "\n *****  Starting simulation (forward model) for shot %d of %d  ********** \n", ishot, nshots);
		fprintf(FP,"\n==================================================================================\n");
	}
			    
	/* initialize AC wavefields with zero */
	/*if (L){
		zero_denise_visc_PSV(-nd+1,NY+nd,-nd+1,NX+nd, (*wavePSV).pvx,(*wavePSV).pvy,(*wavePSV).psxx,(*wavePSV).psyy,(*wavePSV).psxy,
                                 (*wavePSV).ux,(*wavePSV).uy,(*wavePSV).uxy,(*wavePSV).pvxp1, (*wavePSV).pvyp1,(*wavePSV_PML).psi_sxx_x,(*wavePSV_PML).psi_sxy_x,
                                 (*wavePSV_PML).psi_vxx,(*wavePSV_PML).psi_vyx,(*wavePSV_PML).psi_syy_y,(*wavePSV_PML).psi_sxy_y,(*wavePSV_PML).psi_vyy,(*wavePSV_PML).psi_vxy,
                                 (*wavePSV_PML).psi_vxxs,(*wavePSV).pr,(*wavePSV).pp,(*wavePSV).pq);
	}else{*/	
		zero_denise_acoustic_AC(-nd+1,NY+nd,-nd+1,NX+nd,(*waveAC).pvx,(*waveAC).pvy,(*waveAC).p,
                            (*waveAC).ux,(*waveAC).pvxp1, (*waveAC).pvyp1,(*waveAC_PML).psi_p_x,
                            (*waveAC_PML).psi_vxx,(*waveAC_PML).psi_p_y,(*waveAC_PML).psi_vyy,(*waveAC_PML).psi_vxxs);	
		if(FRE_SNAP)
		{
			for (i=1;i<=NX;i++){
		   		for (j=1;j<=NY;j++){
					   (*waveAC).P_real[j][i] = 0.0;
					   (*waveAC).P_imag[j][i] = 0.0;
				   }
		   }
		}
		
	/*} */                                                        
	     
	/*----------------------  loop over timesteps (forward model) ------------------*/

	lsnap=iround(TSNAP1/DT);  
	lsamp=NDT;
	nsnap=0;

	for (nt=1;nt<=NT;nt++){     
		        
		/* Check if simulation is still stable */
		/*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
		if (isnan((*waveAC).pvy[NY/2][NX/2])) {
		   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,(*waveAC).pvy[NY/2][NX/2]);
		   err(" Simulation is unstable !");}

	   infoout = !(nt%10000);

	   if (MYID==0){
	      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
	      /*time3=MPI_Wtime();*/
	   }

	      /* update of particle velocities */
		update_v_PML_AC(1, NX, 1, NY, nt, (*waveAC).pvx, (*waveAC).pvxp1, (*waveAC).pvxm1, (*waveAC).pvy, (*waveAC).pvyp1, (*waveAC).pvym1, (*waveAC).uttx, (*waveAC).utty, (*waveAC).p,       
						(*matAC).prip, (*matAC).prjp, (*acq).srcpos_loc,(*acq).signals,(*acq).signals,nsrc_loc,(*waveAC_PML).absorb_coeff,hc,infoout, 0, (*waveAC_PML).K_x, (*waveAC_PML).a_x, 
						(*waveAC_PML).b_x, (*waveAC_PML).K_x_half, (*waveAC_PML).a_x_half, (*waveAC_PML).b_x_half, (*waveAC_PML).K_y, (*waveAC_PML).a_y, (*waveAC_PML).b_y, (*waveAC_PML).K_y_half, 
						(*waveAC_PML).a_y_half, (*waveAC_PML).b_y_half, (*waveAC_PML).psi_p_x, (*waveAC_PML).psi_p_y);
           
		/* exchange of particle velocities between PEs */
		exchange_v_AC((*waveAC).pvx,(*waveAC).pvy, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);
		                                                       
	   	update_s_acoustic_PML_AC(1, NX, 1, NY, (*waveAC).pvx, (*waveAC).pvy, (*waveAC).ux, (*waveAC).p, (*matAC).ppi, 
                                     (*waveAC_PML).absorb_coeff, (*matAC).prho, hc, infoout, (*waveAC_PML).K_x, (*waveAC_PML).a_x, (*waveAC_PML).b_x, (*waveAC_PML).K_x_half, (*waveAC_PML).a_x_half, 
                                     (*waveAC_PML).b_x_half, (*waveAC_PML).K_y, (*waveAC_PML).a_y, (*waveAC_PML).b_y, (*waveAC_PML).K_y_half, (*waveAC_PML).a_y_half, (*waveAC_PML).b_y_half, (*waveAC_PML).psi_vxx,  
                                     (*waveAC_PML).psi_vyy, 0);  


	   /* explosive source */
	   if (QUELLTYP==1){
	   
		psource_AC(nt,(*waveAC).p,(*acq).srcpos_loc,(*acq).signals,nsrc_loc,0);

           }

 	 
	   if ((FREE_SURF) && (POS[2]==0)){
	   	/* if (L) */   /* visco-acoustic */
			/*surface_visc_PML_PSV(1, (*wavePSV).pvx, (*wavePSV).pvy, (*wavePSV).psxx, (*wavePSV).psyy, (*wavePSV).psxy, (*wavePSV).pp, (*wavePSV).pq, (*matPSV).ppi, 
					    (*matPSV).pu, (*matPSV).prho, (*matPSV).ptaup, (*matPSV).ptaus, (*matPSV).etajm, (*matPSV).peta, hc, (*wavePSV_PML).K_x, (*wavePSV_PML).a_x, 
					    (*wavePSV_PML).b_x, (*wavePSV_PML).psi_vxxs);
		else */     /* acoustic */
	   		
                surface_acoustic_PML_AC(1, (*waveAC).p);
	   }

	   /* stress exchange between PEs */
	   exchange_p_AC((*waveAC).p, (*mpiPSV).bufferlef_to_rig, (*mpiPSV).bufferrig_to_lef, (*mpiPSV).buffertop_to_bot, (*mpiPSV).bufferbot_to_top, req_send, req_rec);

		/* store amplitudes at receivers in section-arrays */
		if (SEISMO){
			seismo_AC(nt, ntr, (*acq).recpos_loc, (*seisPSV).sectionvx, (*seisPSV).sectionvy, 
				(*seisPSV).sectionp, (*seisPSV).sectioncurl, (*seisPSV).sectiondiv, 
				(*waveAC).pvx, (*waveAC).pvy, (*waveAC).p, (*matAC).ppi, (*matAC).ppi, (*matAC).prho, hc);
			/*lsamp+=NDT;*/
		}

	   /* WRITE SNAPSHOTS TO DISK */
	   if ((SNAP) && (ishot == SNAP_SHOT) && (nt==lsnap) && (nt<=iround(TSNAP2/DT))){

	      snap_AC(FP,nt,++nsnap,(*waveAC).pvx,(*waveAC).pvy,(*waveAC).p,(*matAC).ppi,(*matAC).ppi,hc);

	      lsnap=lsnap+iround(TSNAPINC/DT);
	   }

	   /* WRITE FREQUENCY DOMAIN SNAPSHOTS TO DISK*/
	   if(FRE_SNAP)
	   {
		   tt = tt + DT;
		   for (i=1;i<=NX;i++){
		   		for (j=1;j<=NY;j++){
					   (*waveAC).P_real[j][i] = (*waveAC).P_real[j][i] + (*waveAC).p[j][i]*cos(2.*PI*freq0*tt);
					   (*waveAC).P_imag[j][i] = (*waveAC).P_imag[j][i] + (*waveAC).p[j][i]*sin(2.*PI*freq0*tt);
				   }
		   }

	   }


	   }/*--------------------  End  of loop over timesteps ----------*/	

		if(FRE_SNAP)
		{
			sprintf(filename,"%s.real.p.shot%i",SNAP_FILE,ishot);
			writemod(filename,(*waveAC).P_real,3);
			MPI_Barrier(SHOT_COMM);
			if (MYID_SHOT==0) {
				mergemod(filename,3);	
				sprintf(outfiles,"rm %s*.*",filename);
				printf("%s\n",outfiles);
				system(outfiles);
			}

			sprintf(filename,"%s.imag.p.shot%i",SNAP_FILE,ishot);
			writemod(filename,(*waveAC).P_imag,3);
			MPI_Barrier(SHOT_COMM);
			if (MYID_SHOT==0) {
				mergemod(filename,3);	
				sprintf(outfiles,"rm %s*.*",filename);
				printf("%s\n",outfiles);
				system(outfiles);
			}
			
		}


}
