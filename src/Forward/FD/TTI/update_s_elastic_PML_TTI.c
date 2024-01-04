/*  Update_s_elastic_PML_TTI
 *
 *  updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *  by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *  and second order accuracy in time for the elastic PSV TTI problem
 *   
 *  Daniel Koehn
 *  Kiel, 05.02.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_PML_TTI(int nx1, int nx2, int ny1, int ny2,
	float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
	float **   sxy, float ** d11,  float ** d13, float ** d15, float ** d15h, float ** d33, float ** d35, float ** d35h, 
        float ** d55h, float ** absorb_coeff, float *hc, int infoout, float * K_x, float * a_x, float * b_x, float * K_x_half, 
	float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx, int mode){


	int i,j, m, fdoh, h, h1;
	float  vxx, vyy, vxy, vyx;
	float  dhi;	
	extern float DT, DH;
	extern int MYID, FDORDER, FW;
        extern int FREE_SURF;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;
	
	

	dhi = DT/DH;
	fdoh = FDORDER/2;

	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_s (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}
	


	switch (FDORDER){

	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i]))*dhi; 

        /* left boundary */                                         
        if((POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
         }

	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
        }			       				
				
			      /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       

			}
		}
		break;

	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i]))*dhi; 

        /* left boundary */                                         
        if((POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
         }

	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
        }


			      /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       

			}
		}
		break;

	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i]))*dhi; 

        /* left boundary */                                         
        if((POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
         }

	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
        }


			    /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       

			}
		}
		break;

	case 8:

    for (j=ny1;j<=ny2;j++){
	for (i=nx1;i<=nx2;i++){

			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4]))*dhi;
			
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3]))*dhi;

                        vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i]))*dhi;

                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i]))*dhi; 

        /* left boundary */                                         
        if((POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
         }

	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
        }

 			      /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       
   }}
		break;

	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4])
				       + hc[5]*(vx[j][i+4]-vx[j][i-5])
				      )*dhi;
				      
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i])
				       + hc[5]*(vy[j+4][i]-vy[j-5][i])
				      )*dhi;
				      
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3])
				       + hc[5]*(vy[j][i+5]-vy[j][i-4])
				      )*dhi;
				      
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i])
				       + hc[5]*(vx[j+5][i]-vx[j-4][i])
				      )*dhi;
				      
	/* left boundary */                                         
        if((POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
         }

	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
        }

			       /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       
			
				      
			
			}
		}
		break;
		
	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
				vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4])
				       + hc[5]*(vx[j][i+4]-vx[j][i-5])
				       + hc[6]*(vx[j][i+5]-vx[j][i-6])
				      )*dhi;
				      
				vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i])
				       + hc[5]*(vy[j+4][i]-vy[j-5][i])
				       + hc[6]*(vy[j+5][i]-vy[j-6][i])
				      )*dhi;
				      
				vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
				       + hc[2]*(vy[j][i+2]-vy[j][i-1])
				       + hc[3]*(vy[j][i+3]-vy[j][i-2])
				       + hc[4]*(vy[j][i+4]-vy[j][i-3])
				       + hc[5]*(vy[j][i+5]-vy[j][i-4])
				       + hc[6]*(vy[j][i+6]-vy[j][i-5])
				      )*dhi;
				      
				vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
				       + hc[2]*(vx[j+2][i]-vx[j-1][i])
				       + hc[3]*(vx[j+3][i]-vx[j-2][i])
				       + hc[4]*(vx[j+4][i]-vx[j-3][i])
				       + hc[5]*(vx[j+5][i]-vx[j-4][i])
				       + hc[6]*(vx[j+6][i]-vx[j-5][i])
				      )*dhi;
				
	 /* left boundary */                                         
        if((POS[1]==0) && (i<=FW)){
                        
                        psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        vxx = vxx / K_x[i] + psi_vxx[j][i];

                        psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                        vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
         }

        /* right boundary */                                         
        if((POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
		
                        h1 = (i-nx2+2*FW);
                        h = i;
                        
                        psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 

                        /*psi_vyx[j][h] = b_x_half[h] * psi_vyx[j][h] + a_x_half[h] * vyx;
                        vyx = vyx / K_x_half[h] + psi_vyx[j][h];*/
                        
                        psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
         }

	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                        psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                        psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                     
                        vyy = vyy / K_y[j] + psi_vyy[j][i];
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
                        vyy = vyy / K_y[h1] + psi_vyy[h1][i];
                        
                        /*psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * vxy;
                        vxy = vxy / K_y_half[j] + psi_vxy[j][i];*/
                        
                        psi_vxy[h1][i] = b_y_half[h1] * psi_vxy[h1][i] + a_y_half[h1] * vxy;
			vxy = vxy / K_y_half[h1] + psi_vxy[h1][i];
        
        }

 			      /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       				
			}
		}
		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = 0.0;
				vyy = 0.0;
				vyx = 0.0;
				vxy = 0.0;
				for (m=1; m<=fdoh; m++) {
					vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
					vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
					vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
					vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
				}	
	
			       /* save time derivative of forward wavefield for gradient calculation */
			      /*if((mode==0)&&(GRAD_FORM==2)){
	                        ux[j][i] = (c11[j][i] * vxx + c13[j][i] * vyy) / DT;
				uy[j][i] = (c13[j][i] * vxx + c33[j][i] * vyy) / DT;
				uxy[j][i] = c44h[j][i] * (vyx+vxy) / DT;
			      } */			      			      			      
			      
			      sxx[j][i] += d11[j][i] * vxx + d13[j][i] * vyy + d15h[j][i] * (vyx+vxy);
			      syy[j][i] += d13[j][i] * vxx + d33[j][i] * vyy + d35h[j][i] * (vyx+vxy);
			      sxy[j][i] += d15[j][i] * vxx + d35[j][i] * vyy + d55h[j][i] * (vyx+vxy);			       

			}
		}
		break;
		
	} /* end of switch(FDORDER) */


	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
