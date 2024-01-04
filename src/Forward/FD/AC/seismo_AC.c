/*------------------------------------------------------------------------
 *   store amplitudes (particle velocities or pressure or curl and div) 
 *    at receiver positions in arrays
 *   last update 27/12/01, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo_AC(int lsamp, int ntr, int **recpos, float **sectionvx, 
float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
float **vx, float **vy, float **p, float **pi, float **u, float **rho, float *hc){ 
		
	extern int NDT, SEISMO, FDORDER;	
	extern float DH, DT;
	int i,j, itr, ins, nxrec, nyrec, m, fdoh;
	float dh24, dhi, vxx, vyy, vxy, vyx, muss, lamss;


	dh24=1.0/(DH*24.0);
	dhi = DT/DH;
	fdoh = FDORDER/2;

	/*ins=lsamp/NDT;*/
	ins=lsamp;
	for (itr=1;itr<=ntr;itr++){
		nxrec=recpos[1][itr];
		nyrec=recpos[2][itr];
		switch (SEISMO){
		case 1 : 
			sectionvx[itr][ins]=vx[nyrec][nxrec];
			sectionvy[itr][ins]=vy[nyrec][nxrec];
			break;
		
		case 2 : 
			sectionp[itr][ins]=p[nyrec][nxrec];
			break;
		
		case 3 :				
			i=nxrec; j=nyrec;
			
			vxx = 0;
			vyy = 0;
			vyx = 0;
			vxy = 0;
			for (m=1; m<=fdoh; m++) {
				vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
				vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
				vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
				vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
			}
			vxx *= dhi;
			vyy *= dhi;
			vyx *= dhi;
			vxy *= dhi;
			
			muss = u[j][i] * u[j][i] * rho[j][i];
			lamss = pi[j][i] * pi[j][i] * rho[j][i];
			sectiondiv[itr][ins]=(vxx+vyy)*sqrt(lamss);
			sectioncurl[itr][ins]=(vxy-vyx)*sqrt(muss);
			
			break;
		
		case 4 :				
			i=nxrec; j=nyrec;

			vxx = 0;
			vyy = 0;
			vyx = 0;
			vxy = 0;
			for (m=1; m<=fdoh; m++) {
				vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
				vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
				vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
				vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
			}
			vxx *= dhi;
			vyy *= dhi;
			vyx *= dhi;
			vxy *= dhi;

			muss = u[j][i] * u[j][i] * rho[j][i];
			lamss = pi[j][i] * pi[j][i] * rho[j][i];
			sectiondiv[itr][ins]=(vxx+vyy)*sqrt(lamss);
			sectioncurl[itr][ins]=(vxy-vyx)*sqrt(muss);

			
			sectionvx[itr][ins]=vx[nyrec][nxrec];
			sectionvy[itr][ins]=vy[nyrec][nxrec];			
			sectionp[itr][ins]=p[nyrec][nxrec];
			break;

		}

	}
}
