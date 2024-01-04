/**/
/*------------------------------------------------------------------------
 *   Read FD-Parameters from Stdin                           
 *   last update 05/04/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

/* reading input-parameter from input-file or stdin for
viscoelastic finite-differnce modelling with fdveps */

#include "fd.h"

void read_par(FILE *fp_in){

/* declaration of extern variables */
extern int   NX, NY, FDORDER, MAXRELERROR, QUELLART, SNAP, SNAP_FORMAT, SNAP_SHOT;
extern float DH, TIME, DT, TS, DAMPING;
extern float FPML;
extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, READREC, SRCREC, RUN_MULTIPLE_SHOTS;
extern int LOG, FW, PHYSICS;
extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
extern int  NPROC, NPROCX, NPROCY, MYID, IDX, IDY; 
extern float FC_SPIKE_1, FC_SPIKE_2;

extern float npower, k_max_PML;

extern int N_STREAMER, TOPO;
extern float REC_INCR_X, REC_INCR_Y;

extern int WRITE_STF, WRITEMOD;

extern int FRE_SNAP;
extern float freq0;

/* definition of local variables */
char s[120];
int  c=0, lineno=0, l;
 
   while ((c=getc(fp_in)) != EOF){
      if ((c=='\n') && (getc(fp_in)!='#')){     
	 lineno++;
	 /* printf(" reading line %d \n",lineno);*/
	 switch (lineno){
     case 1 :
	    fscanf(fp_in,"%s =%i",s,&PHYSICS);
	    break;
	 case 2 :
	    fscanf(fp_in,"%s =%i",s,&NPROCX);
	    break;
	 case 3 :
	    fscanf(fp_in,"%s =%i",s,&NPROCY);
	    break;     		
	 case 4 :
	    fscanf(fp_in,"%s =%i",s,&FDORDER);
	    break;
	 case 5 :
	    fscanf(fp_in,"%s =%i",s,&MAXRELERROR);
	    break;
	 case 6 :
	    fscanf(fp_in,"%s =%i",s,&NX);
	    break;
	 case 7 :
	    fscanf(fp_in,"%s =%i",s,&NY);
	    break;
	 case 8 :
	    fscanf(fp_in,"%s =%f",s,&DH);
	    break;
	 case 9 :
	    fscanf(fp_in,"%s =%f",s,&TIME);
	    break;
	 case 10 :
	    fscanf(fp_in,"%s =%f",s,&DT);
	    break;
	 case 11 :
	    fscanf(fp_in,"%s =%i",s,&QUELLART);
	    break;
	 case 12 :
	    fscanf(fp_in,"%s =%s",s,SIGNAL_FILE);
	    break;
	 case 13 :
	    fscanf(fp_in,"%s =%f",s,&TS);
	    break;
	 case 14 :
	    fscanf(fp_in,"%s =%i",s,&SRCREC);
	    break;
	 case 15 :
	    fscanf(fp_in,"%s =%s",s,SOURCE_FILE);
	    break;
	 case 16 :
            fscanf(fp_in,"\n%s =%i",s,&RUN_MULTIPLE_SHOTS);
            break;
	 case 17 :
            fscanf(fp_in,"%s =%f",s,&FC_SPIKE_1);
            break;
     case 18 :
            fscanf(fp_in,"%s =%f",s,&FC_SPIKE_2);
            break;
	 case 19 :
            fscanf(fp_in,"%s =%i",s,&WRITE_STF);
            break;
	 case 20 :
	    fscanf(fp_in,"%s =%s",s,MFILE);
	    break;
	 case 21 :
            fscanf(fp_in,"%s =%i",s,&WRITEMOD);
            break;
	 case 22 :
	    fscanf(fp_in,"%s =%i",s,&FREE_SURF);
	    break;
	 case 23 :
	    fscanf(fp_in,"%s =%i",s,&FW);
	    if (FW<0) FW=0;
	    break;
	 case 24 :
	    fscanf(fp_in,"%s =%f",s,&DAMPING);
	    break;
	 case 25 :
            fscanf(fp_in,"%s =%f",s,&FPML);
            break;
	 case 26 :
            fscanf(fp_in,"%s =%f",s,&npower);
            break;
	 case 27 :
            fscanf(fp_in,"%s =%f",s,&k_max_PML);
            break;     			
	 case 28 :
	    fscanf(fp_in,"%s =%i",s,&SNAP);
	    break;
	 case 29 :
	    fscanf(fp_in,"%s =%i",s,&SNAP_SHOT);
	    break;
	 case 30 :
	    fscanf(fp_in,"%s =%f",s,&TSNAP1);
	    break;
	 case 31 :
	    fscanf(fp_in,"%s =%f",s,&TSNAP2);
	    break;
	 case 32 :
	    fscanf(fp_in,"%s =%f",s,&TSNAPINC);
	    break;
	 case 33 :
	    fscanf(fp_in,"%s =%i",s,&IDX);
	    break;
	 case 34 :
	    fscanf(fp_in,"%s =%i",s,&IDY);
	    break;
	 case 35 :
	    fscanf(fp_in,"%s =%i",s,&SNAP_FORMAT);
	    break;
	 case 36 :
	    fscanf(fp_in,"%s =%s",s,SNAP_FILE);
	    break;
	 case 37 :
	    fscanf(fp_in,"%s =%i",s,&SEISMO);
	    break;
	 case 38 :
	    fscanf(fp_in,"%s =%i",s,&READREC);
	    break;
	 case 39 :
	    fscanf(fp_in,"%s =%s",s,REC_FILE);
	    break;
	 case 40 :
	    fscanf(fp_in,"%s =%f ,%f",s,&REFREC[1],&REFREC[2]);
	    break;
	 case 41 :
	    fscanf(fp_in,"%s =%i",s,&N_STREAMER);
	    break;
	 case 42 :
	    fscanf(fp_in,"%s =%f",s,&REC_INCR_X);
	    break;
	 case 43 :
	    fscanf(fp_in,"%s =%f",s,&REC_INCR_Y);
	    break;
	 case 44 :
	    fscanf(fp_in,"%s =%i",s,&NDT);
	    break;
	 case 45 :
	    fscanf(fp_in,"%s =%i",s,&SEIS_FORMAT);
	    break;
	 case 46 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_VX);
	    break;
	 case 47 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_VY);
	    break;
	 case 48 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_CURL);
	    break;
	 case 49 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_DIV);
	    break;
	 case 50 :
	    fscanf(fp_in,"%s =%s",s,SEIS_FILE_P);
	    break;     			
	 case 51 :
	    fscanf(fp_in,"%s =%s",s,LOG_FILE);
	    break;     			
	 case 52 :
	    fscanf(fp_in,"%s =%i",s,&LOG);
	    break; 
	 case 53 :
	    fscanf(fp_in,"%s =%i",s,&TOPO);
	    break; 
	 case 54 :
		fscanf(fp_in,"%s =%i",s,&FRE_SNAP);
		break;
	 case 55 :
	 	fscanf(fp_in,"%s =%f",s,&freq0);
		 break;
	 default:
	    break;
	 }
	 }
      }

/* define total number of MPI processes*/
NPROC = NPROCX*NPROCY;

fclose(fp_in);

}
