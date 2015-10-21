/* 1D test file for the os3d subroutine */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "sundialstypes.h"   /* realtype, integertype, booleantype defination */
#include "pihm.h"            /* Data Model and Variable Declarations     */
#include "rt.h"              /* Data Model and Variable Declarations for chemical processes */

#define UNIT_C 1440
#define ZERO   1E-20
#define LINE_WIDTH 256
#define WORDS_LINE 20
#define WORD_WIDTH 40

void OS3D(realtype , realtype , Chem_Data *);


void masscheck (double* mass, Chem_Data * CD, int flag){

  /* A mass balance checker */
  /* flag indicate which height to use, either before TR or after TR */
  /* not mass only conserves for tracer species */

  int i, j;
  double temp, temp_mass;

  for (i = 0; i < CD->NumSpc; i ++)
    mass[i] = 0.0;

  for (j = 0; j < CD->NumSpc; j ++){
    temp_mass = 0.0;
    for ( i = 0; i < CD->NumVol; i ++){
      temp =  CD->Vcele[i].t_conc[j] * CD->Vcele[i].area * CD->Vcele[i].height_o * CD->Vcele[i].porosity;
      if (temp < -1)
        fprintf(stderr," !mass_check: %12.8f\t%12.8f\t%12.8f\t%12.8f\t%d\n", CD->Vcele[i].porosity, CD->Vcele[i].area , CD\
		->Vcele[i].t_conc[j], CD->Vcele[i].height_t, i);
      else
	temp_mass += temp;
    }
    mass[j] = temp_mass;
  }
}



int main(int argc, char **argv)
{
  if (argc !=6){
    fprintf(stderr, "%s <n> <length> <velocity> <start_time> <final_time>\n",argv[0]);
    fprintf(stderr, "length   : length of column\n");
    fprintf(stderr, "velocity : velocity of flow\n");    
    fprintf(stderr, "start    : time of start\n");
    fprintf(stderr, "final    : time to end\n");
    exit(1);
  }

  Chem_Data chData;
  Chem_Data *CD;
  CD = &chData;

  FILE* debug = fopen("logfile/os3ddebug.log","w");
  fclose(debug);
 
  int blocks, i, j;
  double t_length, length, velocity, timelps = 0, stepsize = 18, endtime;

  blocks   = atoi(argv[1]);

  fprintf(stderr, "Blocks: %d\n",blocks);

  t_length = atof(argv[2]);
  fprintf(stderr, "Total length: %6.4f\n", t_length);
  velocity = atof(argv[3]);
  fprintf(stderr, "Velocity: %6.4f\n", velocity);
  length   = t_length/blocks;

  
  CD->StartTime = atof(argv[4]);
  timelps = CD->StartTime*1440;
  CD->OutItv = 1;
  CD->Cementation = 1.0 ;
  CD->TVDFlg = 1;
  CD->NumVol = blocks;
  CD->NumSpc = 1;
 
  CD->chemtype = (species*) malloc(CD->NumSpc*sizeof(species));
  CD->Vcele = (vol_conc*) malloc( (CD->NumVol+2) * sizeof(vol_conc));
  CD->Flux  = (face*) malloc( 2* blocks *sizeof(face));
  for ( j = 0 ; j < CD->NumSpc; j ++){
    
    CD->chemtype[j].DiffCoe = 1E-20;
    CD->chemtype[j].DispCoe = 1;
  }
  endtime = atof(argv[5])*1440;
  fprintf(stderr, "StartTime = %6.4f\nEndTime = %6.4f\n",CD->StartTime, endtime);


  for (j = 0; j < (blocks -1); j ++){
    CD->Flux[j].nodeup = j+1;
    CD->Flux[j].nodelo = j+2;
    CD->Flux[j].nodeuu = j > 0 ? j : 0;
    CD->Flux[j].nodell = (j+3) <= blocks ? (j+3) : 0;
    CD->Flux[j].distance = length;
    CD->Flux[j].velocity = velocity;
    CD->Flux[j].s_area   = 1.0;
    CD->Flux[j].flux     = velocity * 1.0;
    CD->Flux[j].BC       = 0;
  }

  for (i = blocks - 1; i < 2 * (blocks -1); i ++){
    CD->Flux[i].nodeup = i+2 - (blocks -1);
    CD->Flux[i].nodelo = i+1 - (blocks -1);
    CD->Flux[i].nodeuu = (i+3 - (blocks -1)) <= blocks? (i+3 - (blocks -1)) : 0;
    CD->Flux[i].nodell = (i -   (blocks -1)) >= 1? (i- (blocks -1)) : 0;
    CD->Flux[i].distance = length;
    CD->Flux[i].velocity = -velocity;
    CD->Flux[i].s_area   = 1.0;
    CD->Flux[i].flux     = -velocity * 1.0;
    CD->Flux[i].BC       = 0;
  }

  i = 2*(blocks - 1);
  CD->Flux[i].nodeup = 1;
  CD->Flux[i].nodelo = CD->NumVol + 1;
  CD->Flux[i].distance = length;
  CD->Flux[i].velocity = -velocity;
  CD->Flux[i].s_area = 1.0;
  CD->Flux[i].flux   = -velocity * 1.0;
  CD->Flux[i].BC     = 1;
  
  i = 2*(blocks -1)+1;
  
  CD->Flux[i].nodeup = CD->NumVol;
  CD->Flux[i].nodelo = CD->NumVol + 2;
  CD->Flux[i].distance = length;
  CD->Flux[i].velocity = velocity;
  CD->Flux[i].s_area = 1.0;
  CD->Flux[i].flux   = velocity * 1.0;
  CD->Flux[i].BC     = 1;


  CD->NumFac = 2* blocks;
  CD->NumVol = blocks + 2;

  for ( i = 0 ; i < CD->NumFac ; i ++){
    fprintf(stderr, "%d\t%d\t%d\t%6.4f\t%6.4f\n", i , CD->Flux[i].nodeup, CD->Flux[i].nodelo,CD->Flux[i].distance, CD->Flux[i].velocity);
  }


  for (i = 0; i < CD->NumVol; i ++){
    CD->Vcele[i].index = i + 1;
    CD->Vcele[i].height_o = 1.0;
    CD->Vcele[i].height_t = 1.0;
    CD->Vcele[i].height_v = 1.0;
    CD->Vcele[i].area     = length * 1.0;
    CD->Vcele[i].porosity = 0.25;
    CD->Vcele[i].q        = 0.0;
    CD->Vcele[i].BC       = 0;
    CD->Vcele[i].t_conc = (double*) malloc ( CD->NumSpc * sizeof(double));
    CD->Vcele[i].p_conc = (double*) malloc ( CD->NumSpc * sizeof(double));
    for ( j = 0 ; j < CD->NumSpc; j ++){
      CD->Vcele[i].t_conc[j] = 0.0;
      CD->Vcele[i].p_conc[j] = 0.0;
    }
    fprintf(stderr,"%d %d %12.8f %12.8f %12.8f\n", i, i+1, velocity, CD->Vcele[i].area, CD->Vcele[i].t_conc[0]);
  }

  CD->Vcele[blocks  ].BC = 1;
  CD->Vcele[blocks+1].BC = 1;
  for ( j = 0 ; j < CD->NumSpc; j ++)
    CD->Vcele[blocks  ].t_conc[j] = 1.0;


  //  for (i = 0; i < CD->NumVol; i ++) fprintf(stderr, "Sus %6.4f\n", CD->Vcele[i].t_conc[0]);


  for ( i = 0 ; i < CD->NumFac; i ++){
    fprintf(stderr, "%d\t%d\t%d\t%d\t%6.4f\t%6.4f\n",CD->Flux[i].nodeup, CD->Flux[i].nodelo, CD->Flux[i].nodeuu, CD->Flux[i].nodell ,CD->Flux[i].distance, CD->Flux[i].velocity);

  }

  fprintf(stderr,"%d %d %d\n",CD->NumVol, CD->NumFac, CD->NumSpc);
  //  for ( i = 0 ; i < CD->NumVol; i ++)
  //  fprintf(stderr, "fail %6.4f\n", CD->Vcele[i].t_conc[0]);
  
  double * mass_t = (double*) malloc (CD->NumSpc * sizeof(double));
  double * mass_o = (double*) malloc (CD->NumSpc * sizeof(double));
  FILE  * outfile = fopen("logfile/testlog","w");
  while ( timelps <= endtime){
    masscheck(mass_o, CD, 0);
    OS3D(timelps, stepsize, CD);
    timelps += stepsize;
    
    masscheck(mass_t, CD, 1);
    if ( ((int)timelps % 60) == 0){
      fprintf(outfile, "%4.2f\t",timelps);
      for ( i = 0; i < CD->NumVol; i ++){
	for ( j = 0; j < CD->NumSpc; j ++)
	  fprintf(outfile, "%6.4f\t", CD->Vcele[i].t_conc[j]);
      }
      
      fprintf(outfile, "\n");
    }
    fprintf(stderr, "mass: %12.8f\t%12.8f\t%12.8f\t%12.8f%\n",mass_o[0],mass_t[0],mass_o[0]-mass_t[0], (mass_o[0]-mass_t[0])/mass_o[0]*100);
  }
  fclose(outfile);
  return (0);
  
}
