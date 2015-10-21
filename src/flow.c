/* 3D flow simulator for validating the React Transport */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "pihm.h"            /* Data Model and Variable Declarations     */
#include "sundialstypes.h"   /* realtype, integertype, booleantype defination */
#include "rt.h"              /* Data Model and Variable Declarations for chemical processes */
#include "file_rt.h"


#define UNIT_C 1440
#define ZERO   1E-20
#define LINE_WIDTH 512
#define MIN_STEP 1E-5
#define min(a,b) (((a)<(b))? (a):(b))

void OS3D(realtype , realtype , Chem_Data *);
int  React(realtype, realtype, Chem_Data *, int, int *);
void ReactCon(realtype, realtype, Chem_Data *, int, int *);
void ReactControl(realtype, realtype, Chem_Data *, int);
void Lookup(FILE* ,Chem_Data *, int);
void Speciation(Chem_Data *, int);
int  SpeciationType(FILE*, char *);
void NewCell(vol_conc* , int , int );
void CopyCell(vol_conc* , vol_conc*);
void masscheck(double*, Chem_Data * CD, int flag);
void AdptTime(Chem_Data * , realtype *, double, double,  double * , int);


int Coordinates(int i, int j, int k, int x_grid, int y_grid, int z_grid){
  // This subroutine finds the index of the cell i, j, k along axis x, y, z;
  int layer = x_grid * y_grid;
  int row   = x_grid;
  int total = x_grid * y_grid * z_grid;
  int index;
  index = i-1 + (j-1) * row + (k-1) * layer;
  if ( i < 0 ) { index = -1; return (index);} // non-existence
  else if ( i == 0 ) {index = total; return (index);}
  if ( i > x_grid+1 ) {index = -1; return (index);}
  else if ( i == x_grid+1) {index = total + 1; return (index);}
  if ( j < 0 ) {index = -1; return (index);}
  else if ( j == 0) {index = total +2;return (index);}
  if ( j > y_grid+1 ) {index = -1;return (index);}
  else if ( j == y_grid+1) {index = total +3;return (index);}
  if ( k < 0 ) {index = -1;return (index);}
  else if ( k == 0) {index = total +4;return (index);}
  if ( k > z_grid+1 ) {index = -1;return (index);}
  else if ( k == z_grid+1) {index = total + 5;return (index);}

  return (index);
}

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
      if (flag == 0)
	temp = CD->Vcele[i].t_conc[j] * CD->Vcele[i].area * CD->Vcele[i].height_o * CD->Vcele[i].porosity;
      if (flag == 1)
	temp = CD->Vcele[i].t_conc[j] * CD->Vcele[i].area * CD->Vcele[i].height_t * CD->Vcele[i].porosity;
      if (temp < -1)
	fprintf(stderr," mass_check: %12.8f\t%12.8f\t%12.8f\t%12.8f\t%d\n", CD->Vcele[i].porosity, CD->Vcele[i].area , CD->Vcele[i].t_conc[j], CD->Vcele[i].height_t, i);
	      else 
      temp_mass += temp;
    }
    mass[j] = temp_mass;
  }
}

void NewCell(vol_conc* newcell, int NumStc, int NumSsc){

  /* Allocate memory for the vol_conc elements */
  /* Elemental operation for future reuse */

  newcell->t_conc   = (double*) malloc(NumStc *sizeof(double));
  newcell->p_conc   = (double*) malloc(NumStc *sizeof(double));
  newcell->p_para   = (double*) malloc(NumStc *sizeof(double));
  newcell->p_actv   = (double*) malloc(NumStc *sizeof(double));
  newcell->p_type   = (int *  ) malloc(NumStc *sizeof(int   ));
  newcell->s_conc   = (double*) malloc(NumSsc *sizeof(double));
  newcell->s_actv   = (double*) malloc(NumSsc *sizeof(double));
}

void CopyCell(vol_conc* fromcell, vol_conc* tocell){
  int i, j, NumStc = fromcell->NumStc, NumSsc = fromcell->NumSsc;
  tocell->index = fromcell->index;
  tocell->i_condition = fromcell->i_condition;
  tocell->BC    = fromcell->BC;
  for ( i = 0; i < NumStc; i ++){
    tocell->t_conc[i] = fromcell->t_conc[i];
    tocell->p_conc[i] = fromcell->p_conc[i];
    tocell->p_actv[i] = fromcell->p_actv[i];
    tocell->p_para[i] = fromcell->p_para[i];
    tocell->p_type[i] = fromcell->p_type[i];
  }
  for ( i = 0; i < NumSsc; i ++){
    tocell->s_conc[i] = fromcell->s_conc[i];
    tocell->s_actv[i] = fromcell->s_actv[i];
  }
  tocell->vol_real = fromcell->vol_real;
  tocell->sat    = fromcell->sat;
  tocell->height_o    = fromcell->height_o;
  tocell->height_t    = fromcell->height_t;
  tocell->height_v    = fromcell->height_v;
  tocell->area    = fromcell->area;
  tocell->porosity    = fromcell->porosity;
  tocell->q    = fromcell->q;
  tocell->temperature    = fromcell->temperature;
  tocell->NumStc = fromcell->NumStc;
  tocell->NumSsc = fromcell->NumSsc;
}


int main(int argc, char **argv)
{
  
  if (argc !=2){
    fprintf(stderr, " %s <filename>\n", argv[0]);
    fprintf(stderr, " Please place input files under input directory and type in filename\n");
    exit(1);
  }
  
  char * filename = (char*) malloc (10*sizeof(char));
  if (argc == 2)
    strcpy(filename, argv[1]);
  

  fprintf(stderr, " Input file is %s\n",filename);

  int i, j, k, num_face =0, num_blocks, num_species, num_mineral, num_ads, num_cex, num_other, num_conditions = 0, line_width = LINE_WIDTH, 
    words_line = WORDS_LINE, word_width = WORD_WIDTH, Global_diff = 0, Global_disp = 0, error_flag = 0,
    speciation_flg = 0;
  double total_flux = 0.0, total_area = 0.0,     tmpval[WORDS_LINE], peclet, rt_step, temp_rt_step;


  char keyword[WORD_WIDTH], line[256], word[WORD_WIDTH];
  char ** tmpstr = (char**) malloc(words_line*sizeof(char*));


  for ( i = 0; i < words_line; i ++)
    tmpstr[i] = (char*) malloc(word_width*sizeof(char));


  Chem_Data chData;
  Chem_Data *CD;
  CD = &chData;

  CD->Keq = (double*) malloc( CD->NumSsc * sizeof(double));
  CD->KeqKinect = (double*) malloc( CD->NumMin * sizeof(double));


  char* gridfn   = (char*)malloc  ((strlen(filename)+12)*sizeof(char));
  sprintf(gridfn, "input/%s.grid",filename);
  FILE* gridfile = fopen(gridfn,"r");
  char* chemfn   = (char*)malloc  ((strlen(filename)+12)*sizeof(char));
  sprintf(chemfn, "input/%s.chem",filename);
  FILE* chemfile = fopen(chemfn,"r");
  char* datafn   = (char*) malloc ((strlen(filename)+12)*sizeof(char));
  sprintf(datafn, "input/%s.cdbs",filename);
  FILE* database = fopen(datafn,"r");

  int x_blocks = 0, y_blocks = 0, z_blocks = 0;
  double x_length = 0.0 , y_length =0.0 ,  z_length = 0.0, x_vel =0.0 , y_vel= 0.0, z_vel=0.0, v_cell=0.0;
  double starttime=0.0, endtime=0.0, stepsize=0.0;

  assert ( chemfile != NULL);
  assert ( database != NULL);
  assert ( gridfile != NULL);
  fprintf(stderr, " Input file grid is %s\n",gridfn);
  fprintf(stderr, " Input file chem is %s\n",chemfn);
  fprintf(stderr, " Input file cdbs is %s\n",datafn);

  rewind(gridfile);
  fgets(line, line_width, gridfile);
  //  while ( keymatch(line, "END", tmpval, tmpstr) !=1)
  while( !feof(gridfile))
    {
      fgets(line, line_width, gridfile);
      if ( keymatch(line, "GRID",tmpval,tmpstr) == 1){
	x_blocks = (int)tmpval[0];
	y_blocks = (int)tmpval[1];
	z_blocks = (int)tmpval[2];
	num_blocks = x_blocks*y_blocks*z_blocks;
	fprintf(stderr, " Total blocks in the system = %d*%d*%d = %d\n",x_blocks, y_blocks, z_blocks, num_blocks);
      }

      if ( keymatch(line, "LENGTH", tmpval, tmpstr) ==1){
	x_length = tmpval[0];
	y_length = tmpval[1];
	z_length = tmpval[2];
	v_cell   = x_length*y_length*z_length;
	fprintf(stderr, " Total volume of the system = %6.4f*%6.4f*%6.4f*%d = %8.4f\n", x_length, y_length, z_length, num_blocks, v_cell*num_blocks);
      }

    }


  // NumVol would equal to num_blocks + 6 as there is 6 boundary layers for 3D cartesian coordinates.
  /*  while ( keymatch(line, "LENGTH", tmpval, tmpstr) !=1)
    fgets(line, line_width, gridfile);
  x_length = tmpval[0];
  y_length = tmpval[1];
  z_length = tmpval[2];
  v_cell   = x_length*y_length*z_length;
  fprintf(stderr, " Total volume of the system = %6.4f*%6.4f*%6.4f*%d = %8.4f\n", x_length, y_length, z_length, num_blocks, v_cell*num_blocks);*/
  rewind(gridfile);
  fgets(line, line_width, gridfile);
  while ( keymatch(line, "VELOCITY", tmpval, tmpstr) !=1)
    fgets(line, line_width, gridfile);
  x_vel = tmpval[0];
  y_vel = tmpval[1];
  z_vel = tmpval[2];
  fprintf(stderr, " Velocity stated for x, y and z axis: %6.4f\t %6.4f\t %6.4f\n", x_vel, y_vel, z_vel);
  rewind(gridfile);
  fgets(line, line_width, gridfile);
  while ( keymatch(line, "SPAN", tmpval, tmpstr) !=1)
    fgets(line, line_width, gridfile);
  starttime = tmpval[0];
  endtime   = tmpval[1];
  fprintf(stderr, " Span of simulation is from %6.4f to %6.4f\n", starttime, endtime);
  rewind(gridfile);
  fgets(line, line_width, gridfile);
  while ( keymatch(line, "POROSITY", tmpval, tmpstr) !=1)
    fgets(line, line_width, gridfile);
  double porosity = tmpval[0];
  fprintf(stderr, " Constant porosity is %6.4f\n",porosity);

  CD->NumVol = num_blocks + 6;
  CD->NumOsv = num_blocks;
  CD->SPCFlg = 1;
  CD->PrpFlg = 0;
  /* index num_blocks + 0/1 x boundary layer;
   * index num_blocks + 2/3 y boundary layer;
   * index num_blocks + 4/5 z boundary layer;
   * note that condition index starts from 1;
   * +1 to convert cell index into condition index;
   */

  double timelps, start_step;
  timelps = UNIT_C * starttime;
  endtime = UNIT_C * endtime;
  stepsize = 100; // 1 min
  start_step = MIN_STEP;  // systematic default minimum time step.
  
  CD->CptFlg = 0;
  CD->TimLst = 0;
  
  rewind(chemfile);
  fgets(line, line_width, chemfile);
  while ( keymatch(line, "RUNTIME", tmpval, tmpstr)!=1 )
    fgets(line, line_width, chemfile);
  while (keymatch(line, "END", tmpval, tmpstr)!=1){
    fgets(line, line_width, chemfile);
    if( keymatch(line, "tvd", tmpval, tmpstr) == 1){
      if ( strcmp(tmpstr[1],"false")==0) CD->TVDFlg = 0;
      if ( strcmp(tmpstr[1],"true" )==0) CD->TVDFlg = 1;
      if ( strcmp(tmpstr[1],"false")&& strcmp(tmpstr[1],"true")) fprintf(stderr,"TVD FLAG INPUT ERROR!\n");
      fprintf(stderr, " Total variation diminishing set to %d %s.\n",CD->TVDFlg,tmpstr[1]);
    }
    if ( keymatch(line, "output",tmpval,tmpstr)==1){
      CD->OutItv = (int) tmpval[0];
      fprintf(stderr, " Output interval set to %d hours.\n",CD->OutItv);
    }
    if ( keymatch(line, "maximum_timestep", tmpval, tmpstr)==1){
      CD->TimMax = (double) tmpval[0];
      stepsize = CD->TimMax;
      fprintf(stderr, " Maximum time step is set to % 6.4f minutes.\n", CD->TimMax);
    }
    if ( keymatch(line, "minimum_timestep", tmpval, tmpstr)==1){
      CD->TimMin = (double) tmpval[0];
      start_step = CD->TimMin;
      fprintf(stderr, " Minimum time step is set to % 6.4g minutes.\n", CD->TimMin);
    }
    if ( keymatch(line, "overwrite_timestep", tmpval, tmpstr)==1){
      CD->TimLst = (double) tmpval[0];
      fprintf(stderr, " RT time step is set to % 6.4g minutes.\n", CD->TimLst);
    }
    if ( keymatch(line, "activity",tmpval,tmpstr)==1){
      CD->ACTmod = (int) tmpval[0];
      fprintf(stderr, " Activity correction is set to %d.\n",CD->ACTmod);
      // 0 for unity activity coefficient and 1 for DH equation update
    }
    if ( keymatch(line, "act_coe_delay",tmpval,tmpstr)==1){
      CD->DHEdel = (int) tmpval[0];
      fprintf(stderr, " Activity coefficient update delay is set to %d.\n",CD->DHEdel);
      // 0 for delay and 1 for no delay (solving together )
    }
    if ( keymatch(line, "thermo",tmpval,tmpstr)==1){
      CD->TEMcpl = (int) tmpval[0];
      fprintf(stderr, " Coupling of thermo modelling is set to %d.\n",CD->DHEdel);
      // 0 for delay and 1 for no delay (solving together )                            
    }
  }

  species Global_type;
  Global_type.ChemName = (char*) malloc(WORD_WIDTH* sizeof(char));
  strcpy(Global_type.ChemName,"GLOBAL");

  rewind(chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "GLOBAL", tmpval, tmpstr)!=1)
    fgets(line, line_width, chemfile);
  while (keymatch(line, "END", tmpval, tmpstr)!=1){
    fgets(line, line_width, chemfile);
    if( keymatch(line, "t_species", tmpval, tmpstr) == 1){
      CD->NumStc = (int) tmpval[0];
      fprintf(stderr, " %d chemical species specified.\n",CD->NumStc);
      /* H2O is always a primary species */
    }
    if( keymatch(line, "s_species", tmpval, tmpstr) == 1){
      CD->NumSsc = (int) tmpval[0];
      fprintf(stderr, " %d secondary species specified.\n", (int)tmpval[0]);
    }
    if( keymatch(line, "minerals", tmpval, tmpstr) == 1){
      CD->NumMin = (int) tmpval[0];
      fprintf(stderr, " %d minerals specified.\n", CD->NumMin);
    }
    if( keymatch(line, "adsorption", tmpval, tmpstr) == 1){
      CD->NumAds = (int) tmpval[0];
      fprintf(stderr, " %d surface complexation specified.\n", CD->NumAds);
    }
    if( keymatch(line, "cation_exchange", tmpval, tmpstr) == 1){
      CD->NumCex = (int) tmpval[0];
      fprintf(stderr, " %d cation exchange specified.\n", CD->NumCex);
    }
    if( keymatch(line, "mineral_kinetic", tmpval, tmpstr) == 1){
      CD->NumMkr = (int) tmpval[0];
      fprintf(stderr, " %d mineral kinetic reaction(s) specified.\n", CD->NumCex);
    }
    if( keymatch(line, "aqueous_kinetic", tmpval, tmpstr) == 1){
      CD->NumAkr = (int) tmpval[0];
      fprintf(stderr, " %d aqueous kinetic reaction(s) specified.\n", CD->NumCex);
    }
    if( keymatch(line, "diffusion", tmpval, tmpstr) == 1){
      fprintf(stderr, " Diffusion coefficient =%6.4f cm2/s.\n", tmpval[0]);
      Global_type.DiffCoe = tmpval[0] * 60.0 * 60.0 * 24.0 / 10000.0;
      Global_diff = 1;
      /* Require unit conversion ! */
    }
    if( keymatch(line, "dispersion", tmpval,tmpstr) == 1){
      fprintf(stderr, " Dispersion coefficient =%6.4f m.\n",tmpval[0]);
      Global_type.DispCoe = tmpval[0];
      Global_disp = 1;
      /* Set global flags to indicate the global values are present */
    }
    if( keymatch(line, "cementation", tmpval, tmpstr) == 1){
      fprintf(stderr, " Cementation factor set to %6.4f. \n", tmpval[0]);
      CD->Cementation = tmpval[0];
    }
    if ( keymatch(line, "temperature", tmpval, tmpstr) ==1){
      CD->Temperature = tmpval[0];
      fprintf(stderr, " Temperature set to %6.4f. \n", CD->Temperature);
    }
  }

  CD->NumSpc = CD->NumStc - (CD->NumMin + CD->NumAds + CD->NumCex); 
  /* the number of species that are mobile, later used in the OS3D subroutine */

  CD->NumSdc = CD->NumSpc;

  CD->Dependency = (double **) malloc (CD->NumSsc * sizeof(double*));
  for ( i = 0; i < CD->NumSsc; i ++)
    CD->Dependency[i] = (double*) malloc( CD->NumSdc * sizeof(double));
  /* convert secondary species as an expression of primary species */

  CD->Dep_kinetic= (double **) malloc ((CD->NumMkr + CD->NumAkr) * sizeof(double*));
  for ( i = 0; i < CD->NumMkr + CD->NumAkr; i ++)
    CD->Dep_kinetic[i] = (double*) malloc( CD->NumStc * sizeof(double));
  /* express kinetic species as function of primary species */

  CD->Dep_kinetic_all= (double **) malloc ((CD->NumMin) * sizeof(double*));
  for ( i = 0; i < CD->NumMin; i ++)
    CD->Dep_kinetic_all[i] = (double*) malloc( CD->NumStc * sizeof(double));
  /* Dependencies of minearls, all */


  CD->Keq = (double*) malloc( CD->NumSsc * sizeof(double));
  CD->KeqKinect = (double*) malloc( (CD->NumMkr + CD->NumAkr) * sizeof(double));
  CD->KeqKinect_all = (double * ) malloc( CD->NumMin * sizeof(double));
  /* Keqs of equilibrium/ kinetic and kinetic all */

  CD->Totalconc  = (double **) malloc( CD->NumStc * sizeof(double*));
  for (i = 0; i < CD->NumStc; i ++)
    CD->Totalconc[i] = (double*) malloc((CD->NumStc + CD->NumSsc)*sizeof(double));
  /* convert total concentration as an expression of all species */

  CD->Totalconck  = (double **) malloc( CD->NumStc * sizeof(double*));
  for (i = 0; i < CD->NumStc; i ++)
    CD->Totalconck[i] = (double*) malloc((CD->NumStc + CD->NumSsc)*sizeof(double));
  /* convert total concentration as an expression of all species */

  
  num_species= CD->NumSpc;

  rewind(chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "INITIAL_CONDITIONS", tmpval, tmpstr)!=1)
    fgets(line, line_width, chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "END", tmpval, tmpstr)!=1){
    if (keymatch(line, " ", tmpval, tmpstr)!=2) {
      num_conditions ++;
    }
    fgets(line, line_width, chemfile);
  }
  fprintf(stderr," %d conditions assigned.\n",num_conditions);
 


  char ** chemcon = (char**) malloc ( num_conditions * sizeof(char*));
  for ( i = 0; i < num_conditions; i ++)
    chemcon[i] = (char*) malloc( word_width * sizeof(char));
  char *** con_chem_name = (char***) malloc ( num_conditions *sizeof(char**));
  for ( i = 0; i < num_conditions; i ++){
    con_chem_name[i]= (char**) malloc ( CD->NumStc * sizeof(char*));
    for ( j = 0; j < CD->NumStc; j ++)
      con_chem_name[i][j] = (char*) malloc( WORD_WIDTH * sizeof(char));
  }





  int * condition_index = (int*) malloc( (CD->NumVol + 1) * sizeof(int)); 
  /* when user assign conditions to blocks, they start from 1 */

  for ( i = 0; i <= CD->NumVol; i ++)
    condition_index[i] = 0;

  vol_conc* Condition_vcele = (vol_conc*) malloc(num_conditions*sizeof(vol_conc));
  for ( i = 0 ; i< num_conditions; i++){
    Condition_vcele[i].index    = i+1;
    Condition_vcele[i].t_conc   = (double*) malloc(CD->NumStc *sizeof(double));
    Condition_vcele[i].p_conc   = (double*) malloc(CD->NumStc *sizeof(double));
    Condition_vcele[i].p_para   = (double*) malloc(CD->NumStc *sizeof(double));
    Condition_vcele[i].p_type   = (int *  ) malloc(CD->NumStc *sizeof(int   ));
    Condition_vcele[i].s_conc   = NULL;
    /* we do not input cocentration for secondary speices in rt */
    for (j = 0 ; j< CD->NumStc; j++){
      Condition_vcele[i].t_conc[j] = ZERO;
      Condition_vcele[i].p_conc[j] = ZERO;
    }
  }

  CD->chemtype = (species*) malloc((CD->NumStc + CD->NumSsc)*sizeof(species));
  for ( i = 0 ; i <CD->NumStc + CD->NumSsc; i++){
    if ( Global_diff == 1)
      CD->chemtype[i].DiffCoe = Global_type.DiffCoe;
    else
      CD->chemtype[i].DiffCoe = ZERO;
    /* in squre m per day */
    if ( Global_disp == 1)
      CD->chemtype[i].DispCoe = Global_type.DispCoe;
    else
      CD->chemtype[i].DispCoe = ZERO;
    CD->chemtype[i].ChemName = (char*) malloc(WORD_WIDTH * sizeof(char));
  }


 
  k = 0;
  rewind(chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "INITIAL_CONDITIONS", tmpval, tmpstr)!=1)
    fgets(line, line_width, chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "END", tmpval, tmpstr)!=1){
    if (keymatch(line, " ", tmpval, tmpstr)!=2) {
      strcpy(chemcon[k++], tmpstr[0]);
      fprintf(stderr, " Condition %s %d assigned to cells %s.\n",chemcon[k-1], k, tmpstr[1]);
      ConditionAssign(k, tmpstr[1], condition_index);
    }
    fgets(line, line_width, chemfile);
  }
  
  
  for ( i = 0; i < num_conditions; i++){
    rewind(chemfile);
    num_species = 0;
    num_mineral = 0;
    num_ads     = 0;
    num_cex     = 0;
    num_other   = 0;
    fgets(line, line_width, chemfile);
    while ((keymatch(line, "Condition", tmpval, tmpstr)!=1)||(keymatch(line,chemcon[i],tmpval,tmpstr) != 1))
      fgets(line, line_width, chemfile);
    if (strcmp(tmpstr[1], chemcon[i]) == 0)
      fprintf(stderr," %s", line);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr)!=1){
      if (keymatch(line, "NULL", tmpval, tmpstr)!=2) {
	if( SpeciationType(database, tmpstr[0]) == 1){
	  num_other = num_mineral+ num_ads + num_cex;
	  Condition_vcele[i].t_conc[num_species - num_other ] = tmpval[0];	
	  strcpy(con_chem_name[i][num_species - num_other], tmpstr[0]);
	  fprintf(stderr, " %s\t%6.4f\n", con_chem_name[i][num_species - num_other], tmpval[0]);
	  Condition_vcele[i].p_type[num_species - num_other ] = 1;
	}
	// arrange the concentration of the primary species in such a way that all the mobile species are at the beginning.
	if( SpeciationType(database, tmpstr[0]) == 4){
	  Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds + CD->NumCex + num_mineral] = tmpval[0];
	  if (strcmp(tmpstr[2], "-ssa")==0)
	    Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds + CD->NumCex + num_mineral] = tmpval[1];
	  strcpy(con_chem_name[i][CD->NumSpc + CD->NumAds + CD->NumCex + num_mineral], tmpstr[0]);
	  fprintf(stderr, " mineral %s\t%6.4f specific surface area %6.4f\n", con_chem_name[i][CD->NumSpc + CD->NumAds + CD->NumCex + num_mineral], tmpval[0], tmpval[1]);
	  Condition_vcele[i].p_type[CD->NumSpc + CD->NumAds + CD->NumCex + num_mineral] = 4;
	  num_mineral ++;
	}
	if( SpeciationType(database, tmpstr[0]) == 2){
	  Condition_vcele[i].t_conc[CD->NumSpc + num_ads] = tmpval[0];
	  Condition_vcele[i].p_type[CD->NumSpc + num_ads] = 2;
	  Condition_vcele[i].p_para[CD->NumSpc + num_ads] = 0;
	  // update when fill in the parameters for adsorption.
	  strcpy(con_chem_name[i][CD->NumSpc + num_ads], tmpstr[0]);
	  fprintf(stderr, " surface complex %s\t %6.4f\n", con_chem_name[i][CD->NumSpc + num_ads],tmpval[0]);
	  num_ads ++;
	  // under construction
	}
	if( SpeciationType(database, tmpstr[0]) == 3){
          Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds + num_cex] = tmpval[0];
	  Condition_vcele[i].p_type[CD->NumSpc + CD->NumAds + num_cex] = 3;
	  Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds + num_cex] = 0;
	  // update when fill in the parameters for cation exchange.
          strcpy(con_chem_name[i][CD->NumSpc + CD->NumAds + num_cex], tmpstr[0]);
          fprintf(stderr, " cation exchange %s\t %6.4f\n", con_chem_name[i][CD->NumSpc + CD->NumAds + num_cex],tmpval[0]);
	  num_cex ++;
	  // under construction
	}
	num_species ++;
      }
      fgets(line, line_width, chemfile);
    }
  }
  
  if (num_species != CD->NumStc) fprintf(stderr, " Number of species does not match indicated value!\n");

  for ( i = 1 ; i < num_conditions; i ++){
    for ( j = 0; j < num_species; j ++){
      if ( strcmp(con_chem_name[i][j], con_chem_name[i-1][j])!=0){
	error_flag = 1;
      }
    }
    if ( error_flag == 1)     fprintf(stderr, " The order of the chemicals in condition <%s> is incorrect!\n",chemcon[i]);
  }
  

  for ( i = 0; i < CD->NumStc; i ++){
    strcpy(CD->chemtype[i].ChemName, con_chem_name[0][i]);
    CD->chemtype[i].itype = Condition_vcele[0].p_type[i];
    fprintf(stderr, " %s\t%d\n",CD->chemtype[i].ChemName, CD->chemtype[i].itype);
  }
  
  
  rewind(chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "SECONDARY_SPECIES", tmpval, tmpstr)!=1)
    fgets(line, line_width, chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "END", tmpval, tmpstr)!=1){
    if (keymatch(line, "NULL", tmpval, tmpstr)!=2) {
      strcpy(CD->chemtype[num_species++].ChemName, tmpstr[0]);
      fprintf(stderr, " %s\n",CD->chemtype[num_species -1].ChemName);
    }
    fgets(line, line_width, chemfile);
  }
  int num_dep = 2;

  CD->kinetics = ( Kinetic_Reaction*) malloc( CD->NumMkr * sizeof( Kinetic_Reaction));
  for ( i = 0 ; i < CD->NumMkr ; i ++){
    CD->kinetics[i].species = (char * ) malloc( WORD_WIDTH * sizeof(char));
    CD->kinetics[i].Label   = (char * ) malloc( WORD_WIDTH * sizeof(char));
    CD->kinetics[i].dep_species = (char ** ) malloc( num_dep *sizeof(char*));
    CD->kinetics[i].dep_power = (double*) malloc( num_dep  * sizeof(double));
    CD->kinetics[i].monod   = (char **) malloc( num_dep    * sizeof(char*));
    CD->kinetics[i].monod_para= (double*) malloc( num_dep  * sizeof(double));
    CD->kinetics[i].inhibition=(char **) malloc( num_dep   * sizeof(char*));
    CD->kinetics[i].inhibition_para = (double*) malloc (num_dep * sizeof(double));
    CD->kinetics[i].dep_position = (int*) malloc( num_dep * sizeof(int));
    for ( j = 0 ; j < num_dep ; j ++){
      CD->kinetics[i].dep_species[j] = (char*) malloc(WORD_WIDTH*sizeof(char));
      CD->kinetics[i].monod[j]       = (char*) malloc(WORD_WIDTH*sizeof(char));
      CD->kinetics[i].inhibition[j]  = (char*) malloc(WORD_WIDTH*sizeof(char));
    }
  }
  k = 0;
  rewind(chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "MINERALS", tmpval, tmpstr)!=1)
    fgets(line, line_width, chemfile);
  fgets(line, line_width, chemfile);
  while (keymatch(line, "END", tmpval, tmpstr)!=1){
    if (keymatch(line, " ", tmpval, tmpstr)!=2) {
      strcpy(CD->kinetics[k].species, tmpstr[0]);
      if ( strcmp(tmpstr[1], "-label") == 0)
	strcpy(CD->kinetics[k].Label, tmpstr[2]);
      k ++;
    }
    fgets(line, line_width, chemfile);
  }
  for ( i = 0; i < k ; i ++)
    fprintf(stderr, " Kinetic reaction on %s is specified, label %s\n",CD->kinetics[i].species, CD->kinetics[i].Label);

  /* Define the boundary conditions here after obtaining the chemical condition labels */

  int boundarylayer[6];

  rewind(gridfile);
  fgets(line, line_width, gridfile);
  while ( keymatch(line, "BOUNDARY", tmpval, tmpstr) !=1)
    fgets(line, line_width, gridfile);
  fgets(line, line_width, gridfile);
  while ( keymatch(line, "END", tmpval, tmpstr) !=1){
    if (keymatch(line, "x_begin", tmpval, tmpstr) == 1){
      for ( i = 0; i < num_conditions; i++)
	if ( strcmp(tmpstr[1],chemcon[i]) == 0)
	  condition_index[num_blocks+1] = i + 1;
      if (strcmp(tmpstr[2],"flux")==0) boundarylayer[0] = 1;
      else boundarylayer[0] = 2;
      fprintf( stderr, " x_begin assigned to condition %s %d at cell %d with type %d\n", tmpstr[1], condition_index[num_blocks+1], num_blocks+1, boundarylayer[0]);
    }
    if (keymatch(line, "x_end", tmpval, tmpstr) == 1){
      for ( i = 0; i < num_conditions; i++)
	if ( strcmp(tmpstr[1],chemcon[i]) == 0)
	  condition_index[num_blocks+2] = i + 1;
      if (strcmp(tmpstr[2],"flux")==0) boundarylayer[0] = 1;
      else boundarylayer[0] = 2;
      fprintf( stderr, " x_end assigned to condition %s %d at cell %d with type %d\n", tmpstr[1], condition_index[num_blocks+2],  num_blocks+2, boundarylayer[0]);
    }
    if (keymatch(line, "y_begin", tmpval, tmpstr) == 1){
      for ( i = 0; i < num_conditions; i++)
	if ( strcmp(tmpstr[1],chemcon[i]) == 0)
	  condition_index[num_blocks+3] = i + 1;
      if (strcmp(tmpstr[2],"flux")==0) boundarylayer[0] = 1;
      else boundarylayer[0] = 2;
      fprintf( stderr, " y_begin assigned to condition %s at cell %d with type %d\n", tmpstr[1], num_blocks+3, boundarylayer[0]);
    }
    if (keymatch(line, "y_end", tmpval, tmpstr) == 1){
      for ( i = 0; i < num_conditions; i++)
	if ( strcmp(tmpstr[1],chemcon[i]) == 0)
	  condition_index[num_blocks+4] = i + 1;
      if (strcmp(tmpstr[2],"flux")==0) boundarylayer[0] = 1;
      else boundarylayer[0] = 2;
      fprintf( stderr, " y_end assigned to condition %s at cell %d with type %d\n", tmpstr[1], num_blocks+4, boundarylayer[0]);
    }
    if (keymatch(line, "z_begin", tmpval, tmpstr) == 1){
      for ( i = 0; i < num_conditions; i++)
	if ( strcmp(tmpstr[1],chemcon[i]) == 0)
	  condition_index[num_blocks+5] = i + 1;
      if (strcmp(tmpstr[2],"flux")==0) boundarylayer[0] = 1;
      else boundarylayer[0] = 2;
      fprintf( stderr, " z_begin assigned to condition %s at cell %d with type %d\n", tmpstr[1], num_blocks+5, boundarylayer[0]);
    }
    if (keymatch(line, "z_end", tmpval, tmpstr) == 1){
      for ( i = 0; i < num_conditions; i++)
	if ( strcmp(tmpstr[1],chemcon[i]) == 0)
	  condition_index[num_blocks+6] = i + 1;
      if (strcmp(tmpstr[2],"flux")==0) boundarylayer[0] = 1;
      else boundarylayer[0] = 2;
      fprintf( stderr, " z_end assigned to condition %s at cell %d with type %d\n", tmpstr[1], num_blocks+6, boundarylayer[0]);
    }
    fgets(line, line_width, gridfile);
  }
  /*
  for ( i = 0; i < CD->NumVol + 1; i ++){
    fprintf(stderr, " %d: %d\n", i, condition_index[i]);
  }
  */
  
  FILE* debug = fopen("logfile/flowdebug.log","w");


 
  CD->StartTime = starttime;
  CD->Vcele = (vol_conc*) malloc( (CD->NumVol) * sizeof(vol_conc));
  for ( i = 0 ; i < CD->NumSpc; i ++)
    if ( strcmp(CD->chemtype[i].ChemName, "pH") == 0){
      strcpy(CD->chemtype[i].ChemName, "H+");
      speciation_flg = 1;
    }
  
  /* Initializing concentration distributions */

  for ( i = 0 ; i< CD->NumVol; i++){
    CD->Vcele[i].index    = i+1;
    CD->Vcele[i].NumStc   = CD->NumStc;
    CD->Vcele[i].NumSsc   = CD->NumSsc;
    NewCell(&(CD->Vcele[i]), CD->NumStc, CD->NumSsc);
    CD->Vcele[i].q        = 0.0 ;
    for (j = 0 ; j< CD->NumStc; j++){
      if ((speciation_flg == 1) && (strcmp(CD->chemtype[j].ChemName, "H+")==0)){
        CD->Vcele[i].p_conc[j] = pow(10,-(Condition_vcele[condition_index[i+1]-1].t_conc[j]));
        CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
	CD->Vcele[i].t_conc[j] = CD->Vcele[i].p_conc[j];
        CD->Vcele[i].p_type[j] = 1;
      }
      else if (CD->chemtype[j].itype == 4){
        CD->Vcele[i].t_conc[j] = Condition_vcele[condition_index[i+1]-1].t_conc[j];
        CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
	CD->Vcele[i].p_actv[j] = 1.0;
        CD->Vcele[i].p_para[j] = Condition_vcele[condition_index[i+1]-1].p_para[j];
        CD->Vcele[i].p_type[j] = Condition_vcele[condition_index[i+1]-1].p_type[j];
      }
      else {
        CD->Vcele[i].t_conc[j] = Condition_vcele[condition_index[i+1]-1].t_conc[j];
	CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j]*0.5;
	CD->Vcele[i].p_para[j] = Condition_vcele[condition_index[i+1]-1].p_para[j];
        CD->Vcele[i].p_type[j] = Condition_vcele[condition_index[i+1]-1].p_type[j];
      }
    }
    for (j = 0 ; j< CD->NumSsc; j++){
      CD->Vcele[i].s_conc[j] = ZERO;
    }
  }
  for (i = 0; i < CD->NumVol; i ++){
    /*    if ( i < num_blocks * 0.5 ) {
      CD->Vcele[i].index    = i + 1;
      CD->Vcele[i].height_o = min(0.50 + 0.01 * i , 1.0);
      CD->Vcele[i].height_t = CD->Vcele[i].height_o;
      CD->Vcele[i].area     = x_length*y_length;
      CD->Vcele[i].porosity = porosity;
      CD->Vcele[i].sat      = CD->Vcele[i].height_t / z_length;
      CD->Vcele[i].sat_o    = CD->Vcele[i].sat;
      CD->Vcele[i].vol_o    = CD->Vcele[i].area * CD->Vcele[i].height_o;
      CD->Vcele[i].vol      = CD->Vcele[i].area * CD->Vcele[i].height_t;
    }
    else {*/
      CD->Vcele[i].index = i + 1;
      CD->Vcele[i].height_o = z_length;
      CD->Vcele[i].height_t = z_length;
      CD->Vcele[i].height_v = z_length;
      CD->Vcele[i].area     = x_length*y_length;
      CD->Vcele[i].porosity = porosity;
      CD->Vcele[i].sat      = 1.0;
      CD->Vcele[i].sat_o    = CD->Vcele[i].sat;
      CD->Vcele[i].vol_o    = CD->Vcele[i].area * CD->Vcele[i].height_o;
      CD->Vcele[i].vol      = CD->Vcele[i].area * CD->Vcele[i].height_t;
      // }
    // Room for future updates. Aim: allow heterogeneous porosity field/ conditional porosity input.

    CD->Vcele[i].q        = 0.0;
    CD->Vcele[i].BC       = 0;
  }
  for ( i = num_blocks; i < CD->NumVol; i++){
    CD->Vcele[i].BC       = 1;
    //boundary layer cells;
    //detailed location please refer to previous codes
  }
  /*
  for ( i = 0; i < num_blocks; i++)
    fprintf(debug, " Condition index! #%d# @ cell %d, sat: %f\n",condition_index[i+1],i, CD->Vcele[i].sat);
  */

  // arrange faces
  // assume that we have x, y, z as the number of blocks along x, y, z axis;
  // therefore we have 2xy + 2xz + 2yz boundary faces, connected to x/y/z_begin/end
  // we have 6xyz - 2(xy + xz + yz) interfaces between cells we are interested
  // total number of faces is 6xyz;
  CD->Flux  = (face*) malloc( 6*x_blocks*y_blocks*z_blocks *sizeof(face));
  int ii, jj, kk, centroid;

  /*  for (i = 1; i <= x_blocks; i ++)
    for ( j = 1; j <= y_blocks; j++)
      for (k = 1; k <= z_blocks; k++){
	// assign 6 connections per block
	fprintf(stderr, " index of (%d,%d,%d) is %d\n", i,j,k, centroid = Coordinates(i,j,k,x_blocks, y_blocks, z_blocks));
      }
  */

  /* 
   * Start of chemical speciation
   */

  Lookup(database, CD, 0);
  for ( i = 0; i < CD->NumMin; i++)
    fprintf(stderr, " KEQKINETIC = %6.4f\n ", CD->KeqKinect[i]);
  
  


    
  num_face = 0;
  for (i = 1; i <= x_blocks; i ++)
    for ( j = 1; j <= y_blocks; j++)
      for (k = 1; k <= z_blocks; k++){
	// assign 6 connections per block
	centroid = Coordinates(i,j,k,x_blocks, y_blocks, z_blocks);
	//	fprintf(stderr, " index of (%d,%d,%d) is %d\n", i,j,k, centroid);
	for ( ii = 0 ; ii < 6; ii ++){
	  CD->Flux[num_face+ii].nodeup = centroid + 1;
	  CD->Flux[num_face+ii].BC     = 0;
	}
	// Face 1 & 2 along x

	CD->Flux[num_face].nodelo   = Coordinates(i-1,j,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodell   = Coordinates(i-2,j,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodeuu   = Coordinates(i+1,j,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].distance = x_length;
	CD->Flux[num_face].velocity = -x_vel;
	CD->Flux[num_face].s_area   = y_length*z_length;
	CD->Flux[num_face].flux     = (-x_vel) * y_length*z_length;

	//	fprintf(stderr, " up: %d lo %d uu %d ll %d\n",CD->Flux[num_face].nodeup, CD->Flux[num_face].nodelo, CD->Flux[num_face].nodeuu, CD->Flux[num_face].nodell);


	num_face ++;
	CD->Flux[num_face].nodelo   = Coordinates(i+1,j,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodell   = Coordinates(i+2,j,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodeuu   = Coordinates(i-1,j,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].distance = x_length;
	CD->Flux[num_face].velocity = x_vel;
	CD->Flux[num_face].s_area   = y_length*z_length;
	CD->Flux[num_face].flux     = (x_vel) * y_length*z_length;

	num_face ++;
	CD->Flux[num_face].nodelo   = Coordinates(i,j-1,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodell   = Coordinates(i,j-2,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodeuu   = Coordinates(i,j+1,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].distance = y_length;
	CD->Flux[num_face].velocity = -y_vel;
	CD->Flux[num_face].s_area   = x_length*z_length;
	CD->Flux[num_face].flux     = (-y_vel) * x_length*z_length;

	num_face ++;
	CD->Flux[num_face].nodelo   = Coordinates(i,j+1,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodell   = Coordinates(i,j+2,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodeuu   = Coordinates(i,j-1,k,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].distance = y_length;
	CD->Flux[num_face].velocity = y_vel;
	CD->Flux[num_face].s_area   = x_length*z_length;
	CD->Flux[num_face].flux     = (y_vel) * x_length*z_length;

	num_face ++;
	CD->Flux[num_face].nodelo   = Coordinates(i,j,k-1,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodell   = Coordinates(i,j,k-2,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodeuu   = Coordinates(i,j,k+1,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].distance = z_length;
	CD->Flux[num_face].velocity = -z_vel;
	CD->Flux[num_face].s_area   = x_length*y_length;
	CD->Flux[num_face].flux     = (-z_vel) * x_length*y_length;

	num_face ++;
	CD->Flux[num_face].nodelo   = Coordinates(i,j,k+1,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodell   = Coordinates(i,j,k+2,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].nodeuu   = Coordinates(i,j,k-1,x_blocks,y_blocks,z_blocks)+1;
	CD->Flux[num_face].distance = z_length;
	CD->Flux[num_face].velocity = z_vel;
	CD->Flux[num_face].s_area   = x_length*y_length;
	CD->Flux[num_face].flux     = (z_vel) * x_length*y_length;

	num_face -=5;
	
	//	for ( ii = 0 ; ii < 6; ii ++){
		  //  fprintf(stderr, " Face: %d\t Up: %d Lo:%d LL: %d UU:%d\n",num_face+ii,CD->Flux[num_face+ii].nodeup, CD->Flux[num_face+ii].nodelo, CD->Flux[num_face+ii].nodell, CD->Flux[num_face+ii].nodeuu);
		  //  fprintf(stderr, " Distance: %6.4f Velocity: %6.4f S-area: %6.4f Flux: %6.4f\n", CD->Flux[num_face+ii].distance, CD->Flux[num_face+ii].velocity, CD->Flux[num_face+ii].s_area, CD->Flux[num_face+ii].flux);

	//	}
       
	num_face +=6;

      }
  int tmpcont=0;
  CD->NumFac = num_face;

  double tot_length = x_length * x_blocks;

  for ( i = 0 ; i < num_face; i ++)
    {
      if (CD->Flux[i].nodelo > num_blocks)
	{
	  CD->Flux[i].nodell = 0;
	  CD->Flux[i].nodeuu = 0;
	  CD->Flux[i].BC     = 1;
	  tmpcont ++;
	}
      fprintf(stderr, " Face: %d\t Up: %d Lo:%d LL: %d UU:%d\n", i ,CD->Flux[i].nodeup, CD->Flux[i].nodelo, CD->Flux[i].nodell, CD->Flux[i].nodeuu);
      fprintf(stderr, " Distance: %6.4f Velocity: %6.4f S-area: %6.4f Flux: %6.4f BC: %d\n ", CD->Flux[i].distance, CD->Flux[i].velocity, CD->Flux[i].s_area, CD->Flux[i].flux, CD->Flux[i].BC);

    }

  for ( i = 0 ; i < num_blocks; i ++){
    for ( j = 0 ; j < CD->NumStc; j ++){
      if (CD->chemtype[j].itype ==4){
        CD->Vcele[i].t_conc[j] = CD->Vcele[i].t_conc[j] * 1000 / CD->chemtype[j].MolarVolume / CD->Vcele[i].porosity / CD->Vcele[i].sat;
        CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
      }
    }

  }
  
  
  for ( i = 0; i < CD->NumVol; i++)
    Speciation(CD, i);
  CD->SPCFlg = 0;

  fprintf(stderr, " Number of boundary faces: %d\n", tmpcont);
  

  double * mass_t = (double*) malloc (CD->NumSpc * sizeof(double));
  double * mass_o = (double*) malloc (CD->NumSpc * sizeof(double));
  FILE  * outfile = fopen("logfile/testlog","w");
  FILE  * concfile = fopen("logfile/conclog.out","w");
  FILE  * actvfile = fopen("logfile/actvlog.out","w");
  FILE  * btcurve = fopen("logfile/btcurve.out","w");


  /*
  fprintf(outfile, "%4.2f\t",timelps);
  for ( j = 0; j < CD->NumStc; j ++)
    fprintf(outfile, " %s\t", CD->chemtype[j].ChemName);
  fprintf(outfile, "\n");
  
  for ( i = 0; i < CD->NumVol; i ++){
    fprintf(outfile, " %6.4f\t",(i+0.5)* x_length);
    for ( j = 0; j < CD->NumStc; j ++)
      fprintf(outfile, "%6.4f\t", log10(CD->Vcele[i].t_conc[j]));
    fprintf(outfile, "\n");
  }
  
  fprintf(outfile, "\n");
  */


  fprintf(concfile, "Total Species: %d\nPrimary Species:%d\nSecondary Species:%d\n",CD->NumStc, CD->NumSpc, CD->NumSsc);

  fprintf(concfile, "%4.2f\t",timelps);
  for ( j = 0; j < CD->NumStc + CD->NumSsc; j ++)
    fprintf(concfile, "%s\t", CD->chemtype[j].ChemName);
  fprintf(concfile, "\n");
  
  for ( i = 0; i < CD->NumVol; i ++){
    fprintf(concfile, " %6.4f\t",(i+0.5)* x_length);
    for ( j = 0; j < CD->NumStc; j ++)
      fprintf(concfile, "%6.4f\t", log10(CD->Vcele[i].p_conc[j]));
    for ( j = 0; j < CD->NumSsc; j ++)
      fprintf(concfile, "%6.4f\t", log10(CD->Vcele[i].s_conc[j]));
    fprintf(concfile, "\n");
  }
  
  fprintf(actvfile, "\n");

  fprintf(actvfile, "%4.2f\t",timelps);
  for ( j = 0; j < CD->NumStc + CD->NumSsc; j ++)
    fprintf(actvfile, " %s\t", CD->chemtype[j].ChemName);
  fprintf(actvfile, "\n");
  
  for ( i = 0; i < CD->NumVol; i ++){
    fprintf(actvfile, " %6.4f\t",(i+0.5)* x_length);
    for ( j = 0; j < CD->NumStc; j ++)
      fprintf(actvfile, "%6.4f\t", log10(CD->Vcele[i].p_actv[j]));
    for ( j = 0; j < CD->NumSsc; j ++)
      fprintf(actvfile, "%6.4f\t", log10(CD->Vcele[i].s_actv[j]));
    fprintf(actvfile, "\n");
  }
  
  fprintf(btcurve, "\n");

  for ( j = 0; j < CD->NumStc; j ++)
    fprintf(btcurve, " %s\t", CD->chemtype[j].ChemName);
  for ( j = 0; j < CD->NumSsc; j ++)
    fprintf(btcurve, " %s\t", CD->chemtype[CD->NumStc + j].ChemName);

  fprintf(btcurve, "\n");

  rt_step = CD->OutItv * 60;

  double t1, t2;
  t1 = timer();

  
  for ( i = 0 ; i < CD->NumFac; i++) CD->Vcele[CD->Flux[i].nodeup - 1].rt_step = 0.0;
  
  for ( i = 0; i < CD->NumFac; i++){
      for ( j = 0 ; j < CD->NumSpc; j++){
        peclet = fabs(CD->Flux[i].velocity * CD->Flux[i].distance / (CD->chemtype[j].DiffCoe + CD->chemtype[j].DispCoe * CD->Flux[i].velocity));
	peclet = MAX(peclet, 1.0E-10);
      }
      temp_rt_step = fabs(CD->Flux[i].flux / CD->Vcele[CD->Flux[i].nodeup - 1].vol/ CD->Vcele[CD->Flux[i].nodeup - 1].porosity) * ( 1 + peclet) / peclet;
      CD->Vcele[CD->Flux[i].nodeup - 1].rt_step += temp_rt_step;
      //      fprintf(stderr, " vel: %f dist: %f, diff: %f disp: %f pec:%f flux/porovolume: %f\n",CD->Flux[i].velocity, CD->Flux[i].distance, CD->chemtype[j].DiffCoe, CD->chemtype[j].DispCoe, peclet, temp_rt_step);
  }
  for ( i = 0 ; i < CD->NumOsv; i ++){
    //    fprintf(stderr, " Cell:%d, rtstep : %f\t", CD->Vcele[i].index, CD->Vcele[i].rt_step);
    CD->Vcele[i].rt_step = 0.8 * UNIT_C / CD->Vcele[i].rt_step;
    CD->Vcele[i].illness = 0;
    //    fprintf(stderr, " rtstep : %f\n", CD->Vcele[i].index, CD->Vcele[i].rt_step);
    if ( rt_step > CD->Vcele[i].rt_step)
      rt_step = CD->Vcele[i].rt_step;
    }

  rt_step     = MIN(CD->TimMax,rt_step);
  if ( CD->TimLst > 0) rt_step = CD->TimLst;

  while ( timelps <= endtime){
    
    fprintf(stderr, " ***************\n Writting data to file.\n ***************\n T: %f: Peclet Number: %6.4g\t Timestep constrained by Peclet: %6.4g\t Max timestep: %6.4g\n", timelps, peclet, rt_step, stepsize);

    // rt_step is the maximum time step Operator Splitting can take

    // fprintf(stderr, "\n Timestep Inherited  %6.4g\n", start_step );

    stepsize   = CD->TimMax;
    start_step = CD->TimMin;

    AdptTime(CD, &timelps, rt_step, stepsize, & start_step, num_blocks);
    //    masscheck(mass_o, CD, 0);
    //    OS3D(timelps, stepsize*0.5, CD);
    //  for ( i = 0 ; i < num_blocks; i++)
    //  ReactControl(timelps, stepsize, CD, i);
    //   OS3D(timelps, stepsize*0.5, CD);
    // timelps += stepsize;
    

    //    masscheck(mass_t, CD, 1);
    if ( ((int)timelps % ( 60 )) == 0){
      fprintf(btcurve, "%6.4f\t", timelps);
      
      for ( j = 0; j < CD->NumStc; j ++){
	fprintf(btcurve, "%6.4g\t", (CD->Vcele[num_blocks-1].p_conc[j]));
      }
      for ( j = 0; j < CD->NumSsc; j ++)
	fprintf(btcurve, "%6.4g\t", (CD->Vcele[num_blocks-1].s_conc[j]));
      
      fprintf(btcurve, "\n");
    }

    if ((int)timelps == 72000){

      //      fprintf(outfile, "%4.2f\n",timelps);
	/*for ( j = 0; j < CD->NumStc; j ++)
	fprintf(outfile, " %s", CD->chemtype[j].ChemName);
	fprintf(outfile, "\n");*/
      for ( i = 0; i < num_blocks; i ++){
	fprintf(outfile, " %6.4f\t",(i+0.5)* x_length);
	for ( j = 0; j < CD->NumStc; j ++)
	  fprintf(outfile, "%6.4f\t", log10(CD->Vcele[i].t_conc[j]));
	fprintf(outfile, "\n");
      }   
      fprintf(outfile, "\n");
      
    }

    if ( ((int)timelps % (CD->OutItv * 60)) == 0){
      fprintf(concfile, "%4.2f\n",timelps);
	/*for ( j = 0; j < CD->NumStc; j ++)
	fprintf(concfile, " %s", CD->chemtype[j].ChemName);
	fprintf(concfile, "\n");*/
      for ( i = 0; i < num_blocks; i ++){
	fprintf(concfile, " %6.4f\t",(i+0.5)* x_length);
	for ( j = 0; j < CD->NumStc; j ++){
	  if (CD->chemtype[j].itype == 4)
	    fprintf(concfile, "%12.10f\t ", CD->Vcele[i].p_conc[j]);
	  else
	    fprintf(concfile, "%6.4f\t", log10(CD->Vcele[i].p_conc[j]));
	}
	for ( j = 0; j < CD->NumSsc; j ++)
	  fprintf(concfile, "%6.4f\t", log10(CD->Vcele[i].s_conc[j]));

	fprintf(concfile, "\n");
      }   
      fprintf(concfile, "\n");

      fprintf(actvfile, "%4.2f\n",timelps);
	/*for ( j = 0; j < CD->NumStc; j ++)
	fprintf(actvfile, " %s", CD->chemtype[j].ChemName);
	fprintf(actvfile, "\n");*/
      for ( i = 0; i < num_blocks; i ++){
	fprintf(actvfile, " %6.4f\t",(i+0.5)* x_length);
	for ( j = 0; j < CD->NumStc; j ++)
	  fprintf(actvfile, "%6.4f\t", log10(CD->Vcele[i].p_actv[j]));
	for ( j = 0; j < CD->NumSsc; j ++)
	  fprintf(actvfile, "%6.4f\t", log10(CD->Vcele[i].s_actv[j]));

	fprintf(actvfile, "\n");
      }   
      fprintf(actvfile, "\n");



    }
    //  fprintf(stderr, " mass: %12.8f\t%12.8f\t%12.8f\t%12.8f\n",mass_o[0],mass_t[0],mass_o[0]-mass_t[0], (mass_o[0]-mass_t[0])/mass_o[0]*100);
  }
  t2 = timer();
  fprintf(stderr, " \n Simulation done in %6.4g seconds. Please check datafiles for results. Thanks for using RT!\n", t2-t1);


  fclose(outfile);
  fclose(debug);
  fclose(concfile);
  fclose(actvfile);
  fclose(chemfile);
  fclose(gridfile);
  fclose(database);
  fclose(btcurve);
  return (0);
  
}
void AdptTimeBak(Chem_Data * CD, realtype * timelps, double rt_step, double total_step, double * start_step, int num_blocks)
{
  double stepsize = * start_step, org_time = * timelps, step_rst = * start_step, t = * timelps;
  int i,j, k, nr_tmp, nr_max, tot_nr;

  /*
  while ( * timelps < org_time + CD->OutItv * 60){
    nr_max = 5;
    if (stepsize > org_time + CD->OutItv * 60 - * (timelps)){
      // before adjusting, write the current timestep to file
      step_rst = stepsize;
      stepsize = org_time + CD->OutItv * 60 - * (timelps);
      fprintf(stderr, " Adjusting Timestep to %6.4f to match output time series", stepsize);

    }
    OS3D(*timelps, stepsize*0.5, CD);
    for ( i = 0 ; i < num_blocks; i++){
      //      Speciation(CD,i);
      React(*(timelps), stepsize, CD, i, &nr_tmp);
      if (nr_tmp > nr_max)  nr_max = nr_tmp;
    }
    OS3D(*timelps, stepsize*0.5, CD);
    //   Speciation(CD,i);
    *(timelps) += stepsize;
    if ( *(timelps) >= org_time + CD->OutItv * 60) break;

    if ( nr_max < 10) stepsize *= 2;
    if ( nr_max < 15) stepsize *= 1.5;
    if ( nr_max > 20) stepsize *= 0.5;
    if (stepsize > min(rt_step, hydro_step)) stepsize = min(rt_step,hydro_step);

  }
  if ( step_rst >= CD->OutItv * 60) step_rst = CD->OutItv * 60;
  * start_step = step_rst;
  */
  stepsize = MIN( rt_step, * start_step);

  while ( t < org_time + total_step){
    if (stepsize > org_time + total_step - t ){
      // before adjusting, write the current timestep to file
      step_rst = stepsize;
      stepsize = org_time + total_step - t ;
      fprintf(stderr, " Adjusting Timestep to %6.4f to match output time series\n", stepsize);
    }
    OS3D(t , stepsize, CD);
    //    for ( i = 0 ; i < num_blocks; i++){
      //      Speciation(CD,i);
    // React(t, stepsize, CD, i, &nr_tmp);
      // if (nr_tmp > nr_max)  nr_max = nr_tmp;
    // }
    //   Speciation(CD,i);
    t += stepsize;
    if ( t >= org_time + total_step) break;

    //    if ( nr_max < 10) stepsize *= 2;
    //    if ( nr_max < 15) stepsize *= 1.5;
    //    if ( nr_max > 20) stepsize *= 0.5;
    //    if (stepsize > min(rt_step, total_step)) stepsize = min(rt_step,total_step);
  }
  if ( step_rst >= total_step ) step_rst = total_step;
  * start_step = step_rst;
  * timelps    = t;

  double substep = 0.0;
  
  for ( i = 0 ; i < num_blocks; i++){
    if (React(org_time , total_step, CD, i, &nr_tmp))
      {
	fprintf(stderr, "  ---> React failed at cell %12d.\t" , CD->Vcele[i].index);
	substep = 0.5 * total_step;
	k       = 2;
	while ( j =  React(org_time, substep, CD, i , &nr_tmp))
	  {
	    substep = 0.5 * substep;
	    k       = 2   * k;
	    if ( substep < 0.5) break;
	  }
	if ( j == 0)
	  {
	    tot_nr  += nr_tmp;
	    fprintf(stderr, " Reaction passed with step equals to %f (1/%d)\n", substep,k);
	    for ( j = 1; j < k ; j ++)
	      {
		React(org_time + j * substep, substep, CD, i, &nr_tmp);
		tot_nr += nr_tmp;
	      }
	  }
	else
	  {
	    fprintf(stderr, " Reaction failed anyway\n");
	  }
      }
  }
}


void AdptTime(Chem_Data * CD, realtype * timelps, double rt_step, double hydro_step, double * start_step, int num_blocks)
{
  double stepsize = * start_step, org_time = * timelps, step_rst = * start_step;
  int i,j, nr_tmp, nr_max, nr_tot;
  while ( * timelps < org_time + CD->OutItv * 60){
    nr_max = 5;
    if (stepsize > org_time + CD->OutItv * 60 - * (timelps)){
      // before adjusting, write the current timestep to file
      step_rst = stepsize;
      stepsize = org_time + CD->OutItv * 60 - * (timelps);
      fprintf(stderr, " Adjusting Timestep to %6.4f to match output time series\n", stepsize);

    }
    nr_tot = 0;
    OS3D(*timelps, stepsize*0.5, CD);
    for ( i = 0 ; i < num_blocks; i++){
      //      Speciation(CD,i);
      React(*(timelps), stepsize, CD, i, &nr_tmp);
      if (nr_tmp > nr_max)  nr_max = nr_tmp;
      nr_tot += nr_tmp;
    }
    OS3D(*timelps, stepsize*0.5, CD);
    //   Speciation(CD,i);
    *(timelps) += stepsize;
    if ( *(timelps) >= org_time + CD->OutItv * 60) break;

    fprintf(stderr, " Average NR iteration taken %f\n",(double)nr_tot/num_blocks);

    if ( nr_max < 10) stepsize *= 5;
    if ( nr_max < 15) stepsize *= 2;
    if ( nr_max > 20) stepsize *= 0.5;
    if (stepsize > min(rt_step, hydro_step)) stepsize = min(rt_step,hydro_step);

  }
  if ( step_rst >= CD->OutItv * 60) step_rst = CD->OutItv * 60;
  * start_step = step_rst;
}
