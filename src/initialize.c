/********************************************************************************
 * File		: initialize.c							*
 * Function	: Initialization of elemental attributes using relational	*
 *		  database							*
 * File		: perturb.c							*
 * Function	: Perturb parameters for Flux-PIHM-EnKF				*
 * Version	: Sep, 2011 (2.0)						*
 * Developer of Flux-PIHM:	Yuning Shi (yshi@psu.edu)			*                            
 * Developer of PIHM2.0:	Mukesh Kumar (muk139@psu.edu)		        * 
 * Developer of PIHM1.0:	Yizhong Qu   (quyizhong@gmail.com)	        * 
 *------------------------------------------------------------------------------*
 *										*
 * For questions or comments, please contact					*
 *      --> Yuning Shi (yshi@psu.edu)						*
 *      --> Mukesh Kumar (muk139@psu.edu)					*
 *      --> Prof. Chris Duffy (cxd11@psu.edu)					*
 * This code is free for research purpose only.					*
 * Please provide relevant references if you use this code in your research work*
 *------------------------------------------------------------------------------*
*********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sundialstypes.h"
#include "nvector_serial.h"
#include "pihm.h"  

#define UNIT_C 1440     /* 60*24 for calculation of yDot in m/min units while forcing is in m/day. */

realtype Interpolation(TSD *Data, realtype t);

realtype FieldCapacity(realtype Alpha, realtype Beta, realtype Kv, realtype ThetaS, realtype ThetaR)
	{
	realtype elemSatn, Ktemp;
	realtype ThetaRef;

	for (elemSatn = 0.005; elemSatn<1; elemSatn = elemSatn+0.001)
		{
		Ktemp = Kv*(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,Beta/(Beta-1)),(Beta-1)/Beta),2);
		if (Ktemp>=0.0005)
			{
			ThetaRef = (1./3.+2./3.*elemSatn)*(ThetaS-ThetaR);
			break;
			}
		}
	return ThetaRef;
	}

void initialize(char *filename, Model_Data DS, Control_Data *CS, N_Vector CV_Y)
	{
  	int i,j,k,tmpBool,BoolBR,BoolR=0;
  	realtype a_x, a_y, b_x, b_y, c_x, c_y,distX,distY;
  	realtype a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax; 
  	realtype tempvalue1,tempvalue2,tempvalue3;
	realtype lapserate,dsoilsum=0;
	realtype AquiferDepth, h, dz[5], hsw;
	realtype dmac, droot;
	int macroporelayer, nroot;
	realtype elemSatn, Ktemp;
	int nlayer, bot;
  	FILE *init_file;
  	char *fn;
  	realtype *zmin_cor;
	realtype Tb = 9.0;
	realtype zb = 3.0;

  	zmin_cor=(realtype *)malloc(DS->NumEle*sizeof(realtype));
  
  	printf("\nInitializing data structure ... ");

	DS->dsoil = (realtype **)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */

  	/* allocate memory storage to flux terms */
  	DS->FluxSurf = (realtype **)malloc(DS->NumEle*sizeof(realtype*));
  	DS->FluxSub = (realtype **)malloc(DS->NumEle*sizeof(realtype*));
  	DS->FluxRiv = (realtype **)malloc(DS->NumRiv*sizeof(realtype*));
	DS->VeloSub = (realtype **)malloc(DS->NumEle*sizeof(realtype*));
	DS->DistSub = (realtype **)malloc(DS->NumEle*sizeof(realtype*));
	DS->AreaSub = (realtype **)malloc(DS->NumEle*sizeof(realtype*));
  	DS->EleET = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleH = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
  	DS->EleG = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
  	DS->EleLE = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->EleOVLbuffer = (realtype *)malloc(DS->NumEle*sizeof(realtype));
	DS->EleGWbuffer = (realtype *)malloc(DS->NumEle*sizeof(realtype));
	DS->Eleunsatbuffer =(realtype *)malloc(DS->NumEle*sizeof(realtype)); 
  	DS->EleEp = (realtype *)malloc((DS->NumEle+DS->NumRiv)*sizeof(realtype));			/* Expanded by Y. Shi */
  	DS->EleSM = (realtype **)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
  	DS->EleDew = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->C_h = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->C_m = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->Tsoil = (realtype **)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->Tsfc = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->Tbot = (realtype *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
	DS->EleBot = (int *)malloc(DS->NumEle*sizeof(realtype));			/* Expanded by Y. Shi */
  	DS->ElePrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleViR = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->Recharge = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleIS = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleISmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleISsnowmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleSnow = (realtype *)malloc((DS->NumEle+DS->NumRiv)*sizeof(realtype));  
  	DS->EleSnowGrnd = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleSnowCanopy = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleTF = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleETloss = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleNetPrep = (realtype *)malloc((DS->NumEle+DS->NumRiv)*sizeof(realtype));

	DS->TF = 1.0*CS->Cal.TF;
	DS->IS = 1.0*CS->Cal.IS;
	DS->Czil = 0.1*CS->Cal.Czil;
	DS->fx_soil = 1.0*CS->Cal.fx_soil;
	DS->fx_canopy = 0.5*CS->Cal.fx_canopy;
	DS->Tref = 24.85*CS->Cal.Tref;

 	for (i=0; i<DS->NumSoil; i++)
		{
		DS->Soil[i].ThetaW = 0.5*(DS->Geol[i].ThetaS-DS->Geol[i].ThetaR)*pow(1/(1+pow(200*DS->Geol[i].Alpha,DS->Geol[i].Beta)),1-1/(DS->Geol[i].Beta));
//		DS->Soil[i].ThetaW = 0.25*(DS->Geol[i].ThetaS-DS->Geol[i].ThetaR)*pow(1/(1+pow(200*DS->Soil[i].Alpha,DS->Soil[i].Beta)),1-1/(DS->Soil[i].Beta))+DS->Geol[i].ThetaR;
		DS->Soil[i].ThetaRef = FieldCapacity(DS->Geol[i].Alpha,DS->Geol[i].Beta,DS->Geol[i].KsatV,DS->Geol[i].ThetaS,DS->Geol[i].ThetaR);
//		DS->Soil[i].ThetaRef = FieldCapacity(DS->Soil[i].Alpha,DS->Soil[i].Beta,DS->Geol[i].KsatV,DS->Geol[i].ThetaS,DS->Geol[i].ThetaR)+DS->Geol[i].ThetaR;
//		DS->Soil[i].ThetaRef = CS->Cal.ThetaRef*0.75*(DS->Geol[i].ThetaS-DS->Geol[i].ThetaR)+DS->Geol[i].ThetaR;

//		printf("\n\n Ks = %f",DS->Geol[i].KsatV);
//		printf("\n Beta = %f",DS->Soil[i].Beta);
//		printf("\n ThetaS = %f",DS->Geol[i].ThetaS);
//		printf("\n ThetaR = %f",DS->Geol[i].ThetaR);
		printf("\n ThetaRef = %f, %f",DS->Soil[i].ThetaRef,FieldCapacity(DS->Soil[i].Alpha,DS->Soil[i].Beta,DS->Geol[i].KsatV,DS->Geol[i].ThetaS,DS->Geol[i].ThetaR)/(DS->Geol[i].ThetaS-DS->Geol[i].ThetaR));
		printf("\n ThetaW = %f, %f\n",DS->Soil[i].ThetaW,0.5*pow(1/(1+pow(200*DS->Soil[i].Alpha,DS->Soil[i].Beta)),1-1/(DS->Soil[i].Beta)));
		}

	for(i=0;i<DS->NumEle;i++)
		{
    		DS->FluxSurf[i] = (realtype *)malloc(3*sizeof(realtype));
    		DS->FluxSub[i] = (realtype *)malloc(3*sizeof(realtype));
		DS->VeloSub[i] = (realtype *)malloc(3*sizeof(realtype));
		DS->DistSub[i] = (realtype *)malloc(3*sizeof(realtype));              /* Modified by C. Bao*/
		DS->AreaSub[i] = (realtype *)malloc(3*sizeof(realtype));
    		DS->EleET[i] = (realtype *)malloc(4*sizeof(realtype));
	  	DS->EleSM[i] = (realtype *)malloc(5*sizeof(realtype));			/* Expanded by Y. Shi */
		DS->Tsoil[i] = (realtype *)malloc(4*sizeof(realtype));			/* Expanded by Y. Shi */
		DS->C_h[i] = 0.0001*UNIT_C*60;			/* Expanded by Y. Shi */
		DS->C_m[i] = 0.0001*UNIT_C*60;			/* Expanded by Y. Shi */
		DS->EleBot[i] = 4;

		for ( j = 0 ; j < 4 ; j ++)
		  DS->EleET[i][j] = 0.0;
	       	for ( j = 0 ; j < 5 ; j ++)
		  DS->EleSM[i][j] = 0.0;

		DS->Recharge[i]    = 0.0;
		DS->EleSnowGrnd[i] = 0.0;
	
    		a_x = DS->Node[DS->Ele[i].node[0]-1].x;
    		b_x = DS->Node[DS->Ele[i].node[1]-1].x;
    		c_x = DS->Node[DS->Ele[i].node[2]-1].x;
    		a_y = DS->Node[DS->Ele[i].node[0]-1].y;
    		b_y = DS->Node[DS->Ele[i].node[1]-1].y;
    		c_y = DS->Node[DS->Ele[i].node[2]-1].y;

    		a_zmin = DS->Node[DS->Ele[i].node[0]-1].zmin;
    		b_zmin = DS->Node[DS->Ele[i].node[1]-1].zmin;
    		c_zmin = DS->Node[DS->Ele[i].node[2]-1].zmin;
    		a_zmax = DS->Node[DS->Ele[i].node[0]-1].zmax;
    		b_zmax = DS->Node[DS->Ele[i].node[1]-1].zmax;
    		c_zmax = DS->Node[DS->Ele[i].node[2]-1].zmax;
    
    		DS->Ele[i].area = 0.5*((b_x - a_x)*(c_y - a_y) - (b_y - a_y)*(c_x - a_x));
    		DS->Ele[i].zmax = (a_zmax + b_zmax + c_zmax)/3.0;
    		DS->Ele[i].zmin = (a_zmin + b_zmin + c_zmin)/3.0; 
    		DS->Ele[i].edge[0] = pow((b_x - c_x), 2) + pow((b_y - c_y), 2);
    		DS->Ele[i].edge[1] = pow((c_x - a_x), 2) + pow((c_y - a_y), 2);
    		DS->Ele[i].edge[2] = pow((a_x - b_x), 2) + pow((a_y - b_y), 2);

		DS->dsoil[i] = (realtype *)malloc(5*sizeof(realtype));
		AquiferDepth=DS->Ele[i].zmax-DS->Ele[i].zmin;
		if (AquiferDepth<=0.1)
			{
			DS->dsoil[i][0] = AquiferDepth;
			for (j=1;j<5;j++)
				{
				DS->dsoil[i][j] = -999.0;
				}
			}
		else if (AquiferDepth<=0.4)
			{
			DS->dsoil[i][0] = 0.1;
			DS->dsoil[i][1] = AquiferDepth - 0.1;
			for (j=2;j<5;j++)
				{
				DS->dsoil[i][j] = -999.0;
				}
			}
		else if (AquiferDepth<=1.0)
			{
			DS->dsoil[i][0] = 0.1;
			DS->dsoil[i][1] = 0.3;
			DS->dsoil[i][2] = AquiferDepth - 0.4;
			for (j=3;j<5;j++)
				{
				DS->dsoil[i][j] = -999.0;
				}
			}
		else if (AquiferDepth<=2.0)
			{
			DS->dsoil[i][0] = 0.1;
			DS->dsoil[i][1] = 0.3;
			DS->dsoil[i][2] = 0.6;
			DS->dsoil[i][3] = AquiferDepth - 1.0;
			DS->dsoil[i][4] = -999.0;
			}
		else
			{
			DS->dsoil[i][0] = 0.1;
			DS->dsoil[i][1] = 0.3;
			DS->dsoil[i][2] = 0.6;
			DS->dsoil[i][3] = 1.0;
			DS->dsoil[i][4] = AquiferDepth - 2.0;
			}

    		/* calculate centroid of triangle */
    		DS->Ele[i].x = (a_x + b_x + c_x)/3.0;
    		DS->Ele[i].y = (a_y + b_y + c_y)/3.0;
    
    
    		/* calculate circumcenter of triangle */
  		/*  DS->Ele[i].x = a_x - ((b_y - a_y)*DS->Ele[i].edge[2] - (c_y - a_y)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
    		DS->Ele[i].y = a_y + ((b_x - a_x)*DS->Ele[i].edge[2] - (c_x - a_x)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area); 
    		*/
    		DS->Ele[i].edge[0] = sqrt(DS->Ele[i].edge[0]);
    		DS->Ele[i].edge[1] = sqrt(DS->Ele[i].edge[1]);
    		DS->Ele[i].edge[2] = sqrt(DS->Ele[i].edge[2]);
    		DS->Ele[i].KsatH = CS->Cal.KsatH*DS->Geol[(DS->Ele[i].geol-1)].KsatH;
    		DS->Ele[i].KsatV = CS->Cal.KsatV*DS->Geol[(DS->Ele[i].geol-1)].KsatV;
    		DS->Ele[i].infKsatV = CS->Cal.infKsatV*DS->Soil[(DS->Ele[i].soil-1)].KsatV;
//    		DS->Ele[i].Porosity = CS->Cal.Porosity*(DS->Soil[(DS->Ele[i].soil-1)].ThetaS - DS->Soil[(DS->Ele[i].soil-1)].ThetaR);
		/* Note above porosity statement should be replaced by geologic porosity (in comments below) if the data is available */
		DS->Ele[i].Porosity = CS->Cal.Porosity*(DS->Geol[(DS->Ele[i].geol-1)].ThetaS - DS->Geol[(DS->Ele[i].geol-1)].ThetaR);

		DS->Ele[i].ThetaS = DS->Geol[(DS->Ele[i].geol-1)].ThetaS;
		DS->Ele[i].ThetaR = DS->Geol[(DS->Ele[i].geol-1)].ThetaR;

		DS->Ele[i].ThetaW = CS->Cal.ThetaW*DS->Soil[(DS->Ele[i].soil-1)].ThetaW+DS->Ele[i].ThetaR;
		DS->Ele[i].ThetaRef = CS->Cal.ThetaRef*DS->Soil[(DS->Ele[i].soil-1)].ThetaRef+DS->Ele[i].ThetaR;

 		if((DS->Ele[i].Porosity>1)&&(DS->Ele[i].Porosity==0))
			{
			printf("Warning: Porosity value out of bounds");
			getchar();
			} 
//		DS->Ele[i].Alpha = CS->Cal.Alpha*DS->Soil[(DS->Ele[i].soil-1)].Alpha;
//    		DS->Ele[i].Beta = CS->Cal.Beta*DS->Soil[(DS->Ele[i].soil-1)].Beta; 
		/* Note above van genuchten statement should be replaced by geologic parameters (in comments below) if the data is available */
		DS->Ele[i].Alpha = CS->Cal.Alpha*DS->Geol[(DS->Ele[i].geol-1)].Alpha;
    		DS->Ele[i].Beta = CS->Cal.Beta*DS->Geol[(DS->Ele[i].geol-1)].Beta; 
    		DS->Ele[i].hAreaF = CS->Cal.hAreaF*DS->Soil[(DS->Ele[i].soil-1)].hAreaF; 
    		DS->Ele[i].vAreaF = CS->Cal.vAreaF*DS->Geol[(DS->Ele[i].geol-1)].vAreaF; 
    		DS->Ele[i].macKsatV = CS->Cal.macKsatV*DS->Soil[(DS->Ele[i].soil-1)].macKsatV; 
    		DS->Ele[i].macKsatH = CS->Cal.macKsatH*DS->Geol[(DS->Ele[i].geol-1)].macKsatH; 
		DS->Ele[i].macD=CS->Cal.macD*DS->Geol[DS->Ele[i].geol-1].macD;
		if (DS->Ele[i].macD > DS->Ele[i].zmax - DS->Ele[i].zmin) DS->Ele[i].macD = DS->Ele[i].zmax - DS->Ele[i].zmin;
    		DS->Ele[i].infD=CS->Cal.infD*DS->Soil[DS->Ele[i].soil-1].infD;
    
    		DS->Ele[i].RzD = CS->Cal.RzD*DS->LandC[DS->Ele[i].LC-1].RzD;  
    		DS->Ele[i].LAImax = DS->LandC[DS->Ele[i].LC-1].LAImax;
    		DS->Ele[i].Rmin = CS->Cal.Rmin*DS->LandC[DS->Ele[i].LC-1].Rmin;
    		DS->Ele[i].Rs_ref = CS->Cal.Rs_ref*DS->LandC[DS->Ele[i].LC-1].Rs_ref;

		if (DS->Ele[i].Macropore==1)
			{	
			dmac = 0;
			j = 0;
			macroporelayer = 0;
			while(dmac<DS->Ele[i].macD && j<5)
				{
				dmac = dmac+DS->dsoil[i][j];
				macroporelayer = j;
				j++;
				}
			DS->Ele[i].MPL = macroporelayer;
			}


		// Looking for the deepest layer that root reaches
		j = 0;
		nroot = 0;
		droot = 0;

		while(droot<DS->Ele[i].RzD)
			{
			droot = droot+DS->dsoil[i][j];
			nroot = j;
			j++;
			}
		DS->Ele[i].rootL = nroot;

//		printf("D0 = %f, D1 = %f, D2 = %f, D3 = %f, D4 = %f, macropore layer = %d, root layer = %d\n",DS->dsoil[i][0],DS->dsoil[i][1],DS->dsoil[i][2],DS->dsoil[i][3],DS->dsoil[i][4],DS->Ele[i].MPL,DS->Ele[i].rootL);

    		DS->Ele[i].Albedo_min = CS->Cal.Albedo*DS->LandC[DS->Ele[i].LC-1].Albedo_min;
    		DS->Ele[i].Albedo_max = CS->Cal.Albedo*DS->LandC[DS->Ele[i].LC-1].Albedo_max;

    		DS->Ele[i].Emiss_min = DS->LandC[DS->Ele[i].LC-1].Emiss_min;
    		DS->Ele[i].Emiss_max = DS->LandC[DS->Ele[i].LC-1].Emiss_max;

    		DS->Ele[i].z0_min = DS->LandC[DS->Ele[i].LC-1].z0_min;
    		DS->Ele[i].z0_max = DS->LandC[DS->Ele[i].LC-1].z0_max;

//		if(DS->Ele[i].Albedo>1)
//			{
//                      printf("Warning: Albedo out of bounds");
//                      getchar();
//			}

    		DS->Ele[i].VegFrac = CS->Cal.VegFrac*DS->LandC[DS->Ele[i].LC-1].VegFrac;                
    		DS->Ele[i].Rough = CS->Cal.Rough*DS->LandC[DS->Ele[i].LC-1].Rough;
    		DS->Ele[i].h_s = CS->Cal.h_s*DS->LandC[DS->Ele[i].LC-1].h_s;
		DS->Ele[i].windH = DS->windH[DS->Ele[i].WindVel-1];



//		printf("\n Rmin = %f, Rs_ref = %f, h_s = %f", DS->Ele[i].Rmin*24*3600,DS->Ele[i].Rs_ref/24/3600,DS->Ele[i].h_s);

		}
//		printf("\n RzD = %f, %f, %f, %f, %f, %f, %f, %f",DS->Ele[i-1].RzD, DS->Ele[i-1].LAImax, DS->Ele[i-1].Rmin, DS->Ele[i-1].Rs_ref, DS->Ele[i-1].Albedo, DS->Ele[i-1].VegFrac, DS->Ele[i-1].Rough, DS->Ele[i-1].windH);

  	for(i=0; i<DS->NumRiv; i++)
  		{
    		DS->FluxRiv[i] = (realtype *)malloc(11*sizeof(realtype));
		for(j=0;j<3;j++)
			{
			/* Note: Strategy to use BC < -4 for river identification */
			if(DS->Ele[DS->Riv[i].LeftEle-1].nabr[j]==DS->Riv[i].RightEle)
				{
				DS->Ele[DS->Riv[i].LeftEle-1].BC[j]=-4*(i+1);
				}
			if(DS->Ele[DS->Riv[i].RightEle-1].nabr[j]==DS->Riv[i].LeftEle)
				{
				DS->Ele[DS->Riv[i].RightEle-1].BC[j]=-4*(i+1);
				}
			}
    		DS->Riv[i].x = (DS->Node[DS->Riv[i].FromNode-1].x + DS->Node[DS->Riv[i].ToNode-1].x)/2;
    		DS->Riv[i].y = (DS->Node[DS->Riv[i].FromNode-1].y + DS->Node[DS->Riv[i].ToNode-1].y)/2;
    		DS->Riv[i].zmax = (DS->Node[DS->Riv[i].FromNode-1].zmax + DS->Node[DS->Riv[i].ToNode-1].zmax)/2;
    		DS->Riv[i].depth = CS->Cal.rivDepth*DS->Riv_Shape[DS->Riv[i].shape-1].depth;
		DS->Riv[i].coeff=CS->Cal.rivShapeCoeff*DS->Riv_Shape[DS->Riv[i].shape - 1].coeff;
    		DS->Riv[i].zmin = DS->Riv[i].zmax - DS->Riv[i].depth;  
    		DS->Riv[i].Length = sqrt(pow(DS->Node[DS->Riv[i].FromNode-1].x -DS->Node[DS->Riv[i].ToNode-1].x, 2) + pow(DS->Node[DS->Riv[i].FromNode-1].y - DS->Node[DS->Riv[i].ToNode-1].y, 2));
		DS->Riv[i].KsatH=CS->Cal.rivKsatH*DS->Riv_Mat[DS->Riv[i].material-1].KsatH;
		DS->Riv[i].KsatV=CS->Cal.rivKsatV*DS->Riv_Mat[DS->Riv[i].material-1].KsatV;
		DS->Riv[i].bedThick=CS->Cal.rivbedThick*DS->Riv_Mat[DS->Riv[i].material-1].bedThick;
		DS->Riv[i].Rough=CS->Cal.rivRough*DS->Riv_Mat[DS->Riv[i].material - 1].Rough;
		/* Initialization for rectangular cells beneath river */
		/* Note: Ideally this data should be read from the decomposition itself */
		/* but it is not supported right now in PIHMgis (Bhatt, G and Kumar, M; 2007) */
		DS->Ele[i+DS->NumEle].zmax=DS->Riv[i].zmin;
		DS->Ele[i+DS->NumEle].zmin=DS->Riv[i].zmax-(0.5*(DS->Ele[DS->Riv[i].LeftEle-1].zmax+DS->Ele[DS->Riv[i].RightEle-1].zmax)-0.5*(DS->Ele[DS->Riv[i].LeftEle-1].zmin+DS->Ele[DS->Riv[i].RightEle-1].zmin));
//		DS->Ele[i+DS->NumEle].zmin=DS->Riv[i].zmax-40;
		DS->Ele[i+DS->NumEle].macD=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macD+DS->Ele[DS->Riv[i].RightEle-1].macD)>DS->Riv[i].depth?0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macD+DS->Ele[DS->Riv[i].RightEle-1].macD)-DS->Riv[i].depth:0;	
		DS->Ele[i+DS->NumEle].macKsatH=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macKsatH+DS->Ele[DS->Riv[i].RightEle-1].macKsatH);
		DS->Ele[i+DS->NumEle].vAreaF=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].vAreaF+DS->Ele[DS->Riv[i].RightEle-1].vAreaF);
		DS->Ele[i+DS->NumEle].KsatH=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].KsatH+DS->Ele[DS->Riv[i].RightEle-1].KsatH);
		DS->Ele[i+DS->NumEle].Porosity=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].Porosity+DS->Ele[DS->Riv[i].RightEle-1].Porosity);
 		}


	for(i=0;i<DS->NumPrep;i++)
		{
		DS->TSD_Prep[i].iCounter = 0;		/* Modified by Y. Shi */

		for(j=0; j<DS->TSD_Prep[i].length; j++)
			{
			DS->TSD_Prep[i].TS[j][1]=CS->Cal.Prep*DS->TSD_Prep[i].TS[j][1];
			}
		}
        for(i=0;i<DS->NumTemp;i++)
                {
		DS->TSD_Temp[i].iCounter = 0;		/* Modified by Y. Shi */

                for(j=0; j<DS->TSD_Temp[i].length; j++)
                        {
                        DS->TSD_Temp[i].TS[j][1]=CS->Cal.Temp*DS->TSD_Temp[i].TS[j][1];
                        }
                }

        for(i=0;i<DS->NumHumidity;i++)
                {
		DS->TSD_Humidity[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumWindVel;i++)
                {
		DS->TSD_WindVel[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumSdown;i++)
                {
		DS->TSD_Sdown[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumLdown;i++)
                {
		DS->TSD_Ldown[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumP;i++)
                {
		DS->TSD_Pressure[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumLC;i++)
                {
		DS->TSD_LAI[i].iCounter = 0;		/* Modified by Y. Shi */
		DS->TSD_RL[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumMeltF;i++)
                {
		DS->TSD_MeltF[i].iCounter = 0;		/* Modified by Y. Shi */
                }
        for(i=0;i<DS->NumSource;i++)
                {
		DS->TSD_Source[i].iCounter = 0;		/* Modified by Y. Shi */
                }


	/* Memory allocation of print variables */	
        for(i=0;i<35;i++)				/* Modified by Y. Shi */
                {
                if(i==0)
                        {
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumEle+DS->NumRiv,sizeof(realtype));
                        }
		else if (i == 6)
			{
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumEle+DS->NumRiv,sizeof(realtype));
                        }
                else if((i>=7)&&(i<19))
                        {
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumRiv,sizeof(realtype));
                        }
                else
                        {
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumEle,sizeof(realtype));
                        }
                }
	/* Debugging artifacts in data created due to coarser resolution of model elements */
	if(CS->Debug==1)
		{
		for(i=0; i<DS->NumEle; i++)
			{
			/* Correction of Surf Elev (artifacts due to coarse scale discretization). Not needed if there is lake feature.*/
			tmpBool=1;
			for(j=0;j<3;j++)
				{
				if(DS->Ele[i].nabr[j]>0)
					{
					tempvalue1=DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax;
					if(DS->Ele[i].zmax-tempvalue1>=0)
						{
						tmpBool=0;
						break;
						}
					}
				}	
			if(tmpBool==1)
				{
				printf("\n Ele %d is sink ",i+1);
				/* Note: Following correction is being applied for debug==1 case only */
				printf("\tBefore: %lf Corrected using:",DS->Ele[i].zmax); 
				tempvalue1=10000000;
				for(j=0;j<3;j++)
					{
					if(DS->Ele[i].nabr[j]>0)
						{
						DS->Ele[i].zmax=(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax);
						tempvalue1=tempvalue1>DS->Ele[i].zmax?DS->Ele[i].zmax:tempvalue1;
						printf("(%d)%lf  ",j+1,(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmax:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax));
						}
					}
				DS->Ele[i].zmax=tempvalue1;
				printf("=(New)%lf  ",DS->Ele[i].zmax);
				} 
			}
		/* Correction of BedRck Elev. Is this needed? */ 
		/*
		printf("\n Do you want to correct Bed Rock Elev too (1[y]/0[n])");
		scanf("%d",&BoolBR);
		*/
		BoolBR=1;
		if(BoolBR==1)
		{
		for(i=0; i<DS->NumEle; i++)
			{
			tmpBool=1;
			for(j=0;j<3;j++)
				{
				if(DS->Ele[i].nabr[j]>0)
					{
					tempvalue1=DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmin:DS->Ele[-(DS->Ele[i].BC[j]/4)-1+DS->NumEle].zmin;
					if(DS->Ele[i].zmin-tempvalue1>=0)
						{
						tmpBool=0;
						break;
						}
					}
				}	
			if(tmpBool==1)
				{
				printf("\n Ele %d is sink ",i+1);
				/* Note: Following correction is being applied for debug==1 case only */
				printf("\tBfore: %lf Corrected using:",DS->Ele[i].zmin); 
				tempvalue1=10000000;
				for(j=0;j<3;j++)
					{
					if(DS->Ele[i].nabr[j]>0)
						{
						DS->Ele[i].zmin=(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmin:DS->Ele[-(DS->Ele[i].BC[j]/4)-1+DS->NumEle].zmin);
						tempvalue1=tempvalue1>DS->Ele[i].zmin?DS->Ele[i].zmin:tempvalue1;
						printf("(%d)%lf  ",j+1,(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].zmin:DS->Ele[-(DS->Ele[i].BC[j]/4)-1+DS->NumEle].zmin));
						}
					}
				DS->Ele[i].zmin=tempvalue1;
				printf("=(New)%lf  ",DS->Ele[i].zmin);
				} 
			}
		}
//		getchar();
		printf("\nHit any key to see more details");
		for(i=0;i<DS->NumRiv;i++)
			{
			if(DS->Riv[i].down>0)
				{
				if(DS->Riv[i].zmin<DS->Riv[DS->Riv[i].down-1].zmin)
					{
					BoolR=1;
					printf("\n Riv %d is lower than downstream Riv %d",i+1,DS->Riv[i].down);
					}
				}
			}
		if(BoolR==1)
			{
			printf("\n\tRiver elevation correction needed");
			getchar();
			}	
		}
  	for(i=0; i<DS->NumEle; i++)
  		{
                a_x = DS->Node[DS->Ele[i].node[0]-1].x;
                b_x = DS->Node[DS->Ele[i].node[1]-1].x;
                c_x = DS->Node[DS->Ele[i].node[2]-1].x;
                a_y = DS->Node[DS->Ele[i].node[0]-1].y;
                b_y = DS->Node[DS->Ele[i].node[1]-1].y;
                c_y = DS->Node[DS->Ele[i].node[2]-1].y;


		for(j=0;j<3;j++)
			{
			/* Note: Assumption here is that the forumulation is circumcenter based */
			switch(j)
				{
				case 0:
                			distX=(DS->Ele[i].x-0.5*(b_x+c_x));
                			distY=(DS->Ele[i].y-0.5*(b_y+c_y));
					break;
                                case 1:
                			distX=(DS->Ele[i].x-0.5*(c_x+a_x));
                			distY=(DS->Ele[i].y-0.5*(c_y+a_y));
					break;
                                case 2:
                			distX=(DS->Ele[i].x-0.5*(a_x+b_x));
                			distY=(DS->Ele[i].y-0.5*(a_y+b_y));
					break;
				}
			DS->Ele[i].surfH[j]=(DS->Ele[i].nabr[j]>0)?(DS->Ele[i].BC[j]>-4?(DS->Ele[DS->Ele[i].nabr[j]-1].zmax):DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax):DS->Ele[i].BC[j]<=-4?DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax:(DS->Ele[i].zmax); 
			DS->Ele[i].surfX[j]=(DS->Ele[i].nabr[j]>0)?(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].x:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].x):(DS->Ele[i].x-2*distX);          
			DS->Ele[i].surfY[j]=DS->Ele[i].nabr[j]>0?(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].y:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].y):(DS->Ele[i].y-2*distY);
                	}
		DS->Ele[i].dhBYdx=-(DS->Ele[i].surfY[2]*(DS->Ele[i].surfH[1]-DS->Ele[i].surfH[0])+DS->Ele[i].surfY[1]*(DS->Ele[i].surfH[0]-DS->Ele[i].surfH[2])+DS->Ele[i].surfY[0]*(DS->Ele[i].surfH[2]-DS->Ele[i].surfH[1]))/(DS->Ele[i].surfX[2]*(DS->Ele[i].surfY[1]-DS->Ele[i].surfY[0])+DS->Ele[i].surfX[1]*(DS->Ele[i].surfY[0]-DS->Ele[i].surfY[2])+DS->Ele[i].surfX[0]*(DS->Ele[i].surfY[2]-DS->Ele[i].surfY[1]));  
		DS->Ele[i].dhBYdy=-(DS->Ele[i].surfX[2]*(DS->Ele[i].surfH[1]-DS->Ele[i].surfH[0])+DS->Ele[i].surfX[1]*(DS->Ele[i].surfH[0]-DS->Ele[i].surfH[2])+DS->Ele[i].surfX[0]*(DS->Ele[i].surfH[2]-DS->Ele[i].surfH[1]))/(DS->Ele[i].surfY[2]*(DS->Ele[i].surfX[1]-DS->Ele[i].surfX[0])+DS->Ele[i].surfY[1]*(DS->Ele[i].surfX[0]-DS->Ele[i].surfX[2])+DS->Ele[i].surfY[0]*(DS->Ele[i].surfX[2]-DS->Ele[i].surfX[1]));
  		}
  	/* initialize state variable */
  	/* relax case */
	if (CS->init_type == 0)
  		{
    		for(i=0; i<DS->NumEle; i++)
    			{
      			DS->EleIS[i] = 0;
			DS->EleSnow[i]=0;
			/* Note Two components can be separately read too */
			DS->EleSnowGrnd[i]=0;
			DS->EleSnowCanopy[i]=0;
      			NV_Ith_S(CV_Y, i) = 0;
      			NV_Ith_S(CV_Y, i + DS->NumEle) = 0.1;
     			NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele[i].zmax - DS->Ele[i].zmin -0.1;
//      			NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele[i].zmax - DS->Ele[i].zmin -0.5;
			DS->DummyY[i] = 0;
			DS->DummyY[i+DS->NumEle] = 0.1;
			DS->DummyY[i+2*DS->NumEle] = DS->Ele[i].zmax - DS->Ele[i].zmin -0.1;
			DS->EleOVLbuffer[i] = 0;
			DS->EleGWbuffer[i] = DS->Ele[i].zmax - DS->Ele[i].zmin -0.1;
			DS->Eleunsatbuffer[i] = 0.1; 
//			DS->DummyY[i+2*DS->NumEle] = DS->Ele[i].zmax - DS->Ele[i].zmin -0.5;

//			DS->Tsfc[i]=DS->TSD_Temp[0].TS[0][1];
			DS->Tsfc[i]= Interpolation(&DS->TSD_Temp[DS->Ele[i].temp-1], CS->StartTime);
			DS->Tbot[i]=Tb;		
			lapserate = (DS->Tsfc[i]-DS->Tbot[i])/zb;
			DS->Tsoil[i][0]=DS->Tsfc[i]-lapserate*DS->dsoil[i][0]/2;

			for(j=1; j<4; j++)
				{
				DS->Tsoil[i][j]=DS->Tsoil[i][j-1]-lapserate*(DS->dsoil[i][j-1]+DS->dsoil[i][j])/2;
	      			}

			for(j=0;j<5;j++)							/* Expanded by Y. Shi */
				{
				DS->EleSM[i][j] = DS->Geol[(DS->Ele[i].geol-1)].ThetaS;
				}
			}

    		for(i=0; i<DS->NumRiv; i++)
    			{
			DS->EleSnow[i+DS->NumEle] = 0;
      			NV_Ith_S(CV_Y, i + 3*DS->NumEle) = 0;
			/* Note once the element beneath river is incorporated in decomposition and .mesh file, initialization should be perfomed based on the location data instead of average of neighbor properties */
			NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv)=(DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin) -0.1;
			DS->DummyY[i+3*DS->NumEle] = 0.1;
			DS->DummyY[i+3*DS->NumEle+DS->NumRiv]=(DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin) -0.1;
    			}
  		}
  	/* data initialization mode */  
  	else if (CS->init_type == 1)
  		{

		for(i=0; i<DS->NumEle; i++)
			{
			DS->Tsfc[i]= Interpolation(&DS->TSD_Temp[DS->Ele[i].temp-1], CS->StartTime);
			DS->Tbot[i]=Tb;		
			lapserate = (DS->Tsfc[i]-DS->Tbot[i])/zb;
			DS->Tsoil[i][0]=DS->Tsfc[i]-lapserate*DS->dsoil[i][0]/2;

			for(j=1; j<4; j++)
				{
				DS->Tsoil[i][j]=DS->Tsoil[i][j-1]-lapserate*(DS->dsoil[i][j-1]+DS->dsoil[i][j])/2;
	      			}
			}

		if(DS->UnsatMode ==1)
			{    
    			}
		if(DS->UnsatMode ==2)
			{    
    			for(i=0; i<DS->NumEle; i++)
    				{
      				DS->EleIS[i] = DS->Ele_IC[i].interception;
				DS->EleSnow[i]=DS->Ele_IC[i].snow;
				/* Note Two components can be separately read too */
      				NV_Ith_S(CV_Y, i) = DS->Ele_IC[i].surf;
				/* Note: delete 0.1 here */
      				NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[i].unsat;
				NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele_IC[i].sat; 
				/* Note: delete line below for general*/
//     				NV_Ith_S(CV_Y, i + 2*DS->NumEle) = 0*DS->Ele_IC[i].sat+(DS->Ele[i].zmax - DS->Ele[i].zmin)*0.1;
      				if ((NV_Ith_S(CV_Y, i + DS->NumEle) + NV_Ith_S(CV_Y, i + 2*DS->NumEle)) >= (DS->Ele[i].zmax - DS->Ele[i].zmin))
      					{
        				NV_Ith_S(CV_Y, i + DS->NumEle) = ((DS->Ele[i].zmax - DS->Ele[i].zmin) - NV_Ith_S(CV_Y, i + 2*DS->NumEle))*0.98;
        				if (NV_Ith_S(CV_Y, i + DS->NumEle) < 0) 
						{
						NV_Ith_S(CV_Y, i + DS->NumEle) = 0; 
						}
      					} 
    				}	  
    			for(i=0; i<DS->NumRiv; i++)
    				{
				DS->EleSnow[i+DS->NumEle] = 0;
      				NV_Ith_S(CV_Y, i + 3*DS->NumEle) = DS->Riv_IC[DS->Riv[i].IC-1].value;
				/* Note once the element beneath river is incorporated in decomposition and .mesh file, initialization should be perfomed based on the location data instead of average of neighbor properties */
//				NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv) = 0.5*(DS->Ele_IC[DS->Riv[i].LeftEle-1].sat+DS->Ele_IC[DS->Riv[i].RightEle-1].sat);
				NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv) = (DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin)-0.1;
    				}
    			}    	
  		}  
  	/* hot start mode */
  	else
  		{
    		fn = (char *)malloc((strlen(filename)+12)*sizeof(char));
    		strcpy(fn, "input/");
		strcat(fn,filename);
    		init_file = fopen(strcat(fn, ".init"), "r");
  
    		if(init_file == NULL)
    			{
      			printf("\n  Fatal Error: %s.init is in use or does not exist!\n", filename);
      			exit(1);
    			}
    		else
    			{
      			for(i=0; i<DS->NumEle; i++)
      				{
        			fscanf(init_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &DS->EleIS[i],&DS->EleSnow[i],&tempvalue1,&tempvalue2,&tempvalue3,&DS->Tsfc[i],&DS->Tsoil[i][0],&DS->Tsoil[i][1],&DS->Tsoil[i][2],&DS->Tsoil[i][3],&DS->EleSM[i][0],&DS->EleSM[i][1],&DS->EleSM[i][2],&DS->EleSM[i][3],&DS->EleSM[i][4]);
				DS->Tbot[i]=Tb;
				NV_Ith_S(CV_Y, i)=tempvalue1;
				NV_Ith_S(CV_Y, i + DS->NumEle)=tempvalue2;
				NV_Ith_S(CV_Y, i + 2*DS->NumEle)=tempvalue3;
				DS->DummyY[i] = tempvalue1;
				DS->DummyY[i+DS->NumEle] = tempvalue2;
				DS->DummyY[i+2*DS->NumEle] = tempvalue3;
				DS->EleOVLbuffer[i] = tempvalue1;
				DS->Eleunsatbuffer[i] = tempvalue2;
				DS->EleGWbuffer[i] = tempvalue3;
				AquiferDepth=DS->Ele[i].zmax-DS->Ele[i].zmin;
				h = tempvalue3>0?tempvalue3:0;

				for (j=0; j<4;j++)
					{
					dz[j] = DS->dsoil[i][j];
					}

				if (AquiferDepth > 2.0+0.1)	
					{
					nlayer = 5;
					dz[4] = AquiferDepth - 2.0;
					}
				else
					{
					nlayer = 4;
					}

				if (AquiferDepth - h<dz[0]+dz[1]/2)
					{
					bot = 0;
					}
				else
					{
					// Looking for the deepest layer above water table

					hsw = dz[0]+dz[1];
					k = 2;

					while(h+hsw+dz[k]/2<AquiferDepth & k < nlayer)
						{
						hsw = hsw+dz[k];
						k++;
						}
					bot = k-1;
					bot = bot>=nlayer?(nlayer-1):bot;
					}
				DS->EleBot[i] = bot;
				}
      			for(i=0; i<DS->NumRiv; i++)
      				{
        			fscanf(init_file, "%lf %lf %lf",&tempvalue1,&tempvalue2,&tempvalue3);
        			NV_Ith_S(CV_Y, i + 3*DS->NumEle) = tempvalue1;
				NV_Ith_S(CV_Y, i + 3*DS->NumEle+DS->NumRiv) =tempvalue2;
				DS->DummyY[i+3*DS->NumEle] = tempvalue1;
				DS->DummyY[i+3*DS->NumEle+DS->NumRiv]=tempvalue2;
				DS->EleSnow[i+DS->NumEle] = tempvalue3;
      				} 
    			}
		free(fn);
    		fclose(init_file); 
  		}
	free(zmin_cor);
  	printf("\n done.\n");
	}
