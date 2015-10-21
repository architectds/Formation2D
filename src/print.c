/*******************************************************************************
 * File        : print.c	                                               *
 * Version     : Nov, 2007 (2.0)                                               *
 * Function    : print out model results output files                          *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) This file is downgraded from Version 1.0, as no ancillary results are    *
 *    being output			                                       *
 * b) Only state variables and flux to/in/accross river and its bed are being  *
 *    output							               *
 * c) Addition of Average Function to output average variables at regular time *
 *    intervals								       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"  
#include "cvode.h" 
#include "cvdense.h"

#define UNIT_C 1440     /* Note 60*24 for calculation of yDot in m/min units while forcing is in m/day. */
#define GRAV 9.8*60*60  /* Note the dependence on physical units */ 

realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool);

/*Temporal average of State vectors */
void avgResults_NV(FILE *fpin,realtype *tmpVarCal,N_Vector tmpNV,int tmpIntv, int tmpNumObj,realtype tmpt,int tmpInitObj)
        {
        int j;
	struct tm *timestamp;
	time_t *rawtime;

	rawtime = (time_t *)malloc(sizeof(time_t));

        for(j=0;j<tmpNumObj;j++)
                {
                tmpVarCal[j]=tmpVarCal[j]+NV_Ith_S(tmpNV,j+tmpInitObj);
                }
        if(((int)tmpt%tmpIntv)==0)      
                {
		*rawtime = (int)tmpt*60;
		timestamp = gmtime(rawtime);
                fprintf(fpin,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
                for(j=0;j<tmpNumObj;j++)
                        {               
                        fprintf(fpin,"%lf\t",tmpVarCal[j]/tmpIntv);
                        tmpVarCal[j]=0; 
                        }
                fprintf(fpin,"\n");     
                fflush(fpin);           
                }
	free(rawtime);
        }

void avgResults_NV_GW(FILE *fpin,realtype *tmpVarCal,N_Vector tmpNV,int tmpIntv, int tmpNumObj_1, int tmpNumObj_2,realtype tmpt,int tmpInitObj)
        {
        int j;
	struct tm *timestamp;		/* Modified by Y. Shi */
	time_t *rawtime;		/* Modified by Y. Shi */

	rawtime = (time_t *)malloc(sizeof(time_t));		/* Modified by Y. Shi */

        for(j=0;j<tmpNumObj_1;j++)
                {
                tmpVarCal[j]=tmpVarCal[j]+NV_Ith_S(tmpNV,j+tmpInitObj);
                }
        for(j=tmpNumObj_1;j<tmpNumObj_1+tmpNumObj_2;j++)
                {
                tmpVarCal[j]=tmpVarCal[j]+NV_Ith_S(tmpNV,j+tmpInitObj+tmpNumObj_2);
                }
        if(((int)tmpt%tmpIntv)==0)      
                {
//		fprintf(fpin,"%lf\t",tmpt);
		*rawtime = (int)tmpt*60;		/* Modified by Y. Shi */
		timestamp = gmtime(rawtime);		/* Modified by Y. Shi */
                fprintf(fpin,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);		/* Modified by Y. Shi */
                for(j=0;j<tmpNumObj_1+tmpNumObj_2;j++)
                        {               
                        fprintf(fpin,"\t%lf",tmpVarCal[j]/tmpIntv);
                        tmpVarCal[j]=0; 
                        }
                fprintf(fpin,"\n");     
                fflush(fpin);           
                }
	free(rawtime);
        }

/* Temporal average of Derived states */
void avgResults_MD(FILE *fpin,realtype *tmpVarCal,Model_Data tmpDS,int tmpIntv, int tmpNumObj,realtype tmpt,int tmpFC, N_Vector CV_Y)
        {
        int j;

	struct tm *timestamp;
	time_t *rawtime;

	rawtime = (time_t *)malloc(sizeof(time_t));
        switch(tmpFC)
                {
                case 3: 
                case 4:
                case 5:                 
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleET[j][tmpFC-3];
                                }
                        break;
                case 6: 
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleIS[j];
                                }
                        break;
                case 7:
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleSnow[j];		/* Edited by Y. Shi */
                                }
                        break;
                case 8:
//			for(j=0;j<tmpNumObj;j++)
 //                               {
//				tmpVarCal[j]=tmpVarCal[j]+tmpDS->FluxRiv[j][tmpFC-8];
//                                }
                case 9:
//			for(j=0;j<tmpNumObj;j++)
//                                {
//				tmpVarCal[j]=tmpVarCal[j]+(NV_Ith_S(CV_Y,j + 3*tmpDS->NumEle)<=0)?0:CS_AreaOrPerem(tmpDS->Riv_Shape[tmpDS->Riv[j].shape - 1].interpOrd,NV_Ith_S(CV_Y,j + 3*tmpDS->NumEle),tmpDS->Riv[j].coeff,1)*sqrt(GRAV*UNIT_C*UNIT_C*NV_Ith_S(CV_Y,j + 3*tmpDS->NumEle));
//                                }			
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->FluxRiv[j][tmpFC-8];
                                }
                        break;
                case 19:
                        for(j=0;j<tmpNumObj;j++)
                                {
                                tmpVarCal[j]=tmpVarCal[j]+tmpDS->Recharge[j];
//				if (j==1) printf("\n print.c Recharge = %f", tmpDS->Recharge[j]);
                                }
                        break;
		/* Modified by Y. Shi */
		case 20:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleG[j];
				}
			break;
		case 21:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleH[j];
				}
			break;
		case 22:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleLE[j];
				}
			break;
		case 23:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->Tsfc[j];
				}
			break;
		case 24:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->Tsoil[j][0];
				}
			break;

		case 25:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->Tsoil[j][1];
				}
			break;
		case 26:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->Tsoil[j][2];
				}
			break;
		case 27:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->Tsoil[j][3];
				}
			break;
		case 28:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleSM[j][0];
				}
			break;
		case 29:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleSM[j][1];
				}
			break;
		case 30:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleSM[j][2];
				}
			break;
		case 31:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleSM[j][3];
				}
			break;

		case 32:
			for(j=0;j<tmpNumObj;j++)
				{
				tmpVarCal[j]=tmpVarCal[j]+tmpDS->EleSM[j][4];
				}
			break;

                default:
                        break;
                }
        if(((int)tmpt%tmpIntv)==0)
                {
		*rawtime = (int)tmpt*60;
		timestamp = gmtime(rawtime);
                fprintf(fpin,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
//                fprintf(fpin,"%lf\t",tmpt);
                for(j=0;j<tmpNumObj;j++)
                        {
                        fprintf(fpin,"\t%lf",tmpVarCal[j]/tmpIntv);
                        tmpVarCal[j]=0;
                        }
                fprintf(fpin,"\n");
                fflush(fpin);
                }
	free(rawtime);
        }
/* print individual states */
void PrintData(FILE **outp,Control_Data *cD, Model_Data DS, N_Vector CV_Y, realtype t)
	{
	int k;
	if(cD->gwD==1)
		{
        	avgResults_NV_GW(outp[0],DS->PrintVar[0],CV_Y,cD->gwDInt,DS->NumEle, DS->NumRiv,t,2*DS->NumEle);
		}
        if(cD->surfD==1)
                {
        	avgResults_NV(outp[1],DS->PrintVar[1],CV_Y,cD->surfDInt,DS->NumEle,t,0*DS->NumEle);
		}
	for(k=0;k<3;k++)
		{
		if(cD->et[k]==1)
			{
			avgResults_MD(outp[2+k],DS->PrintVar[2+k],DS,cD->etInt,DS->NumEle,t,3+k,CV_Y);
			}
		}
        if(cD->IsD==1)
                {
        	avgResults_MD(outp[5],DS->PrintVar[5],DS,cD->IsDInt,DS->NumEle,t,6,CV_Y);
                }
        if(cD->snowD==1)
                {
        	avgResults_MD(outp[6],DS->PrintVar[6],DS,cD->snowDInt,DS->NumEle+DS->NumRiv,t,7,CV_Y);
                }
        for(k=0;k<=10;k++)
        	{
		if(cD->rivFlx[k]==1)
			{
        		avgResults_MD(outp[7+k],DS->PrintVar[k+7],DS,cD->rivFlxInt,DS->NumRiv,t,k+8,CV_Y);
			}
        	}
        if(cD->rivStg==1)
                {
        	avgResults_NV(outp[18],DS->PrintVar[18],CV_Y,cD->rivStgInt,DS->NumRiv,t,3*DS->NumEle);
                }
        if(cD->Rech==1)
                {
       		avgResults_MD(outp[20],DS->PrintVar[20],DS,cD->RechInt,DS->NumEle,t,19,CV_Y);
                }
        if(cD->usD==1)
                {
        	avgResults_NV(outp[19],DS->PrintVar[19],CV_Y,cD->usDInt,DS->NumEle,t,1*DS->NumEle);
                }
	if(cD->lsv==1)
		{
		avgResults_MD(outp[21],DS->PrintVar[21],DS,cD->etInt,DS->NumEle,t,20,CV_Y);
		avgResults_MD(outp[22],DS->PrintVar[22],DS,cD->etInt,DS->NumEle,t,21,CV_Y);
		avgResults_MD(outp[23],DS->PrintVar[23],DS,cD->etInt,DS->NumEle,t,22,CV_Y);
		avgResults_MD(outp[24],DS->PrintVar[24],DS,cD->etInt,DS->NumEle,t,23,CV_Y);
		avgResults_MD(outp[25],DS->PrintVar[25],DS,cD->etInt,DS->NumEle,t,24,CV_Y);
		avgResults_MD(outp[26],DS->PrintVar[26],DS,cD->etInt,DS->NumEle,t,25,CV_Y);
		avgResults_MD(outp[27],DS->PrintVar[27],DS,cD->etInt,DS->NumEle,t,26,CV_Y);
		avgResults_MD(outp[28],DS->PrintVar[28],DS,cD->etInt,DS->NumEle,t,27,CV_Y);
		avgResults_MD(outp[29],DS->PrintVar[29],DS,cD->etInt,DS->NumEle,t,28,CV_Y);
		avgResults_MD(outp[30],DS->PrintVar[30],DS,cD->etInt,DS->NumEle,t,29,CV_Y);
		avgResults_MD(outp[31],DS->PrintVar[31],DS,cD->etInt,DS->NumEle,t,30,CV_Y);
		avgResults_MD(outp[32],DS->PrintVar[32],DS,cD->etInt,DS->NumEle,t,31,CV_Y);
		avgResults_MD(outp[33],DS->PrintVar[33],DS,cD->etInt,DS->NumEle,t,32,CV_Y);
		}
	}
  
