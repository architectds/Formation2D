/********************************************************************************
 * File        : swc.c                                    	                *
 * Function    : for calculation of soil water content		       	        *
 * Version     : May, 2009                                                      *
 * Developer   : Yuning Shi (yshi@psu.edu)                                      *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                   *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)              *
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"      

#define EPSILON 0.05
#define MINpsi	-70
#define multF	2
#define UNIT_C 1440     /* 60*24 for calculation of yDot in m/min units while forcing is in m/day. */

realtype effKV(realtype ksatFunc,realtype gradY,realtype macKV,realtype KV,realtype areaF);

void swc(realtype t, realtype stepsize, void *DS, N_Vector VY, Control_Data *CS)
	{
  	int i,j,k,checkele=1;
	int bot;
	realtype p[5],q[5],a[5],c[5],r[5],b[5],m[5],beta[5],s[5],F[5],D[5],K[5],Eroot[5],z[5],dz[5],dz0[5],temp[5],FX[5],SFX,rtx,denom,Esoil,SM[5];
	realtype mm;
	realtype ThetaR, ThetaS;
	realtype ThetaBar, ThetaTot;
	realtype dPsi_dTheta, elemSatn;
	realtype R, I;
	realtype droot, hsw, h, h1, dmac;
	realtype AquiferDepth;
	realtype *Y;
	realtype Deficit, Avg_Y_Sub, TotalY_Ele, Grad_Y_Sub, effK, satKfunc;
	realtype zsat, zunsat;			// saturated depth and unsaturated depth of the bottom layer

  	Y = NV_DATA_S(VY);
  	Model_Data MD;
  	MD = (Model_Data)DS;
 	stepsize=stepsize/UNIT_C;

  	for(i=0; i<MD->NumEle; i++)
  		{
		AquiferDepth=MD->Ele[i].zmax-MD->Ele[i].zmin;
		if (AquiferDepth<MD->Ele[i].macD) MD->Ele[i].macD = AquiferDepth;

		ThetaS = MD->Ele[i].ThetaS;
		ThetaR = MD->Ele[i].ThetaR;
		elemSatn = (AquiferDepth-Y[i+2*MD->NumEle]<0)?1.0:(Y[i+MD->NumEle]<0?0:Y[i+MD->NumEle]/(AquiferDepth-(Y[i+2*MD->NumEle]>0?Y[i+2*MD->NumEle]:0)));
		elemSatn = elemSatn>1.0?1.0:(elemSatn<0?0:elemSatn); 
		ThetaBar = elemSatn*(ThetaS - ThetaR)+ThetaR;
		h = Y[i+2*MD->NumEle]>0?Y[i+2*MD->NumEle]:0;
//		if (i==checkele) printf("\n AquiferDepth = %f, h = %f, elemSatn = %f", AquiferDepth, h, elemSatn);

		for (j=0; j<5;j++)
			{
			dz0[j] = MD->dsoil[i][j];
			}

		if (elemSatn == 0)		// If soil water is 0, soil moisture equals to residual moisture
			{
			for (j = 0; j<5; j++)
				{
				if (MD->dsoil[i][j]>0)
					{
					MD->EleSM[i][j] = ThetaR;
					bot = j;
					}
				else
					{
					MD->EleSM[i][j] = -999.0;
					}
				}

			MD->EleBot[i] = bot;
			}
		else if (elemSatn == 1.0)
			{
			for (j = 0; j<5; j++)
				{
				if (MD->dsoil[i][j]>0)
					{
					MD->EleSM[i][j] = ThetaS;
					bot = j;
					}
				else
					{
					MD->EleSM[i][j] = -999.0;
					}
				}
			}
		else if (AquiferDepth - h <= dz0[0])
			{
			zunsat = AquiferDepth - h;
			zsat = dz0[0] - zunsat;
			zsat = zsat<0?0:zsat;
			MD->EleSM[i][0] = (zunsat*ThetaBar + zsat*ThetaS)/(zunsat+zsat);
			for (j = 1; j<5; j++)
				{
				if (MD->dsoil[i][j]>0)
					{
					MD->EleSM[i][j] = ThetaS;
					}
				else
					{
					MD->EleSM[i][j] = -999.0;
					}
				}
			bot = 0;
			MD->EleBot[i] = 0;
//			if (i==checkele) printf("\n ThetaBar = %f, MD->EleSM[i][0] = %f", ThetaBar, MD->EleSM[i][0]);
			}
		else if (CS->init_type<3)
			{
			hsw = dz0[0]+dz0[1];
			k = 2;

			while(h+hsw < AquiferDepth & dz0[k]>0)
				{
				hsw = hsw+dz0[k];
				k++;
				}
			bot = k-1;

			bot = bot>4?4:bot;
			zsat = hsw - (AquiferDepth - h);
			zunsat = dz0[bot] - zsat;
			for (j = 0; j<bot; j++)
				{
				MD->EleSM[i][j] = ThetaBar;
				}
			MD->EleSM[i][bot] = (ThetaBar*zunsat + ThetaS*zsat)/(zunsat+zsat);
			for (j = bot+1; j<5; j++)
				{
				if (MD->dsoil[i][j]>0)
					{
					MD->EleSM[i][j] = ThetaS;
					}
				else
					{
					MD->EleSM[i][j] = -999.0;
					}			
				}

			}
		else
			{
/*
			Deficit=AquiferDepth-Y[i+2*MD->NumEle];
			elemSatn=((Y[i+MD->NumEle]/Deficit)>1)?1:((Y[i+MD->NumEle]<=0)?EPSILON/1000.0:Y[i+MD->NumEle]/Deficit);
			elemSatn=(elemSatn<multF*EPSILON)?multF*EPSILON:elemSatn;
			Avg_Y_Sub=(-(pow(pow(1/elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha)<MINpsi)?MINpsi:-(pow(pow(1/elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha);
			TotalY_Ele=Avg_Y_Sub+MD->Ele[i].zmin+AquiferDepth-MD->Ele[i].infD;
			Grad_Y_Sub=(Y[i]+MD->Ele[i].zmax-TotalY_Ele)/MD->Ele[i].infD;
//			if (i==30) printf("\n Grad_Y_Sub = %f, Y[i] = %f", Grad_Y_Sub, Y[i]);
			Grad_Y_Sub=((Y[i]<EPSILON/100)&&(Grad_Y_Sub>0))?0:Grad_Y_Sub;
			satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,MD->Ele[i].Beta/(MD->Ele[i].Beta-1)),(MD->Ele[i].Beta-1)/MD->Ele[i].Beta),2);
//			satKfunc = satKfunc<0.13?0.13:satKfunc;
			effK=(MD->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y_Sub,MD->Ele[i].macKsatV,MD->Ele[i].infKsatV,MD->Ele[i].hAreaF):MD->Ele[i].infKsatV;
			MD->Ele[i].effKV = effK;
//     			MD->EleViR[i] = 0.5*(effK+MD->Ele[i].infKsatV)*Grad_Y_Sub;
			I = 0.5*(effK)*Grad_Y_Sub;
//			if (i==30 && I>0) printf("\n Grad_Y_Sub = %f, satKfunc = %f, effK = %f, I = %f, sfc = %f", Grad_Y_Sub, satKfunc, effK, I, Y[i]);
//			I = I/CS->Cal.Porosity;
//	                MD->Recharge[i] = (elemSatn==0.0)?0:(Deficit<=0)?0:(MD->Ele[i].KsatV*satKfunc*(MD->Ele[i].Alpha*Deficit-2*pow(-1+pow(elemSatn,MD->Ele[i].Beta/(-MD->Ele[i].Beta+1)),1/MD->Ele[i].Beta))/(MD->Ele[i].Alpha*((Deficit+MD->DummyY[i+2*MD->NumEle]*satKfunc))));
			effK=(MD->Ele[i].Macropore==1)?((Y[i+2*MD->NumEle]>AquiferDepth-MD->Ele[i].macD)?effK:MD->Ele[i].KsatV*satKfunc):MD->Ele[i].KsatV*satKfunc;
	                R = (elemSatn==0.0)?0:(Deficit<=0)?0:(MD->Ele[i].KsatV*Y[i+2*MD->NumEle]+effK*Deficit)*(MD->Ele[i].Alpha*Deficit-2*pow(-1+pow(elemSatn,MD->Ele[i].Beta/(-MD->Ele[i].Beta+1)),1/MD->Ele[i].Beta))/(MD->Ele[i].Alpha*pow(Deficit+Y[i+2*MD->NumEle],2));
			R=(R>0 && Y[i+MD->NumEle]<=0)?0:R;//??BHATT
			R=(R<0 && Y[i+2*MD->NumEle]<=0)?0:R;//??BHATT
//			R = R/CS->Cal.Porosity;
*/

//			I = MD->EleViR[i];

			I = (Y[i+MD->NumEle]-MD->Eleunsatbuffer[i])*MD->Ele[i].Porosity/stepsize+MD->Recharge[i]+(MD->EleOVLbuffer[i]<EPSILON/100?MD->EleET[i][2]:0)+((MD->EleGWbuffer[i]<AquiferDepth-MD->Ele[i].RzD)?MD->EleET[i][1]:0);	
			I = I>0?I:0;

			R = MD->Recharge[i];
//			printf("\n 1/stepsize = %f", 1/stepsize);
//			R = R - ThetaS/stepsize*(Y[i+2*MD->NumEle] - MD->EleGWbuffer[i])/(AquiferDepth - MD->EleGWbuffer[i])*MD->Eleunsatbuffer[i];
//			if (I>0) printf("\n Ele %d infiltration  = %f", i, I);

//			if (MD->EleViR[i]>0) printf("\n Ele %d\tTrasfer I = %f\tMass balance I = %f, Recharge = %f", i, MD->EleViR[i], I, R);

//			if (i==30) printf("\n swc.c infiltration = %f, Recharge = %f, adjustment = %f", I*CS->Cal.Porosity, R*CS->Cal.Porosity, - ThetaS/stepsize*(Y[i+2*MD->NumEle] - MD->EleGWbuffer[i])/(AquiferDepth - MD->EleGWbuffer[i])*MD->Eleunsatbuffer[i]);
			// Rescale Recharge and infiltration

//			R = MD->Recharge[i]/CS->Cal.Porosity;
//			I = MD->EleViR[i]/CS->Cal.Porosity;

			// Looking for the deepest layer above water table

			hsw = dz0[0]+dz0[1];
			k = 2;

			while(h+hsw < AquiferDepth & dz0[k]>0)
				{
				hsw = hsw+dz0[k];
				k++;
				}
			bot = k-1;

			bot = bot>4?4:bot;
			zsat = hsw - (AquiferDepth - h);
			zunsat = dz0[bot] - zsat;
//			zsat = zsat<0?0:zsat;

			if (zunsat < 0.025)
				{
				bot = bot - 1;
				zunsat = dz0[bot];
				zsat = 0;
				}
			
			z[0] = 0.5*dz0[0];
			for (j=1; j<bot+1; j++)
				{
				z[j] = z[j-1] + 0.5*dz0[j] + 0.5*dz0[j-1];
				}

			/* Simulate soil moisture */
			for (j = 0; j<bot; j++)
				{
				dz[j] = dz0[j];
				SM[j] = MD->EleSM[i][j];
				SM[j] = SM[j]<=ThetaR?(ThetaR + EPSILON/100.0):SM[j];
				SM[j] = SM[j]>=ThetaS?(ThetaS - EPSILON/100.0):SM[j];
				}
			dz[bot] = zunsat;
			SM[bot] = (MD->EleSM[i][bot]*(zunsat+zsat) - ThetaS*zsat)/zunsat;
			SM[bot] = SM[bot]<=ThetaR?(ThetaR + EPSILON/100.0):SM[bot];
			SM[bot] = SM[bot]>=ThetaS?(ThetaS - EPSILON/100.0):SM[bot];

//			if(i==checkele) printf("\n ThetaS = %f, SM[bot]_bar = %f, SM[bot] = %f", ThetaS, MD->EleSM[i][bot], SM[bot]);

//			if (i==checkele) printf("\n Ele = %d, AquiferDepth = %f, AD - h = %f, bot = %d, hsw = %f, zbot = %f, zsat = %f, zunsat = %f , dz[bot] = %f", i, AquiferDepth, AquiferDepth-h, bot, hsw, dz0[bot], zsat, zunsat, dz[bot]);

			for (j=0;j<5;j++)
				{
				Eroot[j] = 0;
				}

			if(h<AquiferDepth-MD->Ele[i].RzD)
				{
				droot = 0;
				SFX = 0;
				for (j=0; j<=MD->Ele[i].rootL; j++)
					{
					FX[j] = (MD->EleSM[i][j] - MD->Ele[i].ThetaW)/(MD->Ele[i].ThetaRef - MD->Ele[i].ThetaW);
	//				if (i==checkele) printf("\n MD->EleSM[%d][%d] = %f, ThetaRef = %f, ThetaW = %f, FX[%d] = %f", i, j, MD->EleSM[i][j], MD->Ele[i].ThetaRef, MD->Ele[i].ThetaW, j, FX[j]);
	//				FX[i] = (FX[i]>1)?1:(FX[i]<0?0:FX[i]);
					SFX = SFX + FX[j];
					droot = droot+MD->dsoil[i][j];
					}
				SFX = SFX / (MD->Ele[i].rootL + 1);

				denom = 0;
				for (j=0; j<=MD->Ele[i].rootL; j++)
					{
					rtx = (MD->dsoil[i][j]/droot) + FX[j] - SFX;
	//				if (i==checkele) printf("\n droot = %f, rtx = %f", droot, rtx);
					FX[j] = FX[j]*(rtx<0?0:rtx);
					denom = denom + FX[j];
					}
				denom = denom<=0?1:denom;
				
				for (j=0; j<=MD->Ele[i].rootL; j++)
					{
					Eroot[j] = FX[j]*MD->EleET[i][1]/denom;
					}
//				if (i == checkele) printf("\n Ele %d: denom = %f, FX[0] = %f, FX[1] = %f, FX[2] = %f, FX[3] = %f, Tr = %f, Tr[0] = %f, Tr[1] = %f, Tr[2] = %f, Tr[3] = %f", i, denom, FX[0], FX[1], FX[2], FX[3], MD->EleET[i][1], Eroot[0], Eroot[1], Eroot[2], Eroot[3]);
				}	

//			if(h<AquiferDepth-MD->Ele[i].RzD) Eroot[MD->Ele[i].rootL] = MD->EleET[i][1];
//			if(h<AquiferDepth-MD->Ele[i].RzD) Eroot[MD->Ele[i].rootL] = MD->EleET[i][1]/CS->Cal.Porosity;

//			printf("\n droot = %f, nroot = %d", droot, nroot);

//			Esoil = Y[i]<EPSILON/100?MD->EleET[i][2]/CS->Cal.Porosity:0;
			Esoil = Y[i]<(EPSILON/100)?MD->EleET[i][2]:0;

			for (j = 0;j<bot;j++)
				{
				satKfunc=pow((SM[j]-ThetaR)/(ThetaS-ThetaR),0.5)*pow(-1+pow(1-pow((SM[j]-ThetaR)/(ThetaS-ThetaR),MD->Ele[i].Beta/(MD->Ele[i].Beta-1)),(MD->Ele[i].Beta-1)/MD->Ele[i].Beta),2);
//				satKfunc = satKfunc<0.13?0.13:satKfunc;
//				if (i==checkele) printf("\n satKfunc = %f", satKfunc);

//				K[j] = (MD->Ele[i].Macropore==1 && j<= MD->Ele[i].MPL)?MD->Ele[i].effKV:satKfunc*MD->Ele[i].KsatV;

//				if  (i==checkele) printf("\n layer = %d, effK = %f, Kr = %f", j, MD->Ele[i].effKV, satKfunc*MD->Ele[i].KsatV);
//				K[j] = MD->Ele[i].macKsatV*MD->Ele[i].hAreaF*satKfunc+(1-MD->Ele[i].hAreaF)*satKfunc*MD->Ele[i].KsatV;
				K[j] = MD->Ele[i].KsatV*satKfunc;
/*
				Y_Sub1=-(pow(pow((SM[j]-ThetaR)/(ThetaS-ThetaR),-MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha);
				Y_Sub2=-(pow(pow((SM[j+1]-ThetaR)/(ThetaS-ThetaR),-MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha);
				Grad_Y_Sub=(Y_Sub1-Y_Sub2+0.5*dz[j]+0.5*dz[j+1])/dz[j];
				K[j] = (MD->Ele[i].Macropore==1 && j<= MD->Ele[i].MPL)?effKV(satKfunc,Grad_Y_Sub,MD->Ele[i].macKsatV,MD->Ele[i].KsatV,MD->Ele[i].hAreaF):satKfunc*MD->Ele[i].KsatV;
*/
//				K[j] = satKfunc*MD->Ele[i].KsatV*exp(-1.8*(z[j]-1.0));
//				macKV*areaF*ksatFunc+KV*(1-areaF)*ksatFunc
//				if (i==checkele) printf("\n effK = %f, K = %f, macKV = %f, areaF = %f, KV = %f, infKV = %f", MD->Ele[i].effKV, K[j], MD->Ele[i].macKsatV, MD->Ele[i].hAreaF, MD->Ele[i].KsatV, MD->Ele[i].infKsatV);
				mm = (MD->Ele[i].Beta-1)/MD->Ele[i].Beta;
				dPsi_dTheta = (1-mm)/MD->Ele[i].Alpha/mm/(ThetaS-ThetaR)*pow(pow((SM[j]-ThetaR)/(ThetaS-ThetaR),-1.0/mm)-1,0.0-mm)*pow((SM[j]-ThetaR)/(ThetaS-ThetaR),-(1.0+mm)/mm);
				dPsi_dTheta = dPsi_dTheta<0?(-dPsi_dTheta):dPsi_dTheta;
				D[j] = dPsi_dTheta*K[j];
//				D[j] = (1-(MD->Ele[i].Beta-1)/MD->Ele[i].Beta)*MD->Ele[i].KsatV/(MD->Ele[i].Alpha*(ThetaS-ThetaR)*(MD->Ele[i].Beta-1)/MD->Ele[i].Beta)*pow((SM[j]-ThetaR)/(ThetaS-ThetaR),0.5-MD->Ele[i].Beta/(MD->Ele[i].Beta-1));
//				D[j] = D[j]*(pow(1-pow((SM[j]-ThetaR)/(ThetaS-ThetaR),MD->Ele[i].Beta/(MD->Ele[i].Beta-1)),-(MD->Ele[i].Beta-1)/MD->Ele[i].Beta)+pow(1-pow((SM[j]-ThetaR)/(ThetaS-ThetaR),MD->Ele[i].Beta/(MD->Ele[i].Beta-1)),(MD->Ele[i].Beta-1)/MD->Ele[i].Beta)-2);
//				if (i==checkele) printf("\n D = %f", D[j]);
//				D[j] = D[j]<0?-D[j]:D[j];
//				D[j] = K[j]*dPsi_dTheta;
				}
/*
			if (bot>MD->EleBot[i] & MD->Recharge[i]<0)
				{
				F[0] = -K[0]+MD->EleViR[i]-Esoil-Eroot[0];
				for (j = 1; j<bot-1;j++)
					{
					F[j] = K[j-1]-K[j]-Eroot[j];
					}
				F[bot-1] = K[bot-2] - K[bot-1]-Eroot[bot-1]-MD->Recharge[i];
				F[bot] = K[bot-1] + MD->Recharge[i]-Eroot[bot];
				}
			else
				{
				F[0] = -K[0]+MD->EleViR[i]-Esoil-Eroot[0];
				for (j = 1; j<bot;j++)
					{
					F[j] = K[j-1]-K[j]-Eroot[j];
					}
				F[bot] = K[bot-1] - MD->Recharge[i]-Eroot[bot];
				}
*/

			F[0] = -K[0]+I-Esoil-Eroot[0];
			for (j = 1; j<bot;j++)
				{
				F[j] = K[j-1]-K[j]-Eroot[j];
				}
			F[bot] = K[bot-1] - R -Eroot[bot];
			MD->EleBot[i] = bot;

			p[0] = 0;
			q[0] = stepsize*D[0]/2/dz[0]/(0.5*dz[0]+0.5*dz[1]);

			for(j=1;j<bot;j++)
				{
				p[j]=stepsize*D[j-1]/2/dz[j]/(0.5*dz[j-1]+0.5*dz[j]);
				q[j]=stepsize*D[j]/2/dz[j]/(0.5*dz[j]+0.5*dz[j+1]);
				}

			p[bot] = stepsize*D[bot-1]/2/dz[bot]/(0.5*dz[bot-1]+0.5*dz[bot]);
			q[bot] = 0;

			a[0]=p[0];
			c[0]=q[0];
			b[0]=-(p[0]+q[0]+1);
			r[0]=(q[0]-1)*MD->EleSM[i][0]-q[0]*MD->EleSM[i][1]-F[0]*stepsize/dz[0];

			for(j=1;j<bot;j++)
				{
				a[j]=p[j];
				c[j]=q[j];
				b[j]=-(p[j]+q[j]+1);
				r[j]=-p[j]*MD->EleSM[i][j-1]+(p[j]+q[j]-1)*MD->EleSM[i][j]-q[j]*MD->EleSM[i][j+1]-F[j]*stepsize/dz[j];
				}

			a[bot]=p[bot];
			c[bot]=q[bot];
			b[bot]=-(p[bot]+q[bot]+1);
			r[bot]=-p[bot]*MD->EleSM[i][bot-1]+(p[bot]+q[bot]-1)*MD->EleSM[i][bot]-F[bot]*stepsize/dz[bot];

			beta[0]=b[0];
			s[0]=r[0];

			for(j=1;j<=bot;j++)
				{
				m[j]=a[j]/beta[j-1];
				beta[j]=b[j]-m[j]*c[j-1];
				s[j]=r[j]-m[j]*s[j-1];
				}

			temp[bot]=s[bot]/beta[bot];

			for(j=bot-1;j>=0;j--)
				{
				temp[j]=(s[j]-c[j]*temp[j+1])/beta[j];
				}

			ThetaTot = 0;

			for(j = 0;j<=bot;j++)
				{
				ThetaTot = temp[j]*dz[j]/hsw + ThetaTot;
				}
/*
			for(j = 0;j<=bot;j++)
				{
				MD->EleSM[i][j] = temp[j]/ThetaTot*ThetaBar+ThetaR;
//				MD->EleSM[i][j] = temp[j]+ThetaR;
				}
*/
			for (j = 0; j<bot; j++)
				{
				MD->EleSM[i][j] = temp[j];
				}
			MD->EleSM[i][bot] = (temp[bot]*zunsat + ThetaS*zsat)/(zunsat+zsat);
//			MD->EleSM[i][bot] = ((temp[bot] + ThetaR)*zunsat + ThetaS*zsat)/(zunsat+zsat);
//			MD->EleSM[i][bot] = temp[bot];
			for (j = bot+1; j<5; j++)
				{
				if (MD->dsoil[i][j]>0)
					{
					MD->EleSM[i][j] = ThetaS;
					}
				else
					{
					MD->EleSM[i][j] = -999.0;
					}			
				}
			}


		for (j = 0; j<5; j++)
			{
			if (MD->dsoil[i][j]>0)
				{	
				MD->EleSM[i][j] = MD->EleSM[i][j]>ThetaS?ThetaS:MD->EleSM[i][j];
				MD->EleSM[i][j] = MD->EleSM[i][j]<ThetaR?ThetaR:MD->EleSM[i][j];
				}
			}
/*
			if ( i == checkele)
				{
				for (j = 0; j<bot; j++)
					{
					printf("\n K%d,D%d,I%d = %f, %f, %f",j, K[j],j,D[j],j,F[j]);
					}
				printf("\n I%d = %f",bot, F[j]);
				}


			if (i == checkele)
				{
				if (bot>0) printf("\n K = %f, infiltration = %f, Recharge = %f, Eroot = %f", K[bot-1] , MD->EleViR[i], MD->Recharge[i],Eroot[bot]);
				printf("\n Ele = %d, H = %f, s = %f, h = %f, H - h = %f, hsw = %f, zsat = %f, zunsat = %f, unsatSM = %f, bot = %d, bot = %d", i, AquiferDepth, Y[i+MD->NumEle], h, AquiferDepth - h, hsw, zsat, zunsat, temp[bot], MD->EleBot[i], bot);
				printf("\n Ele = %d, eleSatn = %f, SM = %f, %f, %f, %f, %f, ThetaS = %f",i,elemSatn, MD->EleSM[i][0],MD->EleSM[i][1],MD->EleSM[i][2],MD->EleSM[i][3],MD->EleSM[i][4], ThetaS);
				}
*/
		MD->EleOVLbuffer[i] = Y[i];
		MD->EleGWbuffer[i] = Y[i+2*MD->NumEle];
		MD->Eleunsatbuffer[i] = Y[i+MD->NumEle];

		}
	}
