/********************************************************************************
 * File        : seb.c                                    	                *
 * Function    : for calculation of surface energy balance       	        *
 * Version     : May, 2009                                                      *
 * Developer   : Yuning Shi (yshi@psu.edu)                                      *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                   *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)              *
 *------------------------------------------------------------------------------*
 * Contains is_sm_et.c of PIHM2.0                                               *
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"      

#define EPSILON 0.05
#define UNIT_C 1440     /* 60*24 for calculation of yDot in m/min units while forcing is in m/day. */
#define multF2	5	/* Melting factor */
#define GRAV 9.8*60*60	/* Note the dependence on physical units */
#define MAXLOOP 1

#define C_water 4218000.0
#define C_soil 2000000.0
#define C_air 1004.0
#define Lv 2503000.0
#define SIGMA (5.67*pow(10,-8)*60)
#define R_dry 287.04
#define R_v 461.5

realtype Interpolation(TSD *Data, realtype t);

realtype pslmu(realtype z)
	{
	realtype x;
	x = -0.96*log(1.0-4.5*z);
	return x;
	}

realtype pslms(realtype z)
	{
	realtype x;
	x = z/0.183-2.076*(1.0-1.0/(z+1.0));
	return x;
	}

realtype pslhu(realtype z)
	{
	realtype x;
	x = -0.96*log(1-4.5*z);
	return x;
	}

realtype pslhs(realtype z)
	{
	realtype x;
	x = z*0.183/(0.8*0.191*0.191)-2.076*(1.0-exp(-1.2*z));
	return x;
	}

realtype pspmu(realtype z)
	{
	realtype x;
	x = -2.0*log((z+1.0)*0.5)-log((z*z+1)*0.5)+2.0*atan(z)-3.1416/2.0;
	return x;
	}

realtype pspms(realtype z)
	{
	realtype x;
	x = 5.0*z;
	return x;
	}

realtype psphu(realtype z)
	{
	realtype x;
	x = -2.0*log((z*z+1.0)*0.5);
	return x;
	}

realtype psphs(realtype z)
	{
	realtype x;
	x = 5.0*z;
	return x;
	}

void excoef(void *DS, realtype T, realtype Tv, realtype Ts, realtype zlvl, realtype Vel, realtype P, realtype rl, int i);

void seb(realtype t, realtype stepsize, void *DS, N_Vector VY)
	{
  	int i,j,k,kk,inabr,ebloop,nroot;
  	realtype totEvap;
  	realtype Soldown, Fdown, LWdown, G, LE, LEp, T, Tv, theta, Vel, RH, VP, P, LAI, Rn, RivPrep;
	realtype zero_dh, cnpy_h, rl, r_a, r_s, f_r, h_s, FX, Rmax, z, Rr, P_c, qv, qv_sat, droot, VegFrac, VegFrac_min;
	realtype F1,F2,F3,F4, avg_F1, avg_F2, avg_F3, avg_F4;
	realtype alb_brg, alb_veg;
	realtype emiss_brg, emiss_veg, emiss;
  	realtype isval=0,etval=0;
	realtype zlvl;
  	realtype fracSnow,snowRate,MeltRateGrnd,MeltRateCanopy,eltRate,MF,Ts=-3.0,Tr=1.0,To=0.0,ret;
	realtype fraction;
	realtype Terror;
	realtype p[5],q[5],a[5],c[5],r[5],b[5],m[5],beta[5],s[5],kappa[5],Cs[5],dz[5],soiltemp[5],SM[5];
	realtype Pf;
  	realtype Delta, Gamma, oldDelta, oldGamma;
	realtype AquiferDepth;
	realtype ThetaS, ThetaR;
	realtype extcoef, sbeta;

	realtype Tsnew,Tsold;

	realtype ET0,ET1,ET2;

	realtype A1,A2,B1,B2;

  	Model_Data MD;
  	MD = (Model_Data)DS;
 	stepsize=stepsize/UNIT_C;

	extcoef = -1.0;
	sbeta = -4.0;

  	for(i=0; i<MD->NumEle; i++)
  		{
		/* Read level of wind speed measurements */

		zlvl = MD->Ele[i].windH;

		/* Find the deepest layer which root reaches */

		k = 0;
		nroot = 0;
		droot = 0;

		while(droot<MD->Ele[i].RzD)
			{
			droot = droot+MD->dsoil[i][k];
			nroot = k;
			k++;
			}

		/* Read forcing data */

		MD->ElePrep[i] = Interpolation(&MD->TSD_Prep[MD->Ele[i].prep-1], t);
		Soldown = Interpolation(&MD->TSD_Sdown[MD->Ele[i].Sdown-1], t);
		T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp-1], t);	
		Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel-1], t);
		LWdown = Interpolation(&MD->TSD_Ldown[MD->Ele[i].Ldown-1], t);
		RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity-1], t);
		P = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure-1], t);
		LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC-1], t);
    		MF = multF2*Interpolation(&MD->TSD_MeltF[MD->Ele[i].meltF-1], t);
//		rl = Interpolation(&MD->TSD_RL[MD->Ele[i].LC-1], t); 

		/* Calculate green vegetation fraction as a function of LAI (Niu et al. 2011 Noah-MP), parameters are slightly changed */


		VegFrac = MD->Ele[i].VegFrac*(1-exp(extcoef*LAI));
		VegFrac_min = MD->Ele[i].VegFrac*(1-exp(extcoef*0.5));

		if (VegFrac<0)
			{
			printf("\n Error! Vegetation fraction < 0!");
			}

//		VegFrac = MD->Ele[i].VegFrac;

		/* Initialize surface skin temperature */

		Tsnew = MD->Tsfc[i];
		Tsold = MD->Tsfc[i];

		/* Calculate moisture variables and temperature variables */

		VP = 611.2*exp(17.67*T/(T+243.5))*RH;
		qv_sat = 0.622*VP/RH/(P-(1.0-0.622)*VP/RH);
		qv = 0.622*VP/(P-(1.0-0.622)*VP);
		theta = T+0.0098*2.0;
		Tv = (T+273.15)*(1.0+0.61*qv)-273.15;
		/* Albedo, emissivity and roughness length vary from minimum to maximum as a function of GVF */

		fraction = (VegFrac - VegFrac_min)/(MD->Ele[i].VegFrac - VegFrac_min);
		fraction = (fraction>1)?1:(fraction<0?0:fraction);

  		MD->Ele[i].Albedo = (MD->EleSnowGrnd[i] > EPSILON)?0.60:(MD->Ele[i].Albedo_max + fraction*(MD->Ele[i].Albedo_min - MD->Ele[i].Albedo_max));

		emiss = (MD->EleSnowGrnd[i] > EPSILON)?0.95:(MD->Ele[i].Emiss_min + fraction*(MD->Ele[i].Emiss_max - MD->Ele[i].Emiss_min));

		rl = MD->Ele[i].z0_min + fraction*(MD->Ele[i].z0_max - MD->Ele[i].z0_min);

		/* Total downward radiation */

		Fdown = Soldown*(1-MD->Ele[i].Albedo)+LWdown;

		/* Initialzize snow distribution */

		MD->EleSnowGrnd[i] = (1-VegFrac) * MD->EleSnow[i];
		MD->EleSnowCanopy[i] = VegFrac * MD->EleSnow[i];

		/******************************************************************************************/
		/*			    Snow Accumulation/Melt Calculation				  */
		/******************************************************************************************/
    		fracSnow = T<Ts?1.0:(T>Tr?0:(Tr-T)/(Tr-Ts));
    		snowRate = fracSnow*MD->ElePrep[i];
		/* EleSnowGrnd, EleSnowCanopy, EleISsnowmax, MeltRateGrnd,MeltRateCanopy are the average value prorated over the whole elemental area */
    		MD->EleSnowGrnd[i]=MD->EleSnowGrnd[i]+(1-VegFrac)*snowRate*stepsize;
    		MD->EleSnowCanopy[i]=MD->EleSnowCanopy[i]+VegFrac*snowRate*stepsize;
		MD->EleISsnowmax[i]=MD->EleSnowCanopy[i]>0?0.003*LAI*VegFrac:0;
		MD->EleISsnowmax[i]=2*MD->IS*MD->EleISsnowmax[i];	// Shi (temp)
		if(MD->EleSnowCanopy[i]>MD->EleISsnowmax[i])
			{
			MD->EleSnowGrnd[i]=MD->EleSnowGrnd[i]+MD->EleSnowCanopy[i]-MD->EleISsnowmax[i];
			MD->EleSnowCanopy[i]=MD->EleISsnowmax[i];
			}
    		MeltRateGrnd=MeltRateCanopy=(T>To?(T-To)*MF:0);		/* Note the units for MF. */
    		if(MD->EleSnowGrnd[i]>MeltRateGrnd*stepsize)
    			{
    			MD->EleSnowGrnd[i]=MD->EleSnowGrnd[i]-MeltRateGrnd*stepsize;
    			}
    		else
    			{
    			MeltRateGrnd=MD->EleSnowGrnd[i]/stepsize;
    			MD->EleSnowGrnd[i]=0;    	
    			} 
                if(MD->EleSnowCanopy[i]>MeltRateCanopy*stepsize)
                        {
                        MD->EleSnowCanopy[i]=MD->EleSnowCanopy[i]-MeltRateCanopy*stepsize;
                        }
                else
                        {
                        MeltRateCanopy=MD->EleSnowCanopy[i]/stepsize;
                        MD->EleSnowCanopy[i]=0;
                        }

		/************************************************************************/
		/*		ThroughFall and Evaporation from canopy			*/
		/************************************************************************/
		/* EleIS, EleET[0] and ret are prorated for the whole element. Logistics are simpler if assumed in volumetric form by multiplication of Area on either side of equation*/

		MD->EleISmax[i] = MD->IS*MD->ISFactor[MD->Ele[i].LC-1]*LAI*VegFrac;

    		if(LAI>0.0)
    			{	 
	    		MD->EleTF[i]=MD->EleIS[i]<=0?0:(5.65*pow(10,-2)*MD->EleISmax[i]*exp(3.89*(MD->EleIS[i]<0?0:MD->EleIS[i])/MD->EleISmax[i])); /* Note the dependece on physical units*/
			MD->EleTF[i]=MD->TF*MD->EleTF[i];
			}
     		else
     			{
			MD->EleTF[i]=0.0;
     			}

		/************************************************************************/
		/*		Surface energy balance loops    			*/
		/************************************************************************/

		ebloop = 0;
		Terror = 999.0;

		ThetaS = MD->Ele[i].ThetaS;
		ThetaR = MD->Ele[i].ThetaR;

		while(ebloop<MAXLOOP && Terror>0.005)
			{
			ebloop++;

			/* Calculate C_h and C_m */
			excoef(MD, T, Tv, Tsnew, zlvl, Vel, P, rl, i);

			/* Simulate soil temperature */

			/* Define dz of each level for soil temp simulation. */

			for(j=0;j<3;j++)
				{
				dz[j] = (MD->dsoil[i][j]+MD->dsoil[i][j+1])/2;
				}
			dz[3] = MD->dsoil[i][3]/2+1.0;

			/* Calculate thermal conductivity for each layer */

			for(j=0;j<4;j++)
				{
				Cs[j] = MD->EleSM[i][j]*C_water+(1-ThetaS)*C_soil+(ThetaS-MD->EleSM[i][j])*C_air;
				SM[j] = MD->EleSM[i][j];
				SM[j] = (SM[j]<=ThetaR)?(ThetaR+EPSILON/100):SM[j];
				SM[j] = (SM[j]>=ThetaS)?(ThetaS-EPSILON/100):SM[j];
				Pf=log10(1500*pow(pow((ThetaS-ThetaR)/(SM[j]-ThetaR),MD->Ele[i].Beta/(MD->Ele[i].Beta-1))-1,1/MD->Ele[i].Beta)/MD->Ele[i].Alpha);
				kappa[j]=Pf>5.1?0.1744:(420*exp(-(2.7+Pf)));
				kappa[j]=kappa[j]>1.9?1.9:kappa[j];
				kappa[j]=kappa[j]*UNIT_C*60;
				}

			if (VegFrac>0)
				{
				kappa[0] = exp(sbeta*VegFrac)*kappa[0];
//				kappa[0] = (exp(-2.0*LAI/MD->Ele[i].LAImax)*VegFrac+1-VegFrac)*kappa[0];
//				kappa[0] = exp(sbeta*LAI)*kappa[0];
//				printf("\n LAI = %f, LAImax = %f, tau = %f",LAI,MD->Ele[i].LAImax, exp(-1.0*LAI/MD->Ele[i].LAImax));
				}

/*
			if(i == 223)
				{
				printf("\n\nstepsize= %f",stepsize);
				printf("\nkappa= %f",kappa[0]/UNIT_C/60);
				printf("\nkappa= %f",kappa[1]/UNIT_C/60);
				printf("\nkappa= %f",kappa[2]/UNIT_C/60);
				printf("\n SM = %f, %f, %f, %f", MD->EleSM[i][0],MD->EleSM[i][1],MD->EleSM[i][2],MD->EleSM[i][3]);
				printf("\nCs= %f",Cs[0]);
				printf("\nC_water= %f",C_water);
				printf("\nC_soil= %f",C_soil);
				printf("\nC_air= %f",C_air);
				printf("\nPorosity=%f ",MD->Ele[i].ThetaS);
				printf("\nSoil Moisture=%f ",MD->EleSM[i][0]);
				printf("\nBeta=%f ",MD->Soil[(MD->Ele[i].soil-1)].Beta);
				printf("\nPf=%f ",Pf);
				printf("\nCs=%f ",Cs[0]);
				}
*/

			/* Solve soil temperature equation following Cranck-Nicolson method */

			p[0] = stepsize*kappa[0]/Cs[0]/MD->dsoil[i][0]/MD->dsoil[i][0];
			q[0] = stepsize/2/Cs[0]/MD->dsoil[i][0]/dz[0]*((MD->dsoil[i][0]+MD->dsoil[i][1])*kappa[0]*kappa[1]/(MD->dsoil[i][0]*kappa[0]+MD->dsoil[i][1]*kappa[1]));
			a[0] = p[0];
			b[0] = -(p[0]+q[0]+1);
			c[0] = q[0];
			r[0] = -p[0]*MD->Tsfc[i]+(p[0]+q[0]-1)*MD->Tsoil[i][0]-q[0]*MD->Tsoil[i][1]-a[0]*Tsnew;
		
			for(j=1;j<3;j++)
				{
				p[j]=stepsize/2/Cs[j]/MD->dsoil[i][j]/dz[j-1]*((MD->dsoil[i][j-1]+MD->dsoil[i][j])*kappa[j-1]*kappa[j]/(MD->dsoil[i][j-1]*kappa[j-1]+MD->dsoil[i][j]*kappa[j]));
				q[j]=stepsize/2/Cs[j]/MD->dsoil[i][j]/dz[j]*((MD->dsoil[i][j]+MD->dsoil[i][j+1])*kappa[j+1]*kappa[j]/(MD->dsoil[i][j+1]*kappa[j+1]+MD->dsoil[i][j]*kappa[j]));
				a[j]=p[j];
				c[j]=q[j];
				r[j]=-p[j]*MD->Tsoil[i][j-1]+(p[j]+q[j]-1)*MD->Tsoil[i][j]-q[j]*MD->Tsoil[i][j+1];
				b[j]=-(p[j]+q[j]+1);
				}

			p[3]=stepsize/2/Cs[3]/MD->dsoil[i][3]/dz[2]*((MD->dsoil[i][3]+MD->dsoil[i][2])*kappa[3]*kappa[2]/(MD->dsoil[i][2]*kappa[2]+MD->dsoil[i][3]*kappa[3]));
			q[3]=stepsize/2/Cs[3]/MD->dsoil[i][3]/dz[3]*kappa[3];

			a[3] = p[3];
			b[3] = -(p[3]+q[3]+1);
			c[3] = q[3];
			r[3] = -p[3]*MD->Tsoil[i][2]+(p[3]+q[3]-1)*MD->Tsoil[i][3]-q[3]*MD->Tbot[i]-c[3]*MD->Tbot[i];

			beta[0]=b[0];
			s[0]=r[0];
			for(j=1;j<4;j++)
				{
				m[j]=a[j]/beta[j-1];
				beta[j]=b[j]-m[j]*c[j-1];
				s[j]=r[j]-m[j]*s[j-1];
				}

			soiltemp[3]=s[3]/beta[3];

			for(j=2;j>=0;j--)
				{
				soiltemp[j]=(s[j]-c[j]*soiltemp[j+1])/beta[j];
				}

			/* Ground heat flux */

			G = kappa[0]/(MD->dsoil[i][0]/2)*(Tsnew-soiltemp[0]);

			A1 = P/R_dry/(Tv+273.15)*C_air*MD->C_h[i];
			B1 = -P/R_dry/(Tv+273.15)*C_air*MD->C_h[i]*T;

			A2 = kappa[0]/(MD->dsoil[i][0]/2);
			B2 = -kappa[0]/(MD->dsoil[i][0]/2)*soiltemp[0];

			/* Potential evaporation (Troen and Mahrt,1986) */
	
			Gamma = 4*emiss*SIGMA*UNIT_C*R_dry/C_air*pow(T+273.15,4)/(P*MD->C_h[i])+1;
			Delta = Lv*Lv*0.622/R_v/C_air/pow(T+273.15,2)*qv_sat;

			LEp =(Delta*(Fdown-emiss*SIGMA*UNIT_C*pow(T+273.15,4)-G+P/R_dry/(Tv+273.15)*C_air*MD->C_h[i]*(theta-T))+P/R_dry/(Tv+273.15)*Lv*MD->C_h[i]*Gamma*(qv_sat-qv))/(Delta+Gamma);
//			LEp =(Delta*(Fdown-SIGMA*UNIT_C*pow(T+273.15,4)-G+P/R_dry/(Tv+273.15)*C_air*MD->C_h[i]*(theta-T))+P/R_dry/(Tv+273.15)*Lv*MD->C_h[i]*Gamma*(qv_sat-qv))/(Delta+Gamma);
			MD->EleEp[i] = LEp/(1000.0*Lv);

			/* Calculate canopy resistence */

			if(LAI>0.0)
				{
				Rmax = 5000.0/(60*UNIT_C);		/* Unit day_per_m */
				f_r= 1.1*Soldown/(MD->Ele[i].Rs_ref*LAI);
				h_s= MD->Ele[i].h_s;
				F1= (f_r+(MD->Ele[i].Rmin/Rmax))/(1+f_r);
				F1=(F1<0.0001)?0.0001:F1;
			    	F3= 1- 0.0016*pow((MD->Tref-T),2);
				F3=(F3<0.0001)?0.0001:F3;
				F2=1/(1+h_s*(qv_sat - qv));
				F2=(F2<0.01)?0.01:F2;
				F4 = 0;

				for(kk=0;kk<=nroot;kk++)
					{
					FX = (MD->EleSM[i][kk]-MD->Ele[i].ThetaW)/(MD->Ele[i].ThetaRef-MD->Ele[i].ThetaW);
					FX = (FX>1)?1:(FX<0?0:FX);
					F4= F4+FX*MD->dsoil[i][kk]/droot;
					}
				F4=(F4<0.0001)?0.0001:F4;

		    		r_s=((MD->Ele[i].Rmin/(F1*F2*F3*F4*LAI))> Rmax)?Rmax:(MD->Ele[i].Rmin/(F1*F2*F3*F4*LAI));

				P_c = (1+Delta/Gamma)/(1+r_s*MD->C_h[i]+Delta/Gamma);

				FX = (MD->EleSM[i][0]-MD->Ele[i].ThetaW)/(MD->Ele[i].ThetaRef-MD->Ele[i].ThetaW);
				FX = (FX>1)?1:(FX<0?0:FX);
				ET2 = MD->pcCal.Et2*pow(FX,MD->fx_soil)*(1-VegFrac);
		    		ET1 = MD->pcCal.Et1*P_c*VegFrac*(1-pow(((MD->EleIS[i]+MD->EleSnowCanopy[i]<0)?0:(MD->EleIS[i]+MD->EleSnowCanopy[i]))/(MD->EleISmax[i]+MD->EleISsnowmax[i]),MD->fx_canopy));
	    			AquiferDepth=MD->Ele[i].zmax-MD->Ele[i].zmin;//??BHATT
				ET1 = ((MD->DummyY[i+2*MD->NumEle]<(AquiferDepth-MD->Ele[i].RzD))&&MD->DummyY[i+MD->NumEle]<=0)?0:ET1;//??BHATT
		    		ET0 = MD->pcCal.Et0*VegFrac*pow((MD->EleIS[i]<0?0:(MD->EleIS[i]>MD->EleISmax[i]?MD->EleISmax[i]:MD->EleIS[i]))/MD->EleISmax[i],MD->fx_canopy);

				if (LEp>0)
					{
					MD->EleET[i][0] = ET0*LEp/(1000.0*Lv);
					MD->EleDew[i] = 0;
					}
				else
					{
					MD->EleDew[i] = -LEp/(1000.0*Lv);
					MD->EleET[i][0] = 0;
					}

				/* ET0 is constrained by IS, needs to be updated */

				MD->EleIS[i] = MD->EleIS[i] + MD->EleDew[i]*stepsize;
				MD->EleETloss[i] = MD->EleET[i][0];

				if(MD->EleIS[i] >= MD->EleISmax[i])
					{
					if(((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)>=MD->EleET[i][0]+MD->EleTF[i])
						{
//						MD->EleETloss[i] = MD->EleET[i][0];
						ret = MD->EleTF[i]+(MD->EleIS[i]-MD->EleISmax[i])/stepsize+(((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)-(MD->EleETloss[i]+MD->EleTF[i]));
						isval=MD->EleISmax[i];	/* Modified by Y. Shi */
						}
     					else if((((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)<MD->EleET[i][0]+MD->EleTF[i])&&(MD->EleIS[i]+stepsize*((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy-MD->EleET[i][0]-MD->EleTF[i])<=0))
						{
     						MD->EleET[i][0]=(MD->EleETloss[i]/(MD->EleETloss[i]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy));
						ret =(MD->EleTF[i]/(MD->EleETloss[i]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy));  
						isval = 0;	/* Modified by Y. Shi */  
//			     			MD->EleETloss[i] =MD->EleET[i][0];
     						}
   		  			else
     						{
			     			isval = MD->EleIS[i]+stepsize*(((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)-MD->EleETloss[i]-MD->EleTF[i]);	/* Modified by Y. Shi */
						ret =  MD->EleTF[i];   
//     						MD->EleETloss[i] =MD->EleET[i][0];
  						}
	    				}
    				else if((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i])*stepsize) >= MD->EleISmax[i]))
    					{
//     					MD->EleETloss[i] =  MD->EleET[i][0];
      					isval = MD->EleISmax[i];	/* Modified by Y. Shi */
      					ret =  MD->EleTF[i]+(MD->EleIS[i]- MD->EleISmax[i])/stepsize + ((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)-(MD->EleETloss[i]+MD->EleTF[i]);

/*
					if (i == 1)
						{
						printf("\n CASE1CASE1CASE1CASE1CASE2CASE2CASE2CASE2CASE2!! ");
						printf("\n ret = %f", ret);
						printf("\n stepsize = %f", stepsize);
						printf("\n prep = %f", (1-fracSnow+MD->EleDew[i])*MD->ElePrep[i]);
						printf("\n dew = %f", MD->EleDew[i]);
						printf("\n vegfrac = %f", VegFrac);
						printf("\n MeltRateCanopy = %f", MeltRateCanopy);
						printf("\n loss = %f", MD->EleET[i][0]+MD->EleTF[i]);
						}
*/

    					}
	    			else if((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)-MD->EleET[i][0]-MD->EleTF[i])*stepsize) <=0))
    					{
	      				if((MD->EleET[i][0]>0)||(MD->EleTF[i]>0))
						{
       		 				MD->EleET[i][0] = (MD->EleETloss[i]/(MD->EleETloss[i]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy));
        					ret =(MD->EleTF[i]/(MD->EleETloss[i]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)); 
						}
	      				else
						{
						MD->EleET[i][0] = 0;
						ret = 0;
						}
//      				MD->EleETloss[i] =  MD->EleET[i][0];
      					isval = 0;     		/* Modified by Y. Shi */
    					}  
    				else
    					{
      					isval = MD->EleIS[i] + (((1-fracSnow)*MD->ElePrep[i]*VegFrac+MeltRateCanopy)-MD->EleETloss[i]-MD->EleTF[i])*stepsize;	/* Modified by Y. Shi */
//    		  			MD->EleETloss[i] = MD->EleET[i][0];
       					ret =  MD->EleTF[i]; 
    					}
				
				ET0 = MD->EleET[i][0]*(1000.0*Lv)/LEp;
				}
			else
				{
				FX = (MD->EleSM[i][0]-MD->Ele[i].ThetaW)/(MD->Ele[i].ThetaRef-MD->Ele[i].ThetaW);
				FX = (FX>1)?1:(FX<0?0:FX);
				ET2 = MD->pcCal.Et2*pow(FX,MD->fx_soil);
				ET1 = 0.0;
				ET0 = 0.0;
				isval = 0.0;
				ret = 0.0;
				}

			if (LEp>0)
				{
				LE = (ET0+ET1+ET2)*LEp;
				}
			else
				{
				LE = LEp;
				ET0 = 0;
				ET1 = 0;
				ET2 = 0;
				}

			MD->EleET[i][1] = ET1*LEp/(1000.0*Lv);
			MD->EleET[i][2] = ET2*LEp/(1000.0*Lv);

			MD->EleLE[i] = LE/(UNIT_C*60);

			Tsnew=(Fdown-B1-B2-LE-4*emiss*SIGMA*UNIT_C*pow(T+273.15,3)*273.15+3*emiss*SIGMA*UNIT_C*pow(T+273.15,4))/(A1+A2+4*emiss*SIGMA*UNIT_C*pow(T+273.15,3));
//			Tsnew=(Fdown-B1-B2-LE-4*SIGMA*UNIT_C*pow(T+273.15,3)*273.15+3*SIGMA*UNIT_C*pow(T+273.15,4))/(A1+A2+4*SIGMA*UNIT_C*pow(T+273.15,3));
			Terror = Tsnew-Tsold;
			Terror = Terror>0?Terror:(-Terror);
			Tsold = Tsnew;
			}

		for (j=0;j<4;j++)
			{
			MD->Tsoil[i][j] = soiltemp[j];
			}


//		printf("\nDelta = %f, Gamma = %f", Delta, Gamma);

		MD->Tsfc[i] = Tsnew;
		MD->EleG[i] =(A2*MD->Tsfc[i]+B2)/(UNIT_C*60);				/* Expanded by Y. Shi*/
		MD->EleH[i] =(A1*MD->Tsfc[i]+B1)/(UNIT_C*60);		/* Expanded by Y. Shi*/

		MD->EleIS[i] = isval;
    		MD->EleNetPrep[i] = (1-VegFrac)*(1-fracSnow)*MD->ElePrep[i]+ret+MeltRateGrnd;
//		if (i == 1) printf("\n Net precip = %f", MD->EleNetPrep[i]);
		MD->EleTF[i]=ret;
		MD->EleSnow[i] = MD->EleSnowGrnd[i] + MD->EleSnowCanopy[i];
		}


	/* Calculate evaporation from river */

	for(i=0; i<MD->NumRiv; i++)
  		{
		RivPrep = (MD->ElePrep[MD->Riv[i].RightEle-1]+MD->ElePrep[MD->Riv[i].LeftEle-1])/2;

		T = (Interpolation(&MD->TSD_Temp[MD->Ele[MD->Riv[i].RightEle-1].temp-1], t)+Interpolation(&MD->TSD_Temp[MD->Ele[MD->Riv[i].LeftEle-1].temp-1], t))/2;
//		printf("\n riv segment precip = %f, T = %f",RivPrep,T);
		/******************************************************************************************/
		/*			    Snow Accumulation/Melt Calculation				  */
		/******************************************************************************************/
    		fracSnow = T<Ts?1.0:(T>Tr?0:(Tr-T)/(Tr-Ts));
    		snowRate = fracSnow*RivPrep;

		/* EleSnowGrnd, EleSnowCanopy, EleISsnowmax, MeltRateGrnd,MeltRateCanopy are the average value prorated over the whole elemental area */
    		MD->EleSnow[i+MD->NumEle]=MD->EleSnow[i+MD->NumEle]+snowRate*stepsize;
    		MeltRateGrnd=(T>To?(T-To)*MF:0);		/* Note the units for MF. */
    		if(MD->EleSnow[i+MD->NumEle]>MeltRateGrnd*stepsize)
    			{
    			MD->EleSnow[i+MD->NumEle]=MD->EleSnow[i+MD->NumEle]-MeltRateGrnd*stepsize;
    			}
    		else
    			{
    			MeltRateGrnd=MD->EleSnow[i+MD->NumEle]/stepsize;
    			MD->EleSnow[i+MD->NumEle]=0;    	
    			}
		MD->EleEp[i+MD->NumEle] = (MD->EleEp[MD->Riv[i].RightEle-1] + MD->EleEp[MD->Riv[i].LeftEle-1])/2;
    		MD->EleNetPrep[i+MD->NumEle] = (1-fracSnow)*MD->ElePrep[i]+MeltRateGrnd;
  		}

	}

void excoef(void *DS, realtype T, realtype Tv, realtype Ts, realtype zlvl, realtype Vel, realtype P, realtype rl, int i)
	{

	/* Calculate surface exchange coefficients, extracted from Noah code */

	realtype czil;
	int ilech = 0;
	realtype Qs;
	realtype zilfc,dthv,zu,cxch,du2,wstar2,ustar,ustark,zt,zslu,zslt,rlogu,rlogt,rlmo;
	int itr;
	realtype zetalt,zetalu,zetau,zetat,xlu4,xlt4,xu4,xt4,xlu,xlt,xu,xt;
	realtype psmz,simm,pshz,simh,rlmn,rlma;
	realtype epsu2 = 0.0001*60*60,theta0=270.0,hpbl=1000.0,epsust=0.07*60;

  	Model_Data MD;
  	MD = (Model_Data)DS;

	/* Calculate C_h */

	MD->C_h[i] = 0.0001*60*UNIT_C;
	MD->C_m[i] = 0.0001*60*UNIT_C;
	czil = MD->Czil;
	zilfc = -czil*0.4*258.2;
	dthv = T-Ts;
	zu = rl;
	cxch = 0.001*UNIT_C*60/zlvl;
	du2 = pow(Vel,2.0)>epsu2*UNIT_C*UNIT_C?pow(Vel,2.0):epsu2*UNIT_C*UNIT_C;
	if (dthv == 0)
		{
		wstar2 = 0;
		}
	else
		{
		Qs = 1/theta0*GRAV*UNIT_C*UNIT_C*hpbl*MD->C_h[i]*dthv;
		Qs = Qs>0?Qs:(-Qs);	
		wstar2 = pow(1.2,2.0)*pow(Qs,2.0/3.0);
		}
	ustar = sqrt(MD->C_m[i]*sqrt(du2+wstar2))>epsust*UNIT_C?sqrt(MD->C_m[i]*sqrt(du2+wstar2)):epsust*UNIT_C;
	zt = zu*exp(zilfc*sqrt(ustar*zu/UNIT_C/60));
	zslu = zlvl+zu;
	zslt = zlvl+zt;
	rlogu = log(zslu/zu);
	rlogt = log(zslt/zt);
    
	rlmo = 0.4/theta0*GRAV*UNIT_C*UNIT_C*MD->C_h[i]*dthv/pow(ustar,3.0);

	for (itr=1;itr<6;itr++)
		{
		zetalt = zslt*rlmo>-5.0?zslt*rlmo:-5.0;
		rlmo = zetalt/zslt;
		zetalu = zslu*rlmo;
		zetau = zu*rlmo;
		zetat = zt*rlmo;

		if (ilech==0)
			{
			if (rlmo<0)
				{
				xlu4 = 1.0-16.0*zetalu;
				xlt4 = 1.0-16.0*zetalt;
				xu4 = 1.0-16.0*zetau;
				xt4 = 1.0-16.0*zetat;
				xlu = sqrt(sqrt(xlu4));
				xlt = sqrt(sqrt(xlt4));
				xu = sqrt(sqrt(xu4));
				xt = sqrt(sqrt(xt4));
				psmz = pspmu(xu);
				simm = pspmu(xlu)-psmz+rlogu;
				pshz = psphu(xt);
				simh = psphu(xlt)-pshz+rlogt;
				}
			else
				{
				zetalu = zetalu<1?zetalu:1.0;
				zetalt = zetalt<1?zetalt:1.0;
				psmz = pspms(zetau);
				simm = pspms(zetalu)-psmz+rlogu;
				pshz = psphs(zetat);
				simh = psphs(zetalt)-pshz+rlogt;
				}
			}
		else
			{
			if (rlmo<0)
				{
				psmz = pslmu(zetau);
				simm = pslmu(zetalu)-psmz+rlogu;
				pshz = pslhu(zetat);
				simh = pslhu(zetalt)-pshz+rlogt;
				}
			else
				{
				zetalu = zetalu<1?zetalu:1.0;
				zetalt = zetalt<1?zetalt:1.0;
				psmz = pslms(zetau);
				simm = pslms(zetalu)-psmz+rlogu;
				pshz=pslhs(zetat);
				simh=pslhs(zetalt)-pshz+rlogt;
				}
			}

		ustar = sqrt(MD->C_m[i]* sqrt(du2+wstar2))>epsust*UNIT_C?sqrt(MD->C_m[i]* sqrt(du2+wstar2)):epsust*UNIT_C;
		zt = zu*exp(zilfc*sqrt(ustar*zu/UNIT_C/60));

		zslt = zlvl+zt;
		rlogt = log(zslt/zt);
		ustark = ustar*0.4;
		MD->C_m[i] = ustark/simm>cxch?ustark/simm:cxch;
		MD->C_h[i] = ustark/simh>cxch?ustark/simh/simm:cxch;
		Qs = 1/theta0*GRAV*UNIT_C*UNIT_C*hpbl*MD->C_h[i]*dthv;
		Qs = Qs>0?Qs:(-Qs);	
		wstar2 = pow(1.2,2.0)*pow(Qs,2.0/3.0);
		rlmn = 0.4*GRAV*UNIT_C*UNIT_C/theta0*MD->C_h[i]*dthv/pow(ustar,3.0);
		rlma = rlmo*0.15+rlmn*(1.0-0.15);
		rlmo = rlma;
		}
	}
