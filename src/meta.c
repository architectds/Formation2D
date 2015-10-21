/*********************************************************************************
 * File        : meta.c               		                                 *
 * Function    : create meta file for PIHM 2.0                     	         *
 * Version     : Feb, 2011 (2.0)                                                 *
 * Developer of PIHM2.0:	Mukesh Kumar (muk139@psu.edu)		         * 
 * Developer of PIHM1.0:	Yizhong Qu   (quyizhong@gmail.com)	         * 
**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "pihm.h"  


void create_meta(char *filename, Model_Data DS, Control_Data *CS)
	{
  	char *meta_name;
  	FILE *meta_file;	/* Pointer to meta file */
	time_t t, t1, t2;
  	int precip, temp, RH, wind, radiation, pressure, topography, soil, BRD, LC;
	time(&t);
	struct tm *timeinfo;
	timeinfo = (struct tm *)malloc(sizeof(struct tm));
	timeinfo = localtime(&t);

  	printf("\nStart creating meta file ... \n");
  
  	/*========== open meta file ==========*/ 

	meta_name = (char *)malloc((strlen(filename)+40)*sizeof(char));
	strcpy(meta_name,"output/");
  	strcat(meta_name, filename);
	strcat(meta_name,".xml");
	meta_file=fopen(meta_name,"w");

	fprintf(meta_file,"%s","<?xml version=\"2.0\"?>");
	fprintf(meta_file,"%s","\n<PIHM xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\">");
	fprintf(meta_file,"%s","\n\t<File_Information>");
	fprintf(meta_file,"\n\t\t<File_Name>%s</File_Name>",meta_name);
	fprintf(meta_file,"\n\t\t<Type>MODEL</Type>");
	fprintf(meta_file,"\n\t\t<Version>2.2</Version>");
	fprintf(meta_file,"\n\t\t<Created>%4.4d-%2.2d-%2.2d %2.2d:%2.2d</Created>",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour,timeinfo->tm_min);
	fprintf(meta_file,"\n\t</File_Information>");
	fprintf(meta_file,"\n\t<Forcing>");

	while (1)
		{
		printf("\n Which precipitation data is being used? (1[RTHnet]/2[NARR]/3[NLDAS])");
//		scanf("%d",&precip);
		precip = 1;
		if (precip == 1)
			{
			fprintf(meta_file,"\n\t\t<Precip>RTHnet</Precip>");
			break;
			}
		else if (precip == 2)
			{
			fprintf(meta_file,"\n\t\t<Precip>NARR</Precip>");
			break;
			}
		else if (precip == 3)
			{
			fprintf(meta_file,"\n\t\t<Precip>NLDAS</Precip>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}

	while (1)
		{
		printf("\n Which temperature data is being used? (1[RTHnet]/2[NARR]/3[NLDAS])");
//		scanf("%d",&temp);
		temp = 1;
		if (temp == 1)
			{
			fprintf(meta_file,"\n\t\t<Temp>RTHnet</Temp>");
			break;
			}
		else if (temp == 2)
			{
			fprintf(meta_file,"\n\t\t<Temp>NARR</Temp>");
			break;
			}
		else if (temp == 3)
			{
			fprintf(meta_file,"\n\t\t<Temp>NLDAS</Temp>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	while (1)
		{
		printf("\n Which RH data is being used? (1[RTHnet]/2[NARR]/3[NLDAS])");
//		scanf("%d",&RH);
		RH = 1;
		if (RH == 1)
			{
			fprintf(meta_file,"\n\t\t<RH>RTHnet</RH>");
			break;
			}
		else if (RH == 2)
			{
			fprintf(meta_file,"\n\t\t<RH>NARR</RH>");
			break;
			}
		else if (RH == 3)
			{
			fprintf(meta_file,"\n\t\t<RH>NLDAS</RH>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	while (1)
		{
		printf("\n Which wind speed data is being used? (1[NARR]/2[NLDAS]/3[RTHnet])");
//		scanf("%d",&wind);
		wind = 1;
		if (wind == 1)
			{
			fprintf(meta_file,"\n\t\t<Wind>NARR</Wind>");
			break;
			}
		else if (wind == 2)
			{
			fprintf(meta_file,"\n\t\t<Wind>NLDAS</Wind>");
			break;
			}
		else if (wind == 3)
			{
			fprintf(meta_file,"\n\t\t<Wind>RTHnet</Wind>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}

	while (1)
		{
		printf("\n Which radiation data is being used? (1[SURFRAD]/2[NARR]/3[RTHnet])");
//		scanf("%d",&radiation);
		radiation = 1;
		if (radiation == 1)
			{
			fprintf(meta_file,"\n\t\t<Radiation>SURFRAD</Radiation>");
			break;
			}
		else if (radiation == 2)
			{
			fprintf(meta_file,"\n\t\t<Radiation>NARR</Radiation>");
			break;
			}
		else if (radiation == 3)
			{
			fprintf(meta_file,"\n\t\t<Radiation>RTHnet</Radiation>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	while (1)
		{
		printf("\n Which pressure data is being used? (1[NARR]/2[NLDAS]/3[RTHnet])");
//		scanf("%d",&pressure);
		pressure = 1;
		if (pressure == 1)
			{
			fprintf(meta_file,"\n\t\t<Pressure>NARR</Pressure>");
			break;
			}
		else if (pressure == 2)
			{
			fprintf(meta_file,"\n\t\t<Pressure>NLDAS</Pressure>");
			break;
			}
		else if (pressure == 3)
			{
			fprintf(meta_file,"\n\t\t<Pressure>RTHnet</Pressure>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	fprintf(meta_file,"\n\t</Forcing>");
	fprintf(meta_file,"\n\t<GIS_FILE>");

	while (1)
		{
		printf("\n Which topography data is being used? (1[1-m CZO LiDAR]/2[3-m GPS])");
//		scanf("%d",&topography);
		topography = 1;
		if (topography == 1)
			{
			fprintf(meta_file,"\n\t\t<Topography>1-m CZO LiDAR project</Topography>");
			break;
			}
		else if (topography == 2)
			{
			fprintf(meta_file,"\n\t\t<Topography>3-m GPS DEM</Topography>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	while (1)
		{
		printf("\n Which soil data is being used? (1[Lin survey]/2[SSURGO])");
//		scanf("%d",&soil);
		soil = 1;
		if (soil == 1)
			{
			fprintf(meta_file,"\n\t\t<Soil>Lin survey</Soil>");
			break;
			}
		else if (soil == 2)
			{
			fprintf(meta_file,"\n\t\t<Soil>SSURGO</Soil>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	while (1)
		{
		printf("\n Which bedrock depth data is being used? (1[Lin survey]/2[2-m uniform])");
//		scanf("%d",&BRD);
		BRD = 1;
		if (BRD == 1)
			{
			fprintf(meta_file,"\n\t\t<Bed_Rock_Depth>Lin survey</Bed_Rock_Depth>");
			break;
			}
		else if (BRD == 2)
			{
			fprintf(meta_file,"\n\t\t<Bed_Rock_Depth>2-m uniform</Bed_Rock_Depth>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}
	while (1)
		{
		printf("\n Which land cover data is being used? (1[Eissenstat survey]/2[NLCD])");
//		scanf("%d",&LC);
		LC = 1;
		if (LC == 1)
			{
			fprintf(meta_file,"\n\t\t<Land_Cover>Eissenstat survey</Land_Cover>");
			break;
			}
		else if (LC == 2)
			{
			fprintf(meta_file,"\n\t\t<Land_Cover>NLCD</Land_Cover>");
			break;
			}
		else
			{
			printf("\n ERROR!\n");
			}
		}

	fprintf(meta_file,"\n\t</GIS_FILE>");
	fprintf(meta_file,"\n\t<Miscellaneous>");
	fprintf(meta_file,"\n\t\t<Ele>%3.3d</Ele>",DS->NumEle);
	fprintf(meta_file,"\n\t\t<Channel>%2.2d</Channel>",DS->NumRiv);
//	free(timeinfo);
	t1 = (time_t)CS->StartTime*60;
	timeinfo = gmtime(&t1);
	fprintf(meta_file,"\n\t\t<Start>%4.4d-%2.2d-%2.2d %2.2d:%2.2d</Start>",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour,timeinfo->tm_min);
	t2 = (time_t)CS->EndTime*60;
//	free(timeinfo);
	timeinfo = gmtime(&t2);
	fprintf(meta_file,"\n\t\t<End>%4.4d-%2.2d-%2.2d %2.2d:%2.2d</End>",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour,timeinfo->tm_min);
	fprintf(meta_file,"\n\t</Miscellaneous>");
	fprintf(meta_file,"\n</PIHM>");
	fclose(meta_file);
//	free(timeinfo);
        free(meta_name);
	}

