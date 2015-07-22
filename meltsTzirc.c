/******************************************************************************
 * FILE: meltsTzirc.c
 * DESCRIPTION:  
 *   Calculates ... using 
 *   ..
 * AUTHOR: C. Brenhin Keller
 * LAST REVISED: 10/08/13
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <strings.h>
#include <math.h>
#include "arrays.h"
#include "runmelts.h"


// 'SiO2';'TiO2';'Al2O3';'Fe2O3';'Cr2O3';'FeO';  'MnO';  'MgO';   'NiO';  'CoO';  'CaO';  'Na2O'; 'K2O'; 'P2O5';
// SiO2, TiO2, Al2O3, Cr2O3, FeOT, MnO, MgO, Ni, Co, CaO, Na2O, K2O, P2O5


// SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, H2O;
	

int main(int argc, char **argv){

	//Check input arguments
	if (argc != 2) {
		fprintf(stderr,"USAGE: %s <input_filename>\n", argv[0]);
		exit(1);
	}

	uint32_t datarows, datacolumns;
//	uint32_t i, j;

	char prefix[200], cmd_string[500];
	FILE *fp;
	uint32_t i, j, k;





	// Variables that control size and location of the simulation
	/***********************************************************/	
	// Number of simulations to run , and 
//	const int nsims=world_size*sims_per_task, 
//	int n=world_rank-world_size;
	// Location of scratch directory (ideally local scratch for each node)
	// This location may vary on your system - contact your sysadmin if unsure
	const char scratchdir[]="/scratch/";

	// Simulation parameters
	double Pi=9001, fo2Delta=0;
	//Temperature step size in each simulation
	const int deltaT=-10;
	/***********************************************************/	
	const int maxMinerals=40, maxSteps=1700/abs(deltaT), maxColumns=50;


	// import 2-d source data array as a flat double array
	double** const restrict data = csvparse(argv[1],',', &datarows, &datacolumns);

	// Malloc space for the imported melts array
	double **rawMatrix=mallocDoubleArray(maxMinerals*maxSteps,maxColumns);
	double ***melts=malloc(maxMinerals*sizeof(double**));
	char **names=malloc(maxMinerals*sizeof(char*));
	char ***elements=malloc(maxMinerals*sizeof(char**));
	int *meltsrows=malloc(maxMinerals*sizeof(int)), *meltscolumns=malloc(maxMinerals*sizeof(int));
	for (i=0; i<maxMinerals; i++){
		names[i]=malloc(30*sizeof(char));
		elements[i]=malloc(maxColumns*sizeof(char*));
		for (k=0; k<maxColumns; k++){
			elements[i][k]=malloc(30*sizeof(char));
		}
	}
	int minerals;

	//  Variables for finding saturation temperature
	int SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, CO2, H2O;




	for (i=0;i<datarows;i++){
		//Print current simulation
		printf("Simulation #%u\n",i);
		for (j=0;j<datacolumns;j++){
			printf("%g,",data[i][j]);
		}
		printf("\n");

		//Configure working directory
		sprintf(prefix,"%sout%u/", scratchdir, i);
		sprintf(cmd_string,"mkdir -p %s", prefix);
		system(cmd_string);

		//Run MELTS to equilibrate at 5% melt
		runmeltsNoCO2(prefix,data[i],"pMELTS","isobaric","FMQ",fo2Delta,"1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n0.99\n1\n10\n0\n4\n0\n","","!",1700,Pi,deltaT,0,0.005);

		// If simulation failed, clean up scratch directory and move on to next simulation
		sprintf(cmd_string,"%sPhase_main_tbl.txt", prefix);
		if ((fp = fopen(cmd_string, "r")) == NULL) {
			fprintf(stderr, "%u: MELTS equilibration failed to produce output.\n", i);
//			sprintf(cmd_string,"rm -r %s", prefix);
//			system(cmd_string);
			continue;
		}


		// Import results, if they exist
		importmelts(prefix, melts, rawMatrix, meltsrows, meltscolumns, names, elements, &minerals);
		if (minerals<1 | strcmp(names[0],"liquid_0")!=0) {
			fprintf(stderr, "%u: MELTS equilibration failed to calculate liquid composition.\n", i);
//			sprintf(cmd_string,"rm -r %s", prefix);
//			system(cmd_string);
			continue;
		}


		// Find the columns containing useful elements
		for(i=0; i<columns[0]; i++){
			if (strcmp(elements[0][i], "SiO2")==0) SiO2=i;
			else if (strcmp(elements[0][i], "TiO2")==0) TiO2=i;
			else if (strcmp(elements[0][i], "Al2O3")==0) Al2O3=i;
			else if (strcmp(elements[0][i], "Fe2O3")==0) Fe2O3=i;
			else if (strcmp(elements[0][i], "Cr2O3")==0) Cr2O3=i;
			else if (strcmp(elements[0][i], "FeO")==0) FeO=i;
			else if (strcmp(elements[0][i], "MnO")==0) MnO=i;
			else if (strcmp(elements[0][i], "MgO")==0) MgO=i;
			else if (strcmp(elements[0][i], "NiO")==0) NiO=i;
			else if (strcmp(elements[0][i], "CoO")==0) CoO=i;
			else if (strcmp(elements[0][i], "CaO")==0) CaO=i;
			else if (strcmp(elements[0][i], "Na2O")==0) Na2O=i;
			else if (strcmp(elements[0][i], "K2O")==0) K2O=i;
			else if (strcmp(elements[0][i], "P2O5")==0) P2O5=i;
			else if (strcmp(elements[0][i], "CO2")==0) CO2=i;
			else if (strcmp(elements[0][i], "H2O")==0) H2O=i;
		}



	}





	for (i=0;i<datarows;i++){
		for (j=0;j<datacolumns;j++){
			printf("%g,",data[i][j]);
		}
		printf("\n");
	}

return 0;
}

