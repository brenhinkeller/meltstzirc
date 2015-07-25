/******************************************************************************
 * FILE: meltsTzirc.c
 * DESCRIPTION:  
 *   Calculates ... using 
 *   ..
 * AUTHOR: C. Brenhin Keller
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <strings.h>
#include <math.h>
#include "arrays.h"
#include "runmelts.h"

double meltsM(double * const array){
	// Format: SiO2 TiO2 Al2O3 Fe2O3 Cr2O3 FeO MnO MgO NiO CoO CaO Na2O K2O P2O5 H2O
	for (int i=0; i<16; i++){
		if (i==0 || i==2 || i==5 || i==7 || i==10 || i==11 || i==12)  printf("%.1f, ", array[i]);
	}
	double Si=array[0]/(28.0844+15.9994*2);
	double Ti=array[1]/(47.867+15.9994*2);
	double Al=array[2]/(26.9815+15.9994*1.5);
	double Fe=array[3]/(55.845+15.9994*1.5) + array[5]/(55.845+15.9994);
	double Cr=array[4]/(51.9961+15.9994*1.5);
	double Mn=array[6]/(54.9380+15.9994);
	double Mg=array[7]/(24.3050+15.9994);
	double Ni=array[8]/(58.6934+15.9994);
	double Co=array[9]/(58.9332+15.9994);
	double Ca=array[10]/(40.078+15.9994);
	double Na=array[11]/(22.9898+15.9994/2);
	double K=array[12]/(39.0983+15.9994/2);
	double P=array[13]/(30.9738+15.9994*2.5);
	double TotalMoles = Si+Ti+Al+Fe+Cr+Mn+Mg+Ni+Co+Ca+Na+K+P;
	double M = (Na + K + 2*Ca) / (Al * Si) * TotalMoles;
	printf("\t\t\t\t%.2g",M);
	return M;
}

double meltsMmajors(double * const array){
	// Format: SiO2 TiO2 Al2O3 Fe2O3 FeO MgO CaO Na2O K2O H2O
	for (int i=0; i<9; i++){
		printf("%.1f, ", array[i]);
	}
	double Si=array[0]/(28.0844+15.9994*2);
	double Ti=array[1]/(47.867+15.9994*2);
	double Al=array[2]/(26.9815+15.9994*1.5);
	double Fe=array[3]/(55.845+15.9994*1.5) + array[4]/(55.845+15.9994);
	double Mg=array[5]/(24.3050+15.9994);
	double Ca=array[6]/(40.078+15.9994);
	double Na=array[7]/(22.9898+15.9994/2);
	double K=array[8]/(39.0983+15.9994/2);
	double H=array[9]/(1.008+15.9994/2);
	double TotalMoles = Si+Ti+Al+Fe+Mg+Ca+Na+K+H;
	double M = (Na + K + 2*Ca) / (Al * Si) * TotalMoles;
	printf("\t\t\t\t%.2g",M);
	return M;
}


double tzirc(const double M, const double Zr){
	double Tsat = 10108.0 / (log(496000.0/Zr) + 1.16*(M-1) + 1.48) - 273.15; // Temperature in Celcius
	if (Zr<=0){
		Tsat = NAN;
	}
	return Tsat;
}

double tzircZr(const double M, const double T){
	double Zrsat = 496000.0 / exp(10108.0/(T+273.15) - 1.16*(M-1) - 1.48);
	return Zrsat;
}

int main(int argc, char **argv){

	//Check input arguments
	if (argc != 2) {
		fprintf(stderr,"USAGE: %s <input_filename>\n", argv[0]);
		exit(1);
	}

	FILE *fp;
	char prefix[200], cmd_string[500];
	uint32_t datarows, datacolumns;
	uint32_t i, j, k;




	// Simulation parameters
	/**********************************************************/
	// Version to run MELTS in (MELTS or pMELTS)
	const char version[]="pMELTS";
	// Melts mode (isobaric, ptpath, etc)
	const char mode[]="isobaric";

	// fO2 buffer to use (None, FMQ, etc.)
	const char fo2Buffer[]="FMQ";
	// fO2 offset from buffer
	double fo2Delta=1;

	// Initial temperature (Celcius)
	double Ti=1700;
	//Initial Pressure (bar)
	double Pi=600;
	//Temperature step size in each simulation
	const int deltaT=-10;
	// Pressure step size;
	const int deltaP=0;
	
	// Variables that control size and location of the simulation
	/***********************************************************/	
	// Number of simulations to run
//	const int nsims=world_size*sims_per_task;
//	int n=world_rank-world_size;
//
	// Location of scratch directory (ideally local scratch for each node)
	// This location may vary on your system - contact your sysadmin if unsure
	const char scratchdir[]="/scratch/";

	// Variables that determine how much memory to allocate to imported results
	const int maxMinerals=40, maxSteps=1700/abs(deltaT), maxColumns=50;
	/***********************************************************/

	// Import 2-d source data array as a flat double array. Format:
	// SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, H2O, Zr;
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
	int row, P, T, mass, SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, CO2, H2O;
	double M, Tf, Tsat, Zrf, Zrsat, MZr;

	for (i=0;i<datarows;i++){


		//Configure working directory
		sprintf(prefix,"%sout%u/", scratchdir, i);
		sprintf(cmd_string,"mkdir -p %s", prefix);
		system(cmd_string);

		//Set water
//		data[i][15]=0.1;
		//Set CO2
//		data[i][14]=0.1;

		//Print current simulation
		printf("\nSimulation #%u\n",i);
		printf("SiO2, TiO2, Al203, Fe2O3, FeO, MgO, CaO, Na2O, K2O, H2O\n");
		for (j=0;j<datacolumns;j++){
			if (j==0 || j==1 || j==2 || j==3 || j==5 || j==7 || j==10 || j==11 || j==12 || j==15) printf("%.1f, ",data[i][j]);
		}
		printf("\n");

		//Run MELTS
		runmelts(prefix,data[i],version,mode,fo2Buffer,fo2Delta,"1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n0.99\n1\n10\n0\n4\n0\n","","!",Ti,Pi,deltaT,deltaP,0.005);

		// If simulation failed, clean up scratch directory and move on to next simulation
		sprintf(cmd_string,"%sPhase_main_tbl.txt", prefix);
		if ((fp = fopen(cmd_string, "r")) == NULL) {
			fprintf(stderr, "%u: MELTS equilibration failed to produce output.\n", i);
//			sprintf(cmd_string,"rm -r %s", prefix);
//			system(cmd_string);
			continue;
		}

		// Import results, if they exist. Format:
		// Pressure Temperature mass S H V Cp viscosity SiO2 TiO2 Al2O3 Fe2O3 Cr2O3 FeO MnO MgO NiO CoO CaO Na2O K2O P2O5 H2O
		importmelts(prefix, melts, rawMatrix, meltsrows, meltscolumns, names, elements, &minerals);
		if (minerals<1 | strcmp(names[0],"liquid_0")!=0) {
			fprintf(stderr, "%u: MELTS equilibration failed to calculate liquid composition.\n", i);
//			sprintf(cmd_string,"rm -r %s", prefix);
//			system(cmd_string);
			continue;
		}


		// Find the columns containing useful elements
		for(int col=0; col<meltscolumns[0]; col++){
			if (strcmp(elements[0][col], "Pressure")==0) P=col;
			else if (strcmp(elements[0][col], "Temperature")==0) T=col;
			else if (strcmp(elements[0][col], "mass")==0) mass=col;
			else if (strcmp(elements[0][col], "SiO2")==0) SiO2=col;
			else if (strcmp(elements[0][col], "TiO2")==0) TiO2=col;
			else if (strcmp(elements[0][col], "Al2O3")==0) Al2O3=col;
			else if (strcmp(elements[0][col], "Fe2O3")==0) Fe2O3=col;
			else if (strcmp(elements[0][col], "Cr2O3")==0) Cr2O3=col;
			else if (strcmp(elements[0][col], "FeO")==0) FeO=col;
			else if (strcmp(elements[0][col], "MnO")==0) MnO=col;
			else if (strcmp(elements[0][col], "MgO")==0) MgO=col;
			else if (strcmp(elements[0][col], "NiO")==0) NiO=col;
			else if (strcmp(elements[0][col], "CoO")==0) CoO=col;
			else if (strcmp(elements[0][col], "CaO")==0) CaO=col;
			else if (strcmp(elements[0][col], "Na2O")==0) Na2O=col;
			else if (strcmp(elements[0][col], "K2O")==0) K2O=col;
			else if (strcmp(elements[0][col], "P2O5")==0) P2O5=col;
			else if (strcmp(elements[0][col], "CO2")==0) CO2=col;
			else if (strcmp(elements[0][col], "H2O")==0) H2O=col;
		}
		
		// Calculate saturation temperature and minimum necessary zirconium content
		
		Tsat=0;	
		for(row=1; row<(meltsrows[0]-1); row++){
			//Calculate melt M and [Zr]
			M = meltsM(&melts[0][row][SiO2]);
			Zrf = data[i][datacolumns-1]*100/(melts[0][row][mass] + 0.01*(100-melts[0][row][mass])); // Assuming bulk Kd=0.1
			printf("\t%.0f\t%.0f\n", tzirc(M, Zrf), melts[0][row][T]);

			// Check if we've cooled below the saturation temperature yet
		       	if (Tsat==0 && tzirc(M, Zrf) > melts[0][row][T]){
				Tsat = tzirc(M, Zrf);
				printf("Saturation reached\n");
			}
 			// Stop when we get to maximum SiO2
			if (melts[0][row-1][SiO2]>(melts[0][row][SiO2])+0.01){
				break;
			}

			// Or when remaining melt falls below 5%
			if (melts[0][row][mass]<5){
				break;
			}
		}
	
		M = meltsM(&melts[0][row][SiO2]);
		Zrf = data[i][datacolumns-1]*100/(melts[0][row][mass] + 0.01*(100-melts[0][row][mass])); // Final zirconium content, assuming bulk kd=0.1
		printf("\t%.0f\t%.0f\n", tzirc(M, Zrf), melts[0][row][T]);

		Tf = melts[0][row][T];
		Zrsat = tzircZr(M, Tf);
		if (Tsat==0){
			Tsat = tzirc(M, Zrf);
		}

		printf("Final melt percentage: %g, Water Content: %g, M value: %g\n", melts[0][row][mass], melts[0][row][H2O], M);
		printf("Final Zr: %g, Zr at saturation: %g, Saturation temperature: %g, Final T: %g\n", Zrf, Zrsat, Tsat, Tf);

		// Determine how much zircon is saturated
		if (Zrf>Zrsat){
			MZr=melts[0][row][mass]/100*(Zrf-Zrsat);
		} else {
			MZr=0;
		}
		printf("Mass of zircon saturated: %g\n", MZr);

	}





	for (i=0;i<datarows;i++){
		for (j=0;j<datacolumns;j++){
			printf("%g,",data[i][j]);
		}
		printf("\n");
	}

return 0;
}

