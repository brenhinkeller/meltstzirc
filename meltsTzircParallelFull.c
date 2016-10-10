/******************************************************************************
 * FILE: meltsTzircFullParallel.c
 * DESCRIPTION:  
 *   As meltsTzircParallel.c, but prints output differently, returning one line
 * for each cooling step of each magma.
 * AUTHOR: C. Brenhin Keller
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <mpi.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "arrays.h"
#include "runmelts.h"
#include "GERM.h"

double meltsM(double * const array){
	// Format: SiO2 TiO2 Al2O3 Fe2O3 Cr2O3 FeO MnO MgO NiO CoO CaO Na2O K2O P2O5
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
	return M;
}

double meltsMmajors(double * const array){
	// Format: SiO2 TiO2 Al2O3 Fe2O3 FeO MgO CaO Na2O K2O 
	double Si=array[0]/(28.0844+15.9994*2);
	double Ti=array[1]/(47.867+15.9994*2);
	double Al=array[2]/(26.9815+15.9994*1.5);
	double Fe=array[3]/(55.845+15.9994*1.5) + array[4]/(55.845+15.9994);
	double Mg=array[5]/(24.3050+15.9994);
	double Ca=array[6]/(40.078+15.9994);
	double Na=array[7]/(22.9898+15.9994/2);
	double K=array[8]/(39.0983+15.9994/2);
	double TotalMoles = Si+Ti+Al+Fe+Mg+Ca+Na+K;
	double M = (Na + K + 2*Ca) / (Al * Si) * TotalMoles;
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
	uint32_t datarows, datacolumns;
	uint32_t i, j, k;
	int world_size, world_rank, rc;


	//Check input arguments
	if (argc != 2) {
		fprintf(stderr,"USAGE: %s <input_filename>\n", argv[0]);
		exit(1);
	}

	// Start MPI
	rc = MPI_Init(&argc,&argv); 
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc);
	}

	// Get world size (number of MPI processes) and world rank (# of this process)
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);


	if (world_rank==0){
		// Declare variables used only on the root node
		int buf[world_size-1], nextReady;
		MPI_Request reqs[world_size-1];
		MPI_Status stats[world_size-1];

		// Print format of output
		printf("Kv\tT\tF\tM\tSiO2\tZr\tZrsat\tMZrIncr\tTSat\n");

		// Import 2-d source data array as a flat double array. Format:
		// SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, CO2, H2O, Zr, Kv;
		double** const data = csvparse(argv[1],',', &datarows, &datacolumns);

		// Listen for task requests from the worker nodes
		for (i=1; i<world_size; i++){
			//        *buf, count, datatype, dest, tag, comm, *request
			MPI_Irecv(&buf[i-1], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[i-1]);
		}

		// Once any worker asks for a new task, send next task to that worker and keep listening
		for (i=0; i<datarows; i++){
			MPI_Waitany(world_size-1, reqs, &nextReady, stats);
			//       *buf, count, datatype, dest, tag, comm
			MPI_Send(data[i], 18, MPI_DOUBLE, nextReady+1, 1, MPI_COMM_WORLD);
			//        *buf, count, datatype, source, tag, comm, *request
			MPI_Irecv(&buf[nextReady], 1, MPI_INT, nextReady+1, 0, MPI_COMM_WORLD, &reqs[nextReady]);
		}

		// Wait for all workers to complete, then send the stop signal
		MPI_Waitall(world_size-1, reqs, stats);	
		double stop[18] = {-1};
		for (i=1; i<world_size; i++){
			MPI_Send(&stop, 18, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);	
		}
	}

	else {
		// Declare variables used only on the worker nodes
		MPI_Request sReq;
		MPI_Status sStat;
		double ic[18], Kd, iKd;
		FILE *fp;
		char* prefix = malloc(500*sizeof(char));
		char* cmd_string = malloc(1000*sizeof(char));

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
		double Pi=6000;
		//Temperature step size in each simulation
		const int deltaT=-10;
		// Pressure step size;
		const int deltaP=0;

		// Stop simulations at a given percent melt
		const double minPercentMelt=0.1;

		// Variables that control size and location of the simulation
		/***********************************************************/	
		// Location of scratch directory (ideally local scratch for each node)
		// This location may vary on your system - contact your sysadmin if unsure
//		const char scratchdir[]="/scratch/gpfs/cbkeller/";
		const char scratchdir[]="/scratch/";

		// Variables that determine how much memory to allocate to imported results
		const int maxMinerals=100, maxSteps=1700/abs(deltaT), maxColumns=50;
		/***********************************************************/


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
		int row, col, P, T, mass, SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, CO2, H2O;
		int fspCaO, fspNa2O, fspK2O, oxideTiO2, oxideFe2O3, oxideFeO, oxideMnO;

		double M, Tf, Tsat, Zrf, Zrsat, MZr, MZrLast, MZrIncrement;
		double AnKd, AbKd, OrKd, IlmKd, MtKd;


		while (1) {
			// Ask root node for new task
			//       *buf, count, datatype, dest, tag, comm, *request
			MPI_Isend(&world_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &sReq);
			//       *buf, count, datatype, source, tag, comm, *status
			MPI_Recv(&ic, 18, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &sStat);

			// Exit loop if stop signal recieved
			if (ic[0]<0) break;

			//Configure working directory
			sprintf(prefix,"%sout%i_%.0f/", scratchdir, world_rank, ic[17]);
			sprintf(cmd_string,"mkdir -p %s", prefix);
			system(cmd_string);

//			//Set water
			ic[15]=3.0;
//			//Set CO2
//			ic[14]=0.1;
			
			//Run MELTS
			runmelts(prefix,ic,version,"isobaric",fo2Buffer,fo2Delta,"1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n0.99\n1\n10\n0\n4\n0\n","","!",Ti,Pi,deltaT,deltaP,0.005);

			// If simulation failed, clean up scratch directory and move on to next simulation
			sprintf(cmd_string,"%sPhase_main_tbl.txt", prefix);
			if ((fp = fopen(cmd_string, "r")) == NULL) {
				fprintf(stderr, "%s : MELTS equilibration failed to produce output.\n", prefix);
				sprintf(cmd_string,"rm -r %s", prefix);
				system(cmd_string);
				continue;
			}

			// Import results, if they exist. Format:
			// Pressure Temperature mass S H V Cp viscosity SiO2 TiO2 Al2O3 Fe2O3 Cr2O3 FeO MnO MgO NiO CoO CaO Na2O K2O P2O5 H2O
			minerals=maxMinerals;
			importmelts(maxSteps, maxColumns, prefix, melts, rawMatrix, meltsrows, meltscolumns, names, elements, &minerals);
			if (minerals<1 | strcmp(names[0],"liquid_0")!=0) {
				fprintf(stderr, "%s : MELTS equilibration failed to calculate liquid composition.\n", prefix);
				sprintf(cmd_string,"rm -r %s", prefix);
				system(cmd_string);
				continue;
			}
			// Can delete temp files after we've read them
			sprintf(cmd_string,"rm -r %s", prefix);
			system(cmd_string);


			// Find the columns containing useful elements for the melt
			for(col=0; col<meltscolumns[0]; col++){
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

			// Find the columns containing useful elements for other minerals
			for (i=1; i<minerals; i++){
				if (strncasecmp(names[i],"feldspar",8)==0){
					for(col=0; col<meltscolumns[i]; col++){
						if (strcmp(elements[i][col], "CaO")==0) fspCaO=col;
						else if (strcmp(elements[i][col], "Na2O")==0) fspNa2O=col;
						else if (strcmp(elements[i][col], "K2O")==0) fspK2O=col;
					}

				} else if (strncasecmp(names[i],"rhm_oxide",9)==0){
					for(col=0; col<meltscolumns[i]; col++){
						if (strcmp(elements[i][col], "TiO2")==0) oxideTiO2=col;
						else if (strcmp(elements[i][col], "Fe2O3")==0) oxideFe2O3=col;
						else if (strcmp(elements[i][col], "FeO")==0) oxideFeO=col;
						else if (strcmp(elements[i][col], "MnO")==0) oxideMnO=col;
					}

				}
			}




			// Calculate saturation temperature and minimum necessary zirconium content	
			MZrLast=0;
			for(row=0; row<(meltsrows[0]-1); row++){

				// Calculate bulk zircon partition coefficient at present step
				Kd = 0;
				for (i=1; i<minerals; i++){
					// See what minerals might be crystallizing at this temperature
					// so we can find their GERM partition coefficients
					for (j=0; j<meltsrows[i]; j++){
						if (fabs(melts[0][row][T]-melts[i][j][T]) < 0.01){
							if (strncasecmp(names[i],"feldspar",8)==0){
								AnKd = getGERMKd("AnKdorthite","Zr",melts[0][row][SiO2]);
								AbKd = getGERMKd("Albite","Zr",melts[0][row][SiO2]);
								OrKd = getGERMKd("Orthoclase","Zr",melts[0][row][SiO2]);
								if (isnan(AnKd)) AnKd=0;
								if (isnan(OrKd)) OrKd=0;
								if (isnan(AbKd)) AbKd = (AnKd + OrKd)/2;

								iKd = (220.1298+56.18)/56.18*melts[i][j][fspCaO]/100 * AnKd\
								      +(228.2335+30.99)/30.99*melts[i][j][fspNa2O]/100 * AbKd\
								      +(228.2335+47.1)/47.1*melts[i][j][fspK2O]/100 * OrKd;

							} else if (strncasecmp(names[i],"rhm_oxide",9)==0){
								IlmKd = getGERMKd("Ilmenite","Zr",melts[0][row][SiO2]);
								MtKd = getGERMKd("Magnetite","Zr",melts[0][row][SiO2]);
								if (isnan(IlmKd)) IlmKd = 0;
								if (isnan(MtKd)) MtKd = 0;

								iKd = (melts[i][j][oxideTiO2]+melts[i][j][oxideMnO]+(melts[i][j][oxideTiO2]\
									*(71.8444/79.8768)-melts[i][j][oxideMnO]*(71.8444/70.9374)))/100 * AnKd\
								      + (1 - (melts[i][j][oxideTiO2]+melts[i][j][oxideMnO]+(melts[i][j][oxideTiO2]\
									*(71.8444/79.8768)-melts[i][j][oxideMnO]*(71.8444/70.9374)))/100) * MtKd;
							} else {
								iKd = getGERMKd(names[i],"Zr",melts[0][row][SiO2]);
							}

							if (isnan(iKd)){iKd = 0;}
//							printf("Mineral: %s, Kd: %g\n",names[i],iKd);
							Kd += iKd * melts[i][j][mass];
						}
					}		
				}
				Kd = Kd / (100 - melts[0][row][mass]);
//				printf("Bulk Kd: %g\n",Kd);


				// Calculate zircon saturation state at present step
				M = meltsM(&melts[0][row][SiO2]);
				Zrf = ic[16]*100/(melts[0][row][mass] + Kd*(100-melts[0][row][mass])); // Assuming bulk Kd=0.1
				Tsat = tzirc(M, Zrf);
				Tf = melts[0][row][T];
				Zrsat = tzircZr(M, Tf);

				if (Zrf>Zrsat){
					MZr=melts[0][row][mass]/100*(Zrf-Zrsat);
				} else {
					MZr=0;
				}
				MZrIncrement=MZr-MZrLast;
				if (MZr>MZrLast){
					MZrLast=MZr;
				}

//				// Stop when we get to maximum SiO2
//				if (row>0){
//					if (melts[0][row-1][SiO2]>(melts[0][row][SiO2])+0.01){
//						break;
//					}
//				}

				// Or when remaining melt falls below minimum percent
				if (melts[0][row][mass]<minPercentMelt){
					break;
				}

				// Print results. Format:
				// Kv, T, F, M, SiO2, Zr, Zrsat, MZr, TSat
				printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",ic[17], Tf, melts[0][row][mass], M, melts[0][row][SiO2], Zrf, Zrsat, MZrIncrement, Tsat);

			}

		}
	}
	MPI_Finalize();
	return 0;
}

