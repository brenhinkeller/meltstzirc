#include <stdlib.h>
#include <math.h>


// Interpolate the value of a double array at a single non-interger index
double interp1i(const double* array, double i){
	uint32_t f = (uint32_t) floor(i);
	
	return array[f+1]*(i-f) + array[f]*(1 - (i-f));	
}



// Make a 1d array of points spaced by Step between Lower and Upper
double* array(double lower, double step, double upper){
	if (step<=0){
		perror("Step size too small");
	} 
	else {
		int i=0, numsteps=(upper-lower)/step+1;
		double *someArray=malloc(numsteps*sizeof(double));
		while (i<numsteps) {
			someArray[i]=lower+i*step;
			i++;
		}
		return someArray;
	}
	return NULL;
}

int* arrayInt(int lower, int step, int upper){
	if (step<=0){
		perror("Step size too small");
	} 
	else {
		int i=0, numsteps=(upper-lower)/step+1;
		int  *someArray=malloc(numsteps*sizeof(int));
		while (i<numsteps) {
			someArray[i]=lower+i*step;
			i++;
		}
		return someArray;
	}
	return NULL;
}

/* Make a 1-d array of n points linearly spaced between Lower and Upper */
double* linspace(double lower, double upper, double n){
	int i=0;
	double *someArray=malloc(n*sizeof(double));
	while (i<n){
		someArray[i]=(upper-lower)*i/(n-1)+lower;
		i++;
	}
	return someArray;
}

int* linspaceInt(int lower, int upper, int n){
	int i=0, *someArray=malloc(n*sizeof(int));
	while (i<n){
		someArray[i]=(upper-lower)*i/(n-1)+lower;
		i++;
	}
	return someArray;
}

/* Find the minimum of a 1-d array */
double min(const double a, const double b){
	if (b<a) return b;
	return a;
}

int minInt(const int a, const int b){
	if (b<a) return b;
	return a;
}

double minArray(const double * const someArray,const uint32_t length){
	uint32_t i;
	double currentMin=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]<currentMin){currentMin=someArray[i];}
	}
	return currentMin;
}

int minArrayInt(const int * const someArray,const uint32_t length){
	uint32_t i, currentMin=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]<currentMin){currentMin=someArray[i];}
	}
	return currentMin;
}

/* Find the maximum of a 1-d array */
double max(const double a, const double b){
	if (b>a) return b;
	return a;
}

int maxInt(const int a, const int b){
	if (b>a) return b;
	return a;
}

double maxArray(const double * const someArray, const uint32_t length){
	uint32_t i;
	double currentMax=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]>currentMax){currentMax=someArray[i];}
	}
	return currentMax;
}

int maxArrayInt(const int * const someArray, const uint32_t length){
	uint32_t i, currentMax=someArray[0];
	for (i=1; i<length; i++){
		if (someArray[i]>currentMax){currentMax=someArray[i];}
	}
	return currentMax;
}

void copyArray(const double * const A, const uint32_t n, double * const restrict B){
	for (uint32_t i=0; i<n; i++){
		B[i]=A[i];
	}
}


void copyArrayInt(const int * const A, const uint32_t n, int * const restrict B){
	for (uint32_t i=0; i<n; i++){
		B[i]=A[i];
	}
}

void copyArrayUint(const uint32_t * const A, const uint32_t n, uint32_t * const restrict B){
	for (uint32_t i=0; i<n; i++){
		B[i]=A[i];
	}
}


/* Malloc a 2D array of doubles of dimensions Rows x Columns */
double **mallocDoubleArray(const uint32_t rows, const uint32_t columns) {
	uint32_t i; double **someArray;
	someArray = (double **) malloc(rows*sizeof(double *));
	for (i = 0; i < rows; i++){
		someArray[i] = (double *) malloc(columns*sizeof(double));
	}
	return someArray;
} 

/* Free a 2D array of doubles of length rows */
void freeDoubleArray(double **someArray, const uint32_t rows) {
	uint32_t i;
	for (i=0; i<rows; i++) {free(someArray[i]);}
	free(someArray);
}


/* Malloc a 2D array of ints of dimensions Rows x Columns */
int **mallocIntArray(const uint32_t rows, const uint32_t columns) {
	uint32_t i; int **someArray;
	someArray = (int **) malloc(rows*sizeof(int *));
	for (i = 0; i < rows; i++){
		someArray[i] = (int *) malloc(columns*sizeof(int));
	}
	return someArray;
} 

/* Free a 2D array of ints of length rows */
void freeIntArray(int **someArray, const uint32_t rows) {
	uint32_t i;
	for (i=0; i<rows; i++) {free(someArray[i]);}
	free(someArray);
}




/* Parse a CSV file, return pointer to a 2d array of doubles */
double **csvparse(const char filePath[], const char delim, uint32_t * const restrict rows, uint32_t * const restrict columns){

	// File to open
	FILE *fp;
	fp=fopen(filePath,"r");

	if (fp==NULL){
		fprintf(stderr,"WARNING: file does not exist!\n");
		return NULL;

	} else {

		char c;
		uint32_t numColumns=0, maxColumns=0, numChars=0, maxChars=0, numRows=0;

		// Determine maximum number of characters per row, delimiters per row, and rows
		for(c=getc(fp); c != EOF; c=getc(fp)){
			numChars++;
			if (c==delim){numColumns++;}
			if (c=='\n'){
				numRows++;
				//if there is a trailing delimiter, don't add an extra column for it
				fseek(fp,-2,SEEK_CUR);
				if(getc(fp)!=delim){numColumns++;}
				fseek(fp,+1,SEEK_CUR);	
				// See if we have a new maximum, and reset the counters
				if (numChars>maxChars){maxChars=numChars;}
				if (numColumns>maxColumns){maxColumns=numColumns;}
				numChars=0;
				numColumns=0;
			}
		} 
		// If the last line isn't blank, add one more to the row counter
		fseek(fp,-1,SEEK_CUR);
		if(getc(fp)!='\n'){numRows++;}


		// Malloc space for the imported array
		double ** const restrict importedMatrix=mallocDoubleArray(numRows,maxColumns);
		*rows=numRows;
		*columns=maxColumns;				

		// For each line,
		uint32_t i=0, j=0, k, field[maxColumns+2];
		char str[maxChars+2];
		char *endp;
		rewind(fp);
		while(fgets(str,maxChars+2,fp)!=NULL){

			// identify the delimited fields,
			field[0]=0;
			for(k=1;k<maxColumns+2;k++){field[k]=maxChars+1;}
			for(k=0, numColumns=0; str[k]!='\0'; k++){
				if (str[k]==delim){
					str[k]='\0';
					field[numColumns+1]=k+1;
					numColumns++;
				} else if(str[k]=='\n'){str[k]='\0';}
			}

			// and perform operations on each field
			for (j=0;j<maxColumns;j++){
				importedMatrix[i][j]=strtod(&str[field[j]],&endp);
				if(endp==&str[field[j]]){importedMatrix[i][j]=NAN;}
			}
			i++;
		}
		fclose(fp);
		#ifdef DEBUG
			fprintf(stderr,"Maximum number of characters: %d\n", maxChars);
			fprintf(stderr,"Maximum number of delimiters: %d\n", maxColumns);
			fprintf(stderr,"Number of rows: %d\n", numRows);
		#endif 

		return importedMatrix;
	}

}


/* Parse a CSV file, return pointer to a 2d array of doubles */
double *csvparseflat(const char filePath[], const char delim, uint32_t * const restrict rows, uint32_t * const restrict columns){

	// File to open
	FILE *fp;
	fp=fopen(filePath,"r");

	if (fp==NULL){
		fprintf(stderr,"WARNING: file does not exist!\n");
		return NULL;

	} else {

		char c, cl;
		uint32_t numColumns=0, maxColumns=0, numChars=0, maxChars=0, numRows=0;

		// Determine maximum number of characters per row, delimiters per row, and rows
		for(c=getc(fp); c != EOF; c=getc(fp)){
			numChars++;
			if (c==delim){numColumns++;}
			if (c=='\n'){
				numRows++;
				numColumns++;
				fseek(fp,-2,SEEK_CUR);
				cl=getc(fp);
				//if there is a trailing delimiter, don't add an extra column for it
				if (cl==delim){
					numColumns--;
				//if there is an empty line, don't add an extra row for it
				} else if (cl=='\n'){
					numRows--;
				}
				fseek(fp,+1,SEEK_CUR);	
				// See if we have a new maximum, and reset the counters
				if (numChars>maxChars){maxChars=numChars;}
				if (numColumns>maxColumns){maxColumns=numColumns;}
				numChars=0;
				numColumns=0;
			}
		} 
		// If the last line isn't blank, add one more to the row counter
		fseek(fp,-1,SEEK_CUR);
		if(getc(fp)!='\n'){numRows++;}


		// Malloc space for the imported array
		double * const restrict importedMatrix=malloc(sizeof(double)*maxColumns*numRows);
		*rows=numRows;
		*columns=maxColumns;				

		// Read the array
		uint32_t i=0, j=0, k, field[maxColumns+2];
		char str[maxChars+2];
		char *endp;
		rewind(fp);
		// For each line,
		while(fgets(str,maxChars+2,fp)!=NULL){
			// if line is not empty,
			if (str[0]!='\n'){
				// identify the delimited fields,
				field[0]=0;
				for(k=1;k<maxColumns+2;k++){field[k]=maxChars+1;}
				for(k=0, numColumns=0; str[k]!='\0'; k++){
					if (str[k]==delim){
						str[k]='\0';
						field[numColumns+1]=k+1;
						numColumns++;
					} else if(str[k]=='\n'){str[k]='\0';}
				}

				// perform operations on each field
				for (j=0;j<maxColumns;j++){
					importedMatrix[j*numRows + i]=strtod(&str[field[j]],&endp);
					if(endp==&str[field[j]]){importedMatrix[j*numRows + i]=NAN;}
				}
				// and increment the row counter
				i++;
			}
		}
		fclose(fp);

		#ifdef DEBUG
			fprintf(stderr,"Maximum number of characters: %d\n", maxChars);
			fprintf(stderr,"Maximum number of delimiters: %d\n", maxColumns);
			fprintf(stderr,"Number of rows: %d\n", numRows);
		#endif
		return importedMatrix;
	}

}

uint32_t fprintfflat(FILE* fp, const double* data, const char delim, uint32_t rows, uint32_t  columns){
	//Print the array
	for (uint32_t i=0;i<rows;i++){
		for (uint32_t j=0;j<columns;j++){
			fprintf(fp, "%g%c", data[j*rows + i], delim);
		}
		// Delete the trailing delimiter and print a newline at the end of the row
		fseek(fp,-1,SEEK_CUR);
		fprintf(fp,"\n");
	}
	return 0;
}

uint32_t fprintfflatindex(FILE* fp, const double* data, const int index, const char delim, uint32_t rows, uint32_t  columns){
	//Print the array
	for (uint32_t i=0;i<rows;i++){
		fprintf(fp, "%i%c", index, delim);
		for (uint32_t j=0;j<columns;j++){
			fprintf(fp, "%g%c", data[j*rows + i], delim);
		}
		// Delete the trailing delimiter and print a newline at the end of the row
		fseek(fp,-1,SEEK_CUR);
		fprintf(fp,"\n");
	}
	return 0;
}

uint32_t csvwriteflat(const double* data, const char filePath[], const char* mode,  const char delim, uint32_t rows, uint32_t  columns){
	// File to open
	FILE *fp;
	fp=fopen(filePath, mode);

	//Print the array
	for (uint32_t i=0;i<rows;i++){
		for (uint32_t j=0;j<columns;j++){
			fprintf(fp, "%g%c", data[j*rows + i], delim);
		}
		// Delete the trailing delimiter and print a newline at the end of the row
		fseek(fp,-1,SEEK_CUR);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 0;
}


// Compute mean and variance with Knuth's algorithm
void Knuth_var(const double* const x, const uint32_t n, double* const restrict mean, double* const restrict var){
	double delta, mu=0, m2=0;

	for (uint32_t i=0; i<n; i++){
		delta = x[i] - mu;
		mu += delta / (i+1);
		m2 += delta*(x[i] - mu);
	}

	if (n>1){
		*mean = mu;
		*var = m2/(n-1);
	} else {
		*mean = x[0];
		*var = 0;
	}
}

// Compute mean and variance with Knuth's algorithm, ignoring NaN data
void Knuth_nanvar(const double* const x, const uint32_t n, double* const mean, double* var){
	double delta, mu=0, m2=0;
	uint32_t exists=0;

	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			delta = x[i] - mu;
			mu += delta / exists;
			m2 += delta*(x[i] - mu);
		}
	}

	if (exists>1){
		*mean = mu;
		*var = m2/(exists-1);
	} else if (exists==1 && n>1){
		*mean = mu;
		*var = 0;
	} else {
		*mean = NAN;
		*var = NAN;
	}
}

// Compute mean and standard deviation with Knuth's algorithm, ignoring NaN data
void Knuth_nanstd(const double* const x, const uint32_t n, double* const restrict mean, double* const restrict std){
	double delta, mu=0, m2=0;
	uint32_t exists=0;

	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			delta = x[i] - mu;
			mu += delta / exists;
			m2 += delta*(x[i] - mu);
		}
	}

	if (exists>1){
		*mean = mu;
		*std = sqrt(m2/(exists-1));
	} else if (exists==1 && n>1){
		*mean = mu;
		*std = 0;
	} else {
		*mean = NAN;
		*std = NAN;
	}
}

// Compute mean and standard deviation with Knuth's algorithm, ignoring NaN data
void Knuth_nanstderr(const double* const x, const uint32_t n, double* const restrict mean, double* const restrict StdErr){
	double delta, mu=0, m2=0;
	uint32_t exists=0;

	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			delta = x[i] - mu;
			mu += delta / exists;
			m2 += delta*(x[i] - mu);
		}
	}

	if (exists>1){
		*mean = mu;
		*StdErr = sqrt(m2/(exists-1)/exists);
	} else if (exists==1 && n>1){
		*mean = mu;
		*StdErr = 0;
	} else {
		*mean = NAN;
		*StdErr = NAN;
	}
}

// Compute sum of a double array, excluding NaNs
double nansum(const double* const x, const uint32_t n){
	double S=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			S+=x[i];
			exists=1;
		}
	}
	if (exists){
		return S;
	} else {
		return NAN;
	}
}

// Compute sum of the squares of a double array, excluding NaNs
double nansumSq(const double* const x, const uint32_t n){
	double S=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			S+=x[i]*x[i];
			exists=1;
		}
	}
	if (exists){
		return S;
	} else {
		return NAN;
	}
}


// Compute mean of a double array, excluding NaNs
double nanmean(const double* const x, const uint32_t n){
	double S=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S+=x[i];
		}
	}
	return S/exists;
}

// Compute variance of a double array
double nanvar(const double* const x, const uint32_t n){
	double c;
	if (n<1){
		return NAN;
	} else if (n<11) {;
		c=x[0];
	} else {
		c=nanmean(x,10);
	}
	if (isnan(c)){
	       c=0.0;
	}

	double S=0.0, S2=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S += x[i]-c;
			S2 += (x[i]-c) * (x[i]-c);
		}
	}
	return (S2 - (S*S)/exists)/(exists-1);
}

// Compute standard deviation of an array, excluding NaNs
double nanstd(const double* const x, const uint32_t n){
	double c;
	if (n<1){
		return NAN;
	} else if (n<11) {;
		c=x[0];
	} else {
		c=nanmean(x,10);
	}
	if (isnan(c)){
	       c=0.0;
	}

	double S=0.0, S2=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S += x[i]-c;
			S2 += (x[i]-c) * (x[i]-c);
		}
	}
	return sqrt((S2 - (S*S)/exists)/(exists-1));
}

// Compute mean and standard deviation of an array, excluding NaNs
int Offset_nanvar(const double* const x, const uint32_t n, double* const restrict mu, double* const restrict var){
	*mu=nanmean(x,n);
	double S=0.0, S2=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S += x[i]-*mu;
			S2 += (x[i]-*mu) * (x[i]-*mu);
		}
	}
	*var=(S2 - (S*S)/exists)/(exists-1);
	return 0;
}

// Compute mean and standard deviation of an array, excluding NaNs
int Offset_nanstd(const double* const x, const uint32_t n, double* const restrict mu, double* const restrict sigma){
	*mu=nanmean(x,n);
	double S=0.0, S2=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S += x[i]-*mu;
			S2 += (x[i]-*mu) * (x[i]-*mu);
		}
	}
	*sigma=sqrt((S2 - (S*S)/exists)/(exists-1));
	return 0;
}

// Compute mean and standard deviation of an array, excluding NaNs
int Offset_nanstderr(const double* const x, const uint32_t n, double* const restrict mu, double* const restrict StdErr){
	*mu=nanmean(x,n);
	double S=0.0, S2=0.0;
	uint32_t exists=0;
	for (uint32_t i=0; i<n; i++){
		if (!isnan(x[i])){
			exists++;
			S += x[i]-*mu;
			S2 += (x[i]-*mu) * (x[i]-*mu);
		}
	}
	*StdErr=sqrt((S2 - (S*S)/exists)/(exists-1)/exists);
	return 0;
}

// De-mean a double array, excluding NaNs
int normalize(double* restrict x, const uint32_t n){
	double mu=nanmean(x,n);
	for (uint32_t i=0; i<n; i++){
		x[i]-=mu;
	}
	return 0;
}

// Standardize a double array to zero mean and unit variance, excluding NaNs
int standardize(double* restrict x, const uint32_t n){
	double mu, sigma;
	Offset_nanstd(x,n,&mu,&sigma);
	if (sigma>0){
		for (uint32_t i=0; i<n; i++){
			x[i] = (x[i]-mu)/sigma;
		}
	} else {
		for (int i=0; i<n; i++){
			x[i] = (x[i]-mu);
		}
	}
	return 0;
}


//// Calculate an MSWD
//double mswd(const double* x, const double* sigma, const uint32_t n){
//	uint32_t i;
//	double s1 = 0,  s2 = 0, s3 = 0;
//	for(i=0; i<n; i++){
//		s1 += x[i] / (sigma[i]*sigma[i]);
//		s2 += 1 / (sigma[i]*sigma[i]);
//	}
//	double wx = s1/s2;
//
//	for(i=0; i<n; i++){
//		s3 += (x[i]-wx)*(x[i]-wx) / (sigma[i]*sigma[i]);
//	}
//
//	return s3 / (n-1);
//}

// Calculate a weigted mean, including MSWD.
int wmean(const double* x, const double* sigma, const uint32_t n, double* wx, double* wsigma, double* mswd){
	uint32_t i;
	double s0 = 0, s1 = 0,  s2 = 0, s3 = 0;
	for(i=0; i<n; i++){
		s0 += 1 / sigma[i];
		s1 += x[i] / (sigma[i]*sigma[i]);
		s2 += 1 / (sigma[i]*sigma[i]);
	}
	*wx = s1/s2;

	for(i=0; i<n; i++){
		s3 += (x[i] - *wx)*(x[i] - *wx) / (sigma[i]*sigma[i]);
	}

	*mswd = s3 / (n-1);
	*wsigma = sqrt(*mswd/s0);

	return 0;
}



// From wikipedia (public domain)
/* Comparison function. Receives two generic (void) pointers. */
int compare_ints(const void *p, const void *q){
	int x = *(const int *)p;
	int y = *(const int *)q;
	/* to avoid undefined behaviour through signed integer overflow,
	 *         avoid: return x - y; */
	int ret;
	if (x == y){
		ret = 0;
	} else {
		ret = (x < y) ? -1 : 1;
	}
	return ret;
}
 
/* Sort an array of n integers, pointed to by a. */
int sort_ints(int *a, size_t n){
	qsort(a, n, sizeof(int), compare_ints);
	return 0;
}

// Sort an array and delete nonunique elements
uint32_t unique_ints(int *a, const uint32_t np){
	size_t n=np;
	qsort(a, n, sizeof(int), compare_ints);
	uint32_t i, k=0;
	int last = a[0];
	for (i=1; i<n; i++){
		if (a[i]!=last){
			last=a[i];
			k++;
			a[k]=a[i];
		}
	}
	return k+1;
}


// From wikipedia (public domain)
/* Comparison function. Receives two generic (void) pointers. */
int compare_uints(const void *p, const void *q){
	uint32_t x = *(const uint32_t *)p;
	uint32_t y = *(const uint32_t *)q;
	/* to avoid undefined behaviour through signed integer overflow,
	 *         avoid: return x - y; */
	int ret;
	if (x == y){
		ret = 0;
	} else {
		ret = (x < y) ? -1 : 1;
	}
	return ret;
}
 
/* Sort an array of n integers, pointed to by a. */
int sort_uints(int *a, size_t n){
	qsort(a, n, sizeof(uint32_t), compare_uints);
	return 0;
}

// Sort an array and delete nonunique elements
uint32_t unique_uints(uint32_t *a, const uint32_t np){
	size_t n=np;
	qsort(a, n, sizeof(uint32_t), compare_uints);
	uint32_t i, k=0, last = a[0];
	for (i=1; i<n; i++){
		if (a[i]!=last){
			last=a[i];
			k++;
			a[k]=a[i];
		}
	}
	return k+1;
}

// Comparison function for doubles
int compare_doubles(const void *p, const void *q){
	double x = *(const double *)p;
	double y = *(const double *)q;

	int ret;
	if (x == y){
		ret = 0;
	} else {
		ret = (x < y) ? -1 : 1;
	}
	return ret;
}

// Comparison function for doubles
int compare_doubles_descending(const void *p, const void *q){
	double x = *(const double *)p;
	double y = *(const double *)q;

	int ret;
	if (x == y){
		ret = 0;
	} else {
		ret = (x < y) ? 1 : -1;
	}
	return ret;
}


// Sort an array of n doubles, pointed to by a
int sort_doubles(double *a, size_t n){
	qsort(a, n, sizeof(a[0]), compare_doubles);
	return 0;
}

#define sort_doubles_ascending	sort_doubles

// Sort an array of n doubles, pointed to by a
int sort_doubles_descending(double *a, size_t n){
	qsort(a, n, sizeof(a[0]), compare_doubles_descending);
	return 0;
}


