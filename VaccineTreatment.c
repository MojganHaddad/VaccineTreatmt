/***********************************************
File: VaccineTreatment.c
Date: 9/13/99
Coded by: Mojgan Haddad
************************************************/
#include <stdio.h>

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>

// -1 used to initialize variables
#define Init			-1

typedef char *string;

/* enumerated type corresponding to either uniform (U) or 
   triangular (T) pdf */
typedef enum{
	U, T
} distnT;

/* paramT is the type of array that will hold all the info from the 
parameter file that is read in. Each element in the array has a name 
(ex. mu), a min value, a peak value (irrelevant in a uniform dist'n, 
but use 0 as a placeholder in the input file), a max value, and 
a type (either U or T for uniform or triangular) */
   
typedef struct{
	string name;
	double min;
	double peak;
	double max;
	distnT dType;
} paramT;

/* Since we need the parameter values and the outcome variable values 
ranked for the sensitivity analysis, I made a small structure, such 
that each value can easily have a corresponding rank. */
   
typedef struct{
	double value;
	int rank;
} matrixDB;

// All function heads below */
static void Error(string msg);
static int RandomInteger(int low, int high);
static void Randomize(void);
static string CopyString(string s);
static void *GetBlock(size_t nbytes);
static distnT ConvertDistnToType(char dType);
static void FillHypercubeMatrix(matrixDB **lhcMatrix, paramT *pArray, int N, int K);
static void FillHypercube(int N, int element, matrixDB **lhcMatrix, paramT parameter);
static void FillUniformHypercube(int N, matrixDB *putArray, double min, double max);
static void FillTriangularHypercube(int N, matrixDB *putArray, double min, double peak, double max);
static void JumbleMatrixRow(matrixDB *row, int N);
static void GiveRankings(matrixDB* array, int N);
static int FindLowest(matrixDB* array, double least, int N);
static void GiveSameRanking(matrixDB* array, int N);
static void PrintHypercubeSamples(matrixDB **lhcMatrix, paramT *pArray, int N, int K);
static void CleanArray(double **timeMatrix, int nTimes, int boxes);
static void PrintOutTimeArray(double **timeMatrix, int nTimes, matrixDB **lhcMatrix, int sample, double **outcomeMatrix, double **outcomeMatrix1);//, double **outcomeMatrix2);
static void FillTimeArray(int sample, matrixDB **lhcMatrix, double **timeMatrix, int nTimes, int boxes);
static void FillInitialConditions(int sample, matrixDB **lhcMatrix, double *firstrow);
static void FillRestThroughRungeKutta(int sample, matrixDB **lhcMatrix, double **array, int nTimes, int boxes);
static void Runge_Kutta(double *inArray, double *derivArray, int n, double x, double h, double *outArray, int sample, matrixDB **lhcMatrix);
static void GetDerivatives(double *array, double *derivArray, int sample, matrixDB **lhcMatrix);
static void UncertaintyAnalysis(matrixDB **lhcMatrix, int N, int K);
static double GetVariance(matrixDB *col, int N, double mean);
static double GetMedian(matrixDB *col, int N);
static double GetTotal(matrixDB *col, int N);
static void EvaluatePRCC(matrixDB **lhcMatrix, paramT *pArray, int K, int N, int timePt, double **outcomeMatrix, string filename);
static double GetSignif(double prcc, int N);
static double GetPRCC(double **invMatrix, int i, int kPlusOne);
static void InvertMatrix(double **symMatrix, double **invMatrix, int size);
static double FillCijElement(matrixDB **lhcMatrix, int N, int i, int j, double mu);
static double GetCu(matrixDB **lhcMatrix, int timeBox, int sample);
static int CopyPertinentParameters(matrixDB **lhcMatrix, matrixDB **copyMatrix, paramT *pArray, int K, int N, int timePt, double **outcomeMatrix);

/* end of function prototypes */


/* global variables- c: # of sex partners/yr
					 pi: influx of sexually active people into the 
					     population
   Both c and pi are derived for each of the N hypercube samples- 
   they come from the prevalence and popthous parameter values. */
   		 				 
double c;
double pi;


// main program
main()
{
	FILE *infile;
	int i, j;
	char name[20];
	char dType;
	paramT *pArray;
	matrixDB **lhcMatrix;
	double **timeMatrix;
	double **outcomeMatrix;
	double **outcomeMatrix1;
	// double **outcomeMatrix2;
	int K;
	int N; 
	float min, peak, max;
	int nTimes, boxes;
	
	/* number of times Runge-Kutta is called- 
	   since the time increment is 0.0005 (see later in the code),
	   nTimes corresponds to 80000*0.0005 yrs -> 40 yrs */
	nTimes = 80000;
	
	/* the boxes correspond to:
		box0 = X
		box1 = Yv1
		box2 = Yv2
		box3 = Yw
		box4 = Yvw
		box5 = YvT
		box6 = YwT
		box7 = YvwT
		box8 = A
		box9 = ADd  : # AIDS Death
		box10 = Treated : # Treated cases
	*/
	boxes = 11;
	
	// Reading in of the input file with the parameter ranges
	infile = fopen("HIVparameters.txt", "r");
		
	if(infile == NULL) Error("Cannot find file named: HIVparameters");
	
	fscanf(infile, "LHC Samples:%d\n", &N);
	
	fscanf(infile, "Number of parameters:%d\n", &K);
	
	pArray = (paramT *) GetBlock(K*sizeof(paramT));
	
	for(i=0; i< K; i++){
		fscanf(infile, "%s %f %f %f %c\n", name, &min, &peak, &max, &dType);
		pArray[i].name = CopyString(name);
		pArray[i].min = (double) min;
		pArray[i].peak = (double) peak;
		pArray[i].max = (double) max;
		pArray[i].dType = ConvertDistnToType(dType);
	}
	
	fclose(infile);
	/* END of reading in input file */
	
	/* N = number of hypercube samples
	   K = number of parameters to be hypercubed */
	/* The way I stored all the hypercubed samples is via a two-dimensional array:
		(memory allocated below):
		lhcMatrix has K+1 columns (each column holds a parameter), and N rows (each row holds a lhc-ed sample)
		For example lhcMatrix[2][3] holds the third hypercubed sample of the parameter muA.
		Actually it holds the value (lhcMatrix[2][3].value) as well as its rank (lhcMatrix[2][3].rank).
		The nice thing about this setup then is that if you look at any row (across the K columns) you get the point
		in the hypercube (ie. all your parameter values for that run).
	*/ 
	
	lhcMatrix = (matrixDB **) GetBlock((K+1)*sizeof(matrixDB *));
	/* extra space in matrix for outcome variable... so after I choose an outcome variable,
	   I can place it alongside the parameter value and rank it as well */
	
	for(i = 0; i< K+1; i++){
		lhcMatrix[i] = (matrixDB *) GetBlock(N*sizeof(matrixDB));
		for(j=0; j<N; j++){
			lhcMatrix[i][j].value = Init;
			lhcMatrix[i][j].rank = Init;
		}	
	}
	
	// Fill the 2-D matrix with hypercubed values
	FillHypercubeMatrix(lhcMatrix, pArray, N, K);
	
	/* I create a 2-D time matrix- nTimes x boxes,
	   such that for each time, I keep the value of 
	   X, Yv1, Yv2, YvT, Yw1, Yw2, YwT, YvwT, and A */
	   
	timeMatrix = (double **) GetBlock(nTimes*sizeof(double *));

	for(i = 0; i< nTimes; i++){
		timeMatrix[i] = (double *) GetBlock(boxes*sizeof(double));	
	}
	
	/* I use the same timeMatrix for each of the N simulations 
	  (too much memory otherwise), so I like to clean the 
	  timeMatrix before reusing it*/
	
	/* OUTCOME VARIABLE matrix */
	outcomeMatrix = (double **) GetBlock(N*sizeof(double *));

	for(i = 0; i< N; i++){
		outcomeMatrix[i] = (double *) GetBlock(nTimes*sizeof(double));	
	}
	
	outcomeMatrix1 = (double **) GetBlock(N*sizeof(double *));

	for(i = 0; i< N; i++){
		outcomeMatrix1[i] = (double *) GetBlock(nTimes*sizeof(double));	
	}   
	/*
	outcomeMatrix2 = (double **) GetBlock(N*sizeof(double *));

	for(i = 0; i< N; i++){
		outcomeMatrix2[i] = (double *) GetBlock(nTimes*sizeof(double));	
	}  */
	
	for(i=0; i < N; i++){
		CleanArray(timeMatrix, nTimes, boxes);
		
		FillTimeArray(i, lhcMatrix, timeMatrix, nTimes, boxes);
		
		PrintOutTimeArray(timeMatrix, nTimes, lhcMatrix, i, outcomeMatrix, outcomeMatrix1);// outcomeMatrix2);*/
	}
	
	// print out the hypercube samples
	PrintHypercubeSamples(lhcMatrix, pArray, N, K);
	
	for(j=0; j< nTimes; j+=1000){
		EvaluatePRCC(lhcMatrix, pArray, K, N, j, outcomeMatrix, "AnnualPRCC");
	    EvaluatePRCC(lhcMatrix, pArray, K, N, j, outcomeMatrix1, "TreatedPRCC");
	//EvaluatePRCC(lhcMatrix, pArray, K, N, outcomeMatrix2, "EffectivenessPRCC");  */
	}
} // end of main program

static void PrintOutTimeArray(double **timeMatrix, int nTimes, matrixDB **lhcMatrix, int sample, double **outcomeMatrix, double **outcomeMatrix1)//, double **outcomeMatrix2)
{
	FILE* outfile;
	FILE* outfile2;
	FILE* outfile3;
	int j, box;
	double deaths, deathperHunK, tot, NrVaccinated, NrTreated;
	double deathsW, deathsN, CuTrtd, Avrtd, hivp, prevalence;
	
	outfile = fopen("AidsDeathsPer100000", "a");
	if(outfile == NULL) Error("Cannot open 'c' for writing");
	
	fprintf(outfile, "AnnDeathsPer100000: ");
	for(j=0; j< nTimes; j+=1000){
		deaths = timeMatrix[j][8]*lhcMatrix[2][sample].value; // A*muA
		if (deaths < 0)  deaths = 0;  // avoid negative death rate
		// total population
		tot = 0;
		for (box=0; box<9; box++) tot += timeMatrix[j][box];
		deathperHunK = deaths*100000/tot; 
		fprintf(outfile, "%g ", deathperHunK);
		outcomeMatrix[sample][j] = deathperHunK;
	};
	fprintf(outfile, "\n");
	fclose(outfile);
	
	outfile = fopen("TotalNrVaccinated", "a");
	if(outfile == NULL) Error("Cannot open 'TotalNrVaccinated' for writing");
	outfile2 = fopen("TotalNrTreated", "a");
	if(outfile2 == NULL) Error("Cannot open 'TotalNrTreated' for writing");
	outfile3 = fopen("Prevalence", "a");
	if(outfile3 == NULL) Error("Cannot open 'Prevalence' for writing");
	
	
	fprintf(outfile, "TotalNrVaccinated: ");
	fprintf(outfile2, "TotalNrTreated: ");
	fprintf(outfile3, "Prevalence: ");
	
	for(j=0; j< nTimes; j+=1000){
			NrVaccinated = timeMatrix[j][1];
			NrTreated = timeMatrix[j][10];
			fprintf(outfile, "%g ", NrVaccinated);
			if (NrTreated < 0)  NrTreated = 0;
			fprintf(outfile2, "%g ", NrTreated*1000);
			outcomeMatrix1[sample][j] = NrTreated*1000;
			
			if (lhcMatrix[0][sample].value == 0) {
			   hivp = timeMatrix[j][3] + timeMatrix[j][6] + timeMatrix[j][8];
			   tot = timeMatrix[j][0] + timeMatrix[j][8] + timeMatrix[j][3] + timeMatrix[j][6];
			   prevalence = hivp / tot;
			   fprintf(outfile3, "%g ", prevalence);
			   }
	}
	fprintf(outfile, "\n");
	fclose(outfile);
	fprintf(outfile2, "\n");
	fclose(outfile2);
	fprintf(outfile3, "\n");
	fclose(outfile3);
	
	outfile = fopen("CumDeathsPrevented", "a");
	if(outfile == NULL) Error("Cannot open 'CumDeathsPrevented' for writing");
	outfile2 = fopen("CumAIDSdeathsNVaccine", "a");
	if(outfile2 == NULL) Error("Cannot open 'CumAIDSdeathsNTrtmnt' for writing");
	outfile3 = fopen("CumDthPrvPTrd", "a");
	if(outfile3 == NULL) Error("Cannot open 'CumDthPrvPTrd' for writing");
	
	fprintf(outfile, "CumDeathsPrevented: ");
	fprintf(outfile2, "CumDeathsNTrtmnt: ");
	fprintf(outfile3, "CumDthPrvPerTrdCase: ");
	
	for(j=0; j< nTimes; j+=1000){
			deathsW = timeMatrix[j][9];
			deathsN = GetCu(lhcMatrix, j, sample);
			CuTrtd = timeMatrix[j][5] + timeMatrix[j][6] + timeMatrix[j][7];
			Avrtd = deathsW-deathsN;
			fprintf(outfile, "%g ", Avrtd);
			fprintf(outfile2, "%g ", deathsN);
			if (CuTrtd != 0) fprintf(outfile3, "%g ", Avrtd/CuTrtd);
	}
	fprintf(outfile, "\n");
	fclose(outfile);
	fprintf(outfile2, "\n");
	fclose(outfile2);
	fprintf(outfile3, "\n");
	fclose(outfile3);
		
	/*	
	outfile = fopen("Et", "a");
	if(outfile == NULL) Error("Cannot open 'Et' for writing");
	
	fprintf(outfile, "Et: ");
	for(j=0; j< nTimes; j+=1000){
		ADd = timeMatrix[j][9];
		Cu = GetCu(lhcMatrix, j, sample);
		Et = 1 - ADd/Cu;
		fprintf(outfile, "%g ", Et);
		// outcomeMatrix2[sample][j] = Et;
	};
	fprintf(outfile, "\n");
	fclose(outfile);
	*/
} // end of PrintOutTimeArray


static void FillTimeArray(int sample, matrixDB **lhcMatrix, double **timeMatrix, int nTimes, int boxes)
{
	/* To fill in the time array for a specific run, get the 
	starting equilib values and then run the scenario through time */
	
	FillInitialConditions(sample, lhcMatrix, timeMatrix[0]);
	
	FillRestThroughRungeKutta(sample, lhcMatrix, timeMatrix, nTimes, boxes);
}

static void FillInitialConditions(int sample, matrixDB **lhcMatrix, double *firstrow)
{
	double p, betaW, mu, VwU, Ro, muA, prev, pop;

	// initial conditions are at equilib without vaccination and treatment	
	p = lhcMatrix[0][sample].value;
	mu = lhcMatrix[1][sample].value;
	muA = lhcMatrix[2][sample].value;
	
	betaW = lhcMatrix[4][sample].value;
	VwU = lhcMatrix[15][sample].value;
	
	prev = lhcMatrix[22][sample].value;
	pop = lhcMatrix[23][sample].value;
	
	Ro = 1/(1-prev);
	c = Ro*(mu + VwU)/betaW;
	pi = pop * (mu + VwU*prev);

	printf("\nRunning Simulation # %d....\n", sample + 1);
	
	// Initial condition: mass vaccination, no treatment
	firstrow[0] = (1-p)*(1-prev)*pop;	// X
	firstrow[1] = p*(1-prev)*pop;	// Yv1
	firstrow[2] = 0; 		// Yv2
	firstrow[3] = pop*prev; // Yw
	firstrow[4] = 0; 		// Yvw
	firstrow[5] = 0; 		// YvT
	firstrow[6] = 0; 		// YwT
	firstrow[7] = 0; 		// YvwT
	firstrow[8] = firstrow[3]*VwU/(mu + muA); 	// A
	firstrow[9] = 0; // ADd
	firstrow[10] = 0; // Treated
} // end of FillInitialConditions


static void FillRestThroughRungeKutta(int sample, matrixDB **lhcMatrix, double **array, int nTimes, int boxes)
{
	double *derivArray;
	double x, timeIncrement;
	int i;
	
	// This function calls the Runge-Kutta function for each of the nTimes
	
	derivArray = (double *) GetBlock(boxes*sizeof(double));
	
	x = 0.0;
	
	timeIncrement = 0.0005;
	for(i=0; i<nTimes-1; i++){	
		GetDerivatives(array[i], derivArray, sample, lhcMatrix);
		Runge_Kutta(array[i], derivArray, 11, x, timeIncrement, array[i+1], sample, lhcMatrix);
		x += timeIncrement;
		if(i%1000 == 0) printf(".");
	}
	printf("\n");
	free(derivArray);
}

static void Runge_Kutta(double *inArray, double *derivArray, int n, double x, double h, double *outArray, int sample, matrixDB **lhcMatrix)
{
	int i;
	double xh, hh, h6;
	double *dym;
	double *dyt;
	double *yt;
	
	// Doing the Runge-Kutta- copied from Numerical Recipes...
	
	dym = (double *) GetBlock((n*sizeof(double)));
	dyt = (double *) GetBlock((n*sizeof(double)));
	yt = (double *) GetBlock((n*sizeof(double)));

	hh = h*0.5;
	h6 = h/6.0;
	xh = x+hh;
	
	for(i = 0; i < n; i++){
		yt[i] = inArray[i]+hh*derivArray[i];
	}
	GetDerivatives(yt, dyt, sample, lhcMatrix);
	
	for(i=0; i<n; i++){
		yt[i] = inArray[i]+hh*dyt[i];
	}
	GetDerivatives(yt, dym, sample, lhcMatrix);
	
	for(i=0; i<n; i++){
		yt[i] = inArray[i]+h*dym[i];
		dym[i] = dyt[i]+dym[i];
	}
	GetDerivatives(yt, dyt, sample, lhcMatrix);
	
	for(i=0; i<n; i++){
		outArray[i] = inArray[i]+h6*(derivArray[i]+dyt[i]+2.0*dym[i]);
	}
	
	free(dym);
	free(dyt);
	free(yt);
}


static void GetDerivatives(double *array, double *derivArray, int sample, matrixDB **lhcMatrix)
{
	double X, Yv1, Yv2, Yw, Yvw, YvT, YwT, YvwT, A;
	double dxdt, dyv1dt, dyv2dt, dywdt, dyvwdt, dTrdt;
	double dywtdt, dyvtdt,dyvwtdt, dadt, SIGMAv, SIGMAw, SIGMAvw;
	double BwU, BwT, Bv1, Bv2, BvT, BvwU, BvwT, LaV, LaW;
	double VvAIDS, Vv, VvU, VvT, VwU, VwT, VvwU, VvwT, Fw, Fv, Fvw;
	double Pv, phi, Trd, ADd, dADddt, TOTAL, p, mu, muA;
	
	/* get access to the latin hypercubed samples stored in the 
	   2-d lhcMatrix for simualtion # 'sample' */
	
	p = lhcMatrix[0][sample].value;
	mu = lhcMatrix[1][sample].value;
	muA = lhcMatrix[2][sample].value;
	phi = lhcMatrix[3][sample].value;
	
	// Transmission probabilities : Bs
	BwU = lhcMatrix[4][sample].value;
	BwT = lhcMatrix[5][sample].value;
	Bv1 = lhcMatrix[6][sample].value;
	Bv2 = lhcMatrix[7][sample].value;
	BvT = lhcMatrix[8][sample].value;
	BvwU = lhcMatrix[9][sample].value;
	BvwT = lhcMatrix[10][sample].value;
	
	// Progression rates : Vs
	
	//Vv = lhcMatrix[12][sample].value;
	//VvU = Vv * VvAIDS / (Vv - VvAIDS);
	VvAIDS = lhcMatrix[11][sample].value;
	Pv = lhcMatrix[12][sample].value;
	Vv = - log (1 - Pv) / 25;
	VvU = - log (1 - (VvAIDS / Pv)) / 25;
	lhcMatrix[13][sample].value = VvU;
	VwU = lhcMatrix[15][sample].value;
	VwT = lhcMatrix[16][sample].value * VwU; // VwTAlpha
	VvT = lhcMatrix[14][sample].value * VwT; // VvTAlpha
	VvwU = lhcMatrix[17][sample].value * VwU; // VvwUAlpha
	VvwT = lhcMatrix[18][sample].value * VwT; // VvwTAlpha
	
	// Fraction of treatment : Fs
	SIGMAw = lhcMatrix[19][sample].value;
	Fw = SIGMAw * (mu + VwU) / (1 - SIGMAw);
	SIGMAv = lhcMatrix[20][sample].value;
	Fv = SIGMAv * (mu + VvU) / (1 - SIGMAv);
	SIGMAvw = lhcMatrix[21][sample].value;
	Fvw = SIGMAvw * (mu + VvwU) / (1 - SIGMAvw);
	
	/* get the values in each of the boxes */
	X = array[0];
	Yv1 = array[1];
	Yv2 = array[2];
	Yw = array[3];
	Yvw = array[4];
	YvT = array[5];
	YwT = array[6];
	YvwT = array[7];
	A = array[8];
	ADd = array[9];
	Trd = array[10];
	
	 // total sexually active population
	TOTAL = X + Yv1 + Yv2 + Yw + Yvw + YvT + YwT + YvwT;
	
	// Transmissibility-combined : lambda Vaccine & lambda Wild-type
	LaV = (Bv1*Yv1 + Bv2*Yv2 + BvT*YvT) / TOTAL;
	LaW = (BwU*Yw + BwT*YwT + BvwU*Yvw + BvwT*YvwT) / TOTAL;
	
	
	// the model's differential equations
	dxdt = (1-p)*pi - mu*X - c*LaV*X - c*LaW*X;
	dyv1dt = p*pi + c*LaV*X - mu*Yv1 - (1-phi)*c*LaW*Yv1 - Vv*Yv1;
	dyv2dt = Vv*Yv1 - mu*Yv2 - Fv*Yv2 - (1-phi)*c*LaW*Yv2 - VvU*Yv2;
	dywdt = c*LaW*X - mu*Yw - Fw*Yw - VwU*Yw;
	dyvwdt = (1-phi)*c*LaW*Yv1 + (1-phi)*c*LaW*Yv2 - mu*Yvw - Fvw*Yvw - VvwU*Yvw;
	dyvtdt = Fv*Yv2 - mu*YvT - (1-phi)*c*LaW*YvT - VvT*YvT;
	dywtdt = Fw*Yw - mu*YwT - VwT*YwT;
	dyvwtdt = Fvw*Yvw + (1-phi)*c*LaW*YvT - mu*YvwT - VvwT*YvwT;
	dadt = VwU*Yw + VvU*Yv2 + VvwU*Yvw + VvT*YvT + VwT*YwT + VvwT*YvwT - mu*A - muA*A;
	dADddt = muA*A;
	dTrdt = Fv*Yv2 + Fw*Yw + Fvw*Yvw - VvT*YvT - VwT*YwT - VvwT*YvwT - mu*Trd;

	derivArray[0] = dxdt;
	derivArray[1] = dyv1dt;
	derivArray[2] = dyv2dt;
	derivArray[3] = dywdt;
	derivArray[4] = dyvwdt;
	derivArray[5] = dyvtdt;
	derivArray[6] = dywtdt;
	derivArray[7] = dyvwtdt;
	derivArray[8] = dadt;
	derivArray[9] = dADddt;
	derivArray[10] = dTrdt;

} // end of GetDerivatives


static void CleanArray(double **timeMatrix, int nTimes, int boxes)
{
	int i, j;
	
	for(i = 0; i< nTimes; i++){
		for(j=0; j<boxes; j++){
			timeMatrix[i][j] = Init;
		}	
	}
}

static void *GetBlock(size_t nbytes)
{
	void *result;
	
	result = (void *)malloc(nbytes);
	if(result == NULL) Error("No memory available");
	return(result);
}

static void Error(string msg)
{
	va_list args;
	
	va_start(args, msg);
	fprintf(stderr, "Error: ");
	vfprintf(stderr, msg, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(1);
}

static string CopyString(string s)
{
	string newstr;
	
	if(s==NULL) Error("NULL string passed to CopyString");
	newstr = (string) GetBlock(strlen(s) + 1);
	strcpy(newstr, s);
	return(newstr);
}

static void Randomize(void)
{
	srand((int) time(NULL));
}

static int RandomInteger(int low, int high)
{
	int k;
	double d;
	
	d = (double) rand()/((double) RAND_MAX + 1);
	k = (int) (d*(high - low + 1));
	return(low + k);
}

static distnT ConvertDistnToType(char dType)
{
	switch(dType){
		case 'U': return U;
		case 'T': return T;
		default: Error("Check parameter file: distribution type must be 'T' or 'U'\n"); return(U);
	}
}

static void FillHypercubeMatrix(matrixDB **lhcMatrix, paramT *pArray, int N, int K)
{	
	int element;
	
	// cycle through the parameters and hypercube each one of them
	Randomize();
	for(element = 0; element< K; element++){
		printf("Reading in Parameter %d from input file...\n", element+1);
		FillHypercube(N, element, lhcMatrix, pArray[element]);
	}
}

static void FillHypercube(int N, int element, matrixDB **lhcMatrix, paramT parameter)
{
	/* hypercube by pdf, and then jumble up the samples, 
	   since the samples were initially ordered */
	switch(parameter.dType){
		case U: FillUniformHypercube(N, lhcMatrix[element], parameter.min, parameter.max); break;
		case T: FillTriangularHypercube(N, lhcMatrix[element], parameter.min, parameter.peak, parameter.max); break;
		default: Error("Parameter's distn type must be 'T' or 'U'."); break;
	}
	
	
	JumbleMatrixRow(lhcMatrix[element], N);
	// also give the samples ranks
	if(parameter.min != parameter.max){
		GiveRankings(lhcMatrix[element], N);
	} else {
		GiveSameRanking(lhcMatrix[element], N);
	}
}

static void GiveSameRanking(matrixDB* array, int N)
{
	int i;
	
	for(i=0; i< N; i++){
		array[i].rank = 1;
	}
}


static void GiveRankings(matrixDB* array, int N)
{
	double lastLowest;
	int i, index;
	
	lastLowest = Init;
	
	for(i = 0; i< N; i++){
		index = FindLowest(array, lastLowest, N);
		array[index].rank = i+1;
		lastLowest = array[index].value;
	}
}

static int FindLowest(matrixDB* array, double least, int N)
{
	double max;
	int i, index;
	
	max = 100000;
	for(i=0; i<N; i++){
		if(array[i].value < max && array[i].value > least){
			index = i;
			max = array[i].value;
		}
	}
	return(index);
}


static void FillTriangularHypercube(int N, matrixDB *putArray, double min, double peak, double max)
{
	double base, height, volumeSlice, ximin, ximax;
	double m1, m2, m, b1, b2, b, xintoInverse;
	int i;
	
	/* dices up a triangular distn, and goes from min to max, taking 
	   the midpoint from each of the N diced-up segments */
	base = max - min;
	height = 2.0/base;
	volumeSlice = 1.0/N;
	 
	ximin = min;
	m1 = height/(peak - min);
	m2 = -height/(max - peak);
	
	b1 = -height*min/(peak-min);
	b2 = height*max/(max-peak);
	
	m = m1;
	b = b1;
	
	for(i = 0; i< N; i++){
		xintoInverse = m*ximin*ximin/2 + b*ximin + volumeSlice;
		ximax = (-b + sqrt(b*b + 2*m*xintoInverse))/m;
		putArray[i].value = (ximin + ximax)/2;
		ximin = ximax;
		if(ximin > peak){
			ximin = peak;
			m = m2;
			b = b2;
			peak = 2*max;
		}
	}
}


static void FillUniformHypercube(int N, matrixDB *putArray, double min, double max)
{
	double incr;
	int i; 
	
	/* this one dices into N equivolume parts and picks the midpoints
	   from each segment */
	if(max == min){
		for(i=0; i< N; i++){
			putArray[i].value = min;
		}
	} else {
		incr = (max-min)/N;
		
		min += incr/2;
	
		for(i = 0; i< N; i++){
			putArray[i].value = min;
			min += incr;
		}
	}
}

static void JumbleMatrixRow(matrixDB *row, int N)
{
	int i, loc;
	double temp;
	
	// just jumbles up the ordered samples
	for(i=0; i<N; i++){
		loc = RandomInteger(i, N-1);
		temp = row[i].value;
		row[i].value = row[loc].value;
		row[loc].value = temp;
	}
}

static void PrintHypercubeSamples(matrixDB **lhcMatrix, paramT *pArray, int N, int K)
{
	int i, j, l;
	FILE *outfile;
	
	// just in case you want to see the samples in an output file
	outfile = fopen("LHCmatrixValues", "w");
	if(outfile == NULL) Error("Not enough memory available.");
	
	for(l=0; l<K; l++){
		fprintf(outfile, "%s ", pArray[l].name);
	}
	fprintf(outfile, "\n");
	
	for(i=0; i<N; i++){
		for(j=0; j<K; j++){
			fprintf(outfile, "%g ", lhcMatrix[j][i].value);
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
}


// Uncertainty Analysis
static void UncertaintyAnalysis(matrixDB **lhcMatrix, int N, int K)
{
	FILE* outfile;
	double totalOutcome, meanOutcome, medianOutcome, varOutcome;
	
	GiveRankings(lhcMatrix[K], N);

	totalOutcome = GetTotal(lhcMatrix[K], N);
	meanOutcome = totalOutcome/N;
	medianOutcome = GetMedian(lhcMatrix[K], N);
	varOutcome = GetVariance(lhcMatrix[K], N, meanOutcome);

	outfile = fopen("UncertaintyAnalysis", "w");
	fprintf(outfile, "Mean: %g Median: %g Variance: %g\n", meanOutcome, medianOutcome, varOutcome);

	fclose(outfile);

}

static double GetVariance(matrixDB *col, int N, double mean)
{
	double numer, diff, var;
	int i;
	
	numer = 0;
	for(i=0; i< N; i++){
		diff = col[i].value - mean;
		numer += diff*diff;
	}
	var = numer/N;
	return(var);
}

static double GetMedian(matrixDB *col, int N)
{
	int i;
	
	for(i=0; i< N; i++){
		if(col[i].rank == N/2) return(col[i].value);
	}
	return(Init);
}

static double GetTotal(matrixDB *col, int N)
{
	int i;
	double sum;
	
	sum = 0;
	
	for(i=0; i< N; i++){
		sum += col[i].value;
	}
	return(sum);
}

// The sensitivity analysis...see Sally's LHC paper for specifics as to Cij element etc
static void EvaluatePRCC(matrixDB **lhcMatrix, paramT *pArray, int K, int N, int timePt, double **outcomeMatrix, string filename)
{
	double **symMatrix;
	int n, p, i, j, newK;
	double mu, prcc;
	double **invMatrix;
	FILE* outfile;
	// FILE* outfile2;
	matrixDB **copyMatrix;
	
	
	/* copy the lhcMatrix into another 2-D matrix to get rid of 
	   the parameters that are constant (ex: prev) before doing 
	   the sensitivity analysis */
	
	copyMatrix = (matrixDB **) GetBlock((K+1)*sizeof(matrixDB *));
	
	for(i = 0; i< (K+1); i++){
		copyMatrix[i] = (matrixDB *) GetBlock(N*sizeof(matrixDB));
		for(j=0; j<N; j++){
			copyMatrix[i][j].value = Init;
			copyMatrix[i][j].rank = Init;
		}	
	}
	
	// Copy only the parameters that have a range (ie not constant)
	newK = CopyPertinentParameters(lhcMatrix, copyMatrix, pArray, K, N, timePt, outcomeMatrix);
		
	/*outfile2 = fopen("CopyCheck", "a");
	
	for(i = 0; i< N; i++){
		fprintf(outfile2, "%d ", i);
		for(j=0; j<newK+1; j++){
			fprintf(outfile2, "%g ", copyMatrix[j][i].value);
		}
		fprintf(outfile2, "\n");
	fclose(outfile2); */
	printf("New number of parameters for sensitivity analysis: %d\n", newK);
		
		
		symMatrix = (double**) GetBlock((newK+1)*sizeof(double *));
	
		for(n = 0; n< (newK + 1); n++){
			symMatrix[n] = (double *) GetBlock((newK+1)*sizeof(double));
			for(p=0; p<(newK+1); p++){
				symMatrix[n][p] = Init;
			}	
		}

		for(i = 0; i< newK+1; i++){
			GiveRankings(copyMatrix[i], N);
		}
	
		mu = (double)(1+N)/2; 
	
		for(i = 0; i< newK+1; i++){
			for(j=0; j< newK+1; j++){
				symMatrix[i][j] = FillCijElement(copyMatrix, N, i, j, mu);
			}
		}
	
		invMatrix = (double**) GetBlock((newK+1)*sizeof(double *));
	
		for(n = 0; n< (newK + 1); n++){
			invMatrix[n] = (double *) GetBlock((newK+1)*sizeof(double));
			for(p=0; p<(newK+1); p++){
				invMatrix[n][p] = Init;
			}	
		}
	
	/* see Sally's paper again- inverting the matrix */
		InvertMatrix(symMatrix, invMatrix, newK+1);
		
		outfile = fopen(filename, "a");
		if(outfile == NULL) Error("Can't open prcc file for writing.");
		
		// print a line of parameter names
	if (timePt == 0){
		for(i=0; i< K; i++){
			fprintf(outfile, "%s ", pArray[i].name);
			}
		fprintf(outfile, "\n");
		}

		for(i=0; i< newK; i++){
			prcc = GetPRCC(invMatrix, i, newK);
			fprintf(outfile, "%g ", prcc);
				
		}
		fprintf(outfile, "\n");
		fclose(outfile);
		
} // end of EvaluatePRCC

static int CopyPertinentParameters(matrixDB **lhcMatrix, matrixDB **copyMatrix, paramT *pArray, int K, int N, int timePt, double **outcomeMatrix)
{
	int newK, i, j, m;
	
	newK = 0;
	
	for(i=0; i<K; i++){
		if(pArray[i].min != pArray[i].max){
			for(j=0; j<N; j++){
				copyMatrix[newK][j].value = lhcMatrix[i][j].value;
				}
		newK++;
		}
		else  pArray[i].name = " ";
			
		} 

	
	/* determines which outcome variable */
	for(m=0; m<N; m++){
		copyMatrix[newK][m].value = outcomeMatrix[m][timePt];
	}

	return(newK);
} // end of CopyPertinentParameters

static double GetSignif(double prcc, int N)
{
	double signif;
	
	signif = prcc*(sqrt((N-2)/(1-prcc)));

	return(signif);
}

static double GetPRCC(double **invMatrix, int i, int kPlusOne)
{
	double prccNum, prccDenom, prcc;
	
	prccNum = -1*(invMatrix[i][kPlusOne]);

	prccDenom = sqrt(invMatrix[i][i]*invMatrix[kPlusOne][kPlusOne]);

	prcc = prccNum/prccDenom;
	
	return(prcc);

}

static void InvertMatrix(double **symMatrix, double **invMatrix, int size)
{
	/* hayley's code */
 	int r, c, k;
    double tem, tem2, tem3, temp[100][100];
    FILE *f2;

	for (r=0; r<size; r++){
   		for (c=0; c<size; c++){
      		temp[r][c] = symMatrix[r][c];
       		if (r==c){
       			invMatrix[r][c]=1.0;
       		} else {
       			invMatrix[r][c]=0.0;
       		}
   		 }
   	}
   	
 	for (k=0; k<size; k++){
	    tem=temp[k][k];
   		for (c=0; c<size; c++){
      		temp[k][c]/=tem;
      		invMatrix[k][c]/=tem;
   		}
    	if (k!= size-1){
      		for (r=k+1; r<size; r++){
	        	tem=temp[r][k];
   				for (c=0; c<size; c++){
            		tem2=temp[k][c];
               		tem3=invMatrix[k][c];
   					temp[r][c]-=tem2*tem;
         			invMatrix[r][c]-=tem3*tem;
      			}
        	}
     	}
  	}
  	
 	for (k=size-1; k>=1; k--){
  		for (r=k-1; r>=0; r--){
      		tem=temp[r][k];
      		for (c=0; c<size; c++){
          		tem2=temp[k][c];
           		tem3=invMatrix[k][c];
   				temp[r][c]-=tem*tem2;
      	   		invMatrix[r][c]-=tem*tem3;
        	}
     	}
  	}
	
   	f2=fopen("inv_mat.txt","w");
   	
  	for (r=0; r<size; r++){
   		for (c=0; c< size; c++){
      		fprintf(f2,"%4.2f ",symMatrix[r][c]);
         	if (c== size-1) fprintf(f2," | ");
     	}
     	for (c=0; c<size; c++){
      		fprintf(f2,"%4.2f ",invMatrix[r][c]);
        	if (c== size-1) fprintf(f2,"\n");
     	}
 	}
   	fclose(f2);

} // end of InvertMatrix


static double FillCijElement(matrixDB **lhcMatrix, int N, int i, int j, double mu)
{
	double sum, numerator, sum1, sum2, denominator, Cij;
	int t, s;
	
	/* numerator of Cij expression */
	sum = 0;
	for(t = 0; t< N; t++){
		sum += (lhcMatrix[i][t].rank - mu)*(lhcMatrix[j][t].rank - mu);
	}
	numerator = sum;
	
	/* denominator of Cij expression */
	sum1 = 0;
	for(t=0; t< N; t++){
		sum1 += (lhcMatrix[i][t].rank - mu)*(lhcMatrix[i][t].rank - mu);
	}
	sum2 = 0;
	for(s=0; s< N; s++){
		sum2 += (lhcMatrix[j][s].rank - mu)*(lhcMatrix[j][s].rank - mu);
	}

	denominator = sqrt(sum1*sum2);
	
	/* division */
	Cij = numerator/denominator;
	
	return(Cij);
	
}

static double GetCu(matrixDB **lhcMatrix, int timeBox, int sample)
{
	double VwA, mu, muA, yw, a;
	double dCudt, Cut, prev, pop;
	
	mu = lhcMatrix[1][sample].value;
	muA = lhcMatrix[2][sample].value;
	VwA = lhcMatrix[14][sample].value;
	prev = lhcMatrix[23][sample].value;
	pop = lhcMatrix[24][sample].value;
	
	yw = prev * pop;
	a = yw*VwA/(mu + muA);
	
	dCudt = muA*a;

	Cut = dCudt*timeBox*0.0005;
	return(Cut);

}