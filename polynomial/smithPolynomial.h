int M, N, rank, consistent;
int maxDegree;
int maxEntry;
float precision = 0.0001;
float precision2 = 0.1;

void leastEntryAlgo(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);

void initializeSize(int *, char (*));
void changeMaxDegree(int *, int *);


// void matxvecDiff(float (*)[N][maxDegree+1], float (*)[maxDegree+1], float (*)[N][maxDegree+1]);

void matNNxmatNN(float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1]);
void matMMxmatMM(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);

void matNNxmatNM(float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);
void matNMxmatMM(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);

/***********************************MATRIX_OPERATIONS_H***************************************/

void printAll(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], 
			  float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], 
			  float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], 
			  float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], 
			  float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], 
			  float (*)[maxDegree+1]);

void init(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], 
		  float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], 
		  float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], 
		  float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], 
		  float (*)[maxDegree+1]);

void initializeA(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);
void initializeP(float (*)[N][maxDegree+1]);
void initializeQ(float (*)[M][maxDegree+1]);
void initializeb(float (*));

void rowOperations1(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int, float);
void columnOperations2(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], int, int, float (*));
void rowOperations2(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int, int, float(*));
void rowOperations3(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int, int);
void columnOperations3(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], int, int);

void type1rowM(float (*)[M][maxDegree+1], int, float);
void type2rowM(float (*)[M][maxDegree+1], int, int, float (*), int);
void type3rowM(float (*)[M][maxDegree+1], int, int);

void type1rowN(float (*)[N][maxDegree+1], int, float);
void type2rowN(float (*)[N][maxDegree+1], int, int, float (*), int);
void type3rowN(float (*)[N][maxDegree+1], int, int);

void type1col(float (*)[M][maxDegree+1], int, int, float);
void type2col(float (*)[M][maxDegree+1], int, int, int, float(*), int);
void type3col(float (*)[M][maxDegree+1], int, int, int);

void transposeN(float (*)[N][maxDegree+1]);
void transposeM(float (*)[M][maxDegree+1]);

void print2ArrayM(float (*)[M][maxDegree+1], int);
void print2ArrayN(float (*)[N][maxDegree+1], int);
void printArray(float (*)[maxDegree+1], int);

void printDiff(float (*)[N][maxDegree+1], int);
void printDiffRow(float (*)[maxDegree+1]);

void makeAllMonic(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1]);


/***********************************POLYNOMIAL_OPERATIONS_H***************************************/

void clearZeroes(float (*));
void clearZeroes2(float (*));
void add(float (*), float (*));
void subtract(float (*), float (*));
void multiply(float (*), float (*), float (*));

void scale(float (*), float);
int equalsZero(float (*));
void setZero(float (*));
void setEquals(float (*), float (*));

void polyTimesXn(float (*), int);
void printPoly(float (*), int);

int contains(int(*), int, int);
int done(int (*), int (*), int, int);

void updateFinishedRows(float (*)[M][maxDegree+1], int (*));
void updateFinishedColumns(float (*)[M][maxDegree+1], int (*));

int dividesRowAndCol(float (*)[M][maxDegree+1], int, int);
void eucDiv(float (*), float (*), float (*));

int getRank(float (*)[M][maxDegree+1]);
int min(int, int);

void orderDiagonals(float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);
void orderHelper(float (*)[M][maxDegree+1], int *, int *, int);

void findLeastEntry(float (*)[M][maxDegree+1], int (*), int (*), int *, int *, int *);
