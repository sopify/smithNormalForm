int M, N, rank, consistent, numEDs;
int maxED = 100;

void leastEntryAlgo(int (*)[M], int (*)[N], int (*)[N], int (*)[M], int (*)[M]);

void initializeSize(int *, char (*));

void initializeA(int (*)[M], int (*)[M]);
void initializeP(int (*)[N]);
void initializeQ(int (*)[M]);
void initializeb(int (*));

void matxvecN(int (*)[N], int (*), int (*));
void matxvecM(int (*)[M], int (*), int (*));
void matxvecNM(int (*)[M], int (*), int (*));

void matNNxmatNN(int (*)[N], int (*)[N], int (*)[N]);
void matMMxmatMM(int (*)[M], int (*)[M], int (*)[M]);

void matNNxmatNM(int (*)[N], int (*)[M], int (*)[M]);
void matNMxmatMM(int (*)[M], int (*)[M], int (*)[M]);

void calcy(int (*)[N], int (*), int (*));

void type1rowM(int (*)[M], int, int);
void type2rowM(int (*)[M], int, int, int);
void type3rowM(int (*)[M], int, int);

void type1rowN(int (*)[N], int, int);
void type2rowN(int (*)[N], int, int, int);
void type3rowN(int (*)[N], int, int);

void type1col(int (*)[M], int, int, int);
void type2col(int (*)[M], int, int, int, int);
void type3col(int (*)[M], int, int, int);

void rowOperations1(int (*)[M], int (*)[N], int (*)[N], int, int);

void columnOperations2(int (*)[M], int (*)[M], int (*)[M], int, int, int);
void rowOperations2(int (*)[M], int (*)[N], int (*)[N], int, int, int);

void rowOperations3(int (*)[M], int (*)[N], int (*)[N], int, int);
void columnOperations3(int (*)[M], int (*)[M], int (*)[M], int, int);

void transposeN(int (*)[N]);
void transposeM(int (*)[M]);

void print2ArrayM(int (*)[M], int);
void print2ArrayN(int (*)[N], int);
void printArray(int (*), int);
void printEDs(int (*), int);

int contains(int(*), int, int);
int done(int (*), int(*), int, int);
void makeAllDiagsPositive(int (*)[M], int (*)[N], int (*)[N]);

void updateFinishedRows(int (*)[M], int (*));
void updateFinishedColumns(int (*)[M], int (*));

int dividesRowAndCol(int (*)[M], int, int);
int eucDiv(int, int);
int min(int, int);
int comp(const void* a, const void* b);

int checkSmith(int (*)[M]);
int getRank(int (*)[M]);

void smithTransform(int (*)[N], int (*)[N], int (*)[N], int(*)[M], int(*)[M], int);
void orderDiagonals(int (*)[N], int (*)[N], int (*)[N], int(*)[M], int(*)[M]);

void findLeastEntry(int (*)[M], int (*), int (*), int *, int *, int *);
void findLeastEntry2(int (*)[M], int *, int *, int);

void primeFactorize(int, int (*));
void getDiags(int (*)[M], int (*));
void getElementaryDivisors(int (*), int (*));
void initializeZero(int (*), int);