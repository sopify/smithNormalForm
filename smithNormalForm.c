#include <stdio.h>

/* STATUS AS OF 3:30 AM, MONDAY APRIL 9:
- dividesRowAndCol produces "Floating point exception"
- both branches of divides calls a division by zero
*/

/* STATUS AS OF 1:00 PM, MONDAY APRIL 9:
- almost works, just doesn't go through on some last steps
*/

/* STATUS AS OF 8:30 PM, MONDAY APRIL 9:
- P, Q and Pinv, Qinv are working properly
- leastEntryAlgo takes A to a form where there is
	at most one non-zero entry in each row and column 
*/

/* STATUS AS OF 12:50 AM, TUESDAY APRIL 10:
- PAQ does not multiply together to get proper diagonal matrix
- matrix is diagonalized and diagonals are ordered from least to greatest
- still need to put matrix in Smith normal form
*/

/* STATUS AS OF 1:08 AM, TUESDAY APRIL 10:
- now PAQ = a diagonal matrix
- these diagonals are ordered from least to greatest
- still need to put matrix in Smith normal form,
	although sometimes (usually) it happens automatically

TO DO:
- write function to multiply matrices and test automatically the following:
	- whether PAQ = diag
	- whether Q*Qinv = Qinv*Q = I_M
	- whether P*Pinv = Pinv*P = I_N
*/


int M, N;

void leastEntryAlgo(int (*)[M], int (*)[N], int (*)[N], int (*)[M], int (*)[M]);

void initializeA(int (*)[M]);
void initializeP(int (*)[N]);
void initializeQ(int (*)[M]);

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

int contains(int(*), int, int);
int done(int (*), int(*), int, int);

void updateFinishedRows(int (*)[M], int (*));
void updateFinishedColumns(int (*)[M], int (*));

int dividesRowAndCol(int (*)[M], int, int);
int eucDiv(int, int);
int min(int, int);

void findLeastEntry(int (*)[M], int (*), int (*), int *, int *, int *);
void findLeastEntry2(int (*)[M], int *, int *, int);

main() {
	M = 2; //initialize M, # columns
	N = 3; //initialize N, # rows
	int A[N][M]; // matrix to diagonalize
	int P[N][N];
	int Pinv[N][N];
	int Q[M][M];
	int Qinv[M][M];
	// A[0][0] = 2;
	// A[0][1] = -1;
	// A[1][0] = 3;
	// A[1][1] = 1;
	// A[2][0] = 1;
	// A[2][1] = 7;

	initializeA(A);
	initializeP(P); // initialized to identity
	initializeP(Pinv);
	initializeQ(Q); // intiialized to identity
	initializeQ(Qinv);

	printf("A: ");
	print2ArrayM(A, N);

	leastEntryAlgo(A, P, Pinv, Q, Qinv);
	
	printf("diag: ");
	print2ArrayM(A, N);

	printf("P: ");
	print2ArrayN(P, N);

	printf("Q: ");
	print2ArrayM(Q, M);

	printf("Pinv: ");
	print2ArrayN(Pinv, N);
	
	printf("Qinv: ");
	print2ArrayM(Qinv, M);
}


void leastEntryAlgo(int A[][M], int P[][N], int Pinv[][N], int Q[][M], int Qinv[][M]) {
	int n, m; // used in for loops
	int tempN, tempM, tempMin; //used to store temporary locations in A
	int finished = 0; // bool to decide if loop is over
	int divides, allEntriesPos, smith;  // 1 is true, 0 is false
	int tempRowEntry, tempColEntry, tempMult, tempEntry;
	int diag = 0;
	int counter = 0;
	
	// for finishedRows and finishedColumns,
	// set entry i to -1 if row/col i is not finished, or
	// set entry i to i if row/col i is finished
	int finishedRows[N];
	int finishedColumns[M];

	//initializes finishedRows and finishedColumns
	for(n = 0; n < N; ++n) {
		finishedRows[n] = -1;
	}
	for(m = 0; m < M; ++m) {
		finishedColumns[m] = -1;
	}

	// printf("finishedRows: ");
	// printArray(finishedRows, N);
	// printf("finishedColumns: ");
	// printArray(finishedColumns, M);

	while (finished == 0) {

		// finds least non-zero entry, stores it in tempN, tempM
		tempM = -1;
		tempN = -1;
		findLeastEntry(A, finishedRows, finishedColumns, &tempN, &tempM, &finished);
		//printf("Least entry: A[%i][%i] = %i\n", tempN, tempM, A[tempN][tempM]);

		if(finished == 1) { break; }
		
		//checks if entry divides its entire row and column
		divides = dividesRowAndCol(A, tempN, tempM);
		//printf("divides: %i\n", divides);
		/* if entry divides its entire row and column,
		we annihilate the entire row and column with euclidean division */
		if (divides == 1) {
			for(n = 0; n < N; ++n) {
				//we do not operate on the row itself
				if (n != tempN) {
					tempMult = eucDiv(A[n][tempM], A[tempN][tempM]);
					//printf("\nMultiplier: %i\n", tempMult);
				/* if tempMult is 0, this operation is useless
				   so we do not execute it.*/
					if (tempMult != 0) {
						rowOperations2(A, P, Pinv, n, tempN, tempMult);
					}
				}
			}
			for(m = 0; m < M; ++m) {
				// we do not operate on the column itself
				if (m != tempM) {
					tempMult = eucDiv(A[tempN][m], A[tempN][tempM]);
					//printf("\nMultiplier: %i\n", tempMult);
					/* if tempMult is 0, this operation is useless
				   	so we do not execute it. */
					if (tempMult != 0) {
						columnOperations2(A, Q, Qinv, m, tempM, tempMult);
					}
				}
			}
		}
		else {
			tempRowEntry = -1;
			tempColEntry = -1;
			if(contains(finishedColumns, M, tempM) == 0) {
				for(n = 0; n < N; ++n) {
					if(A[n][tempM] != 0 && n != tempN) {
						tempRowEntry = n;
						break;
					}
				}
				if (tempRowEntry == -1) { break; }
				tempMult = eucDiv(A[tempRowEntry][tempM], A[tempN][tempM]);
				rowOperations2(A, P, Pinv, tempRowEntry, tempN, tempMult);
				
			}
			else {
				for(m = 0; m < M; ++m) {
					if(A[tempN][m] != 0 && m != tempM) {
						tempColEntry = m;
						break;
					}
				}
				if (tempColEntry == -1) { break; }
				tempMult = eucDiv(A[tempN][tempColEntry], A[tempN][tempM]);
				columnOperations2(A, Q, Qinv, tempColEntry, tempM, tempMult);
			}
		}

		// updates rows and colums that are finished
		updateFinishedRows(A, finishedRows);
		updateFinishedColumns(A, finishedColumns);
		finished = done(finishedRows, finishedColumns, N, M);
		//++counter;
		//printf("finished: %i\n", finished);
		//print2ArrayM(A, N);

	}

	/* if any entries are negative at this point,
	   multiply whole row by -1 */
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if(A[n][m] < 0) {
				rowOperations1(A, P, Pinv, n, -1);
			}
		}
	}

	/* perform type 3 operations to order the 
	   non-zero elements onto the diagonals
	*/
	while(counter < min(N, M)) {
		tempM = -1;
		tempN = -1;

		//printf("\ncounter: %i", counter);
		//print2ArrayM(A, N);

		findLeastEntry2(A, &tempN, &tempM, counter);
		rowOperations3(A, P, Pinv, tempN, counter);
		columnOperations3(A, Q, Qinv, tempM, counter);

		++counter;
	}

	smith = 1;
	for(n = 0; n < min(N, M) - 1; ++n) {
		if (A[n+1][n+1] % A[n][n] != 0) {
			smith = 0;
			break;
		}
	}

	while(smith == 0) {
		//
		smith = 1;
		for(n = 0; n < min(N, M) - 1; ++n) {
			if (A[n+1][n+1] % A[n][n] != 0) {
				//smith = 0;
				break;
			}
		}
	}

	// computation trick: transpose Q and Pinv so they are correct
	transposeM(Q);
	transposeN(Pinv);
}

int min(int a, int b) {
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

/*
to find inverses, we will use the rule that transpose(transpose(invP)) = invP
*/

// row operations

// multiplies entire row by unit for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type1rowM(int A[][M], int row, int unit) {
	if (unit != 1 && unit != -1) {
		printf("A type 1 row operation did not execute because the input was not a unit.");
		return;
	}
	else {
		int m;
		for(m = 0; m < M; ++m) {
			A[row][m] *= unit;
		}	
	}
}

/* row operation of type 2 will take two matrices, A and P,
and transform them the same way.
For parameters: pointers to both matrices,
the row to add to, the row to add from,
and the multiple of the row we are adding
*/


// adds mult*row addFrom to row addTo for column size M
//transpose of this operation is flipping addTo and addFrom
//inverse of this operation is changing mult to -1*mult
void type2rowM(int A[][M], int addTo, int addFrom, int mult) {
	int m;
	for(m = 0; m < M; ++m) {
		A[addTo][m] += mult*A[addFrom][m];
	}
}


// interchanges row1 and row2 for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type3rowM(int A[][M], int row1, int row2) {
	int m;
	for(m = 0; m < M; ++m) {
		A[row1][m] += A[row2][m];
		A[row2][m] = A[row1][m] - A[row2][m];
		A[row1][m] -= A[row2][m];
	}
}

// transpose of this operation is itself
// inverse of this operation is itself
void type1rowN(int A[][N], int row, int unit) {
	if (unit != 1 && unit != -1) {
		printf("A type 1 row operation did not execute because the input was not a unit.");
		return;
	}
	else {
		int n;
		for(n = 0; n < N; ++n) {
			A[row][n] *= unit;
		}	
	}
}

// adds mult*row addFrom to row addTo for column size M
//transpose of this operation is flipping addTo and addFrom
//inverse of this operation is changing mult to -1*mult
void type2rowN(int A[][N], int addTo, int addFrom, int mult) {
	int n;
	for(n = 0; n < N; ++n) {
		A[addTo][n] += mult*A[addFrom][n];
	}
}


// interchanges row1 and row2 for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type3rowN(int A[][N], int row1, int row2) {
	int n;
	for(n = 0; n < N; ++n) {
		A[row1][n] += A[row2][n];
		A[row2][n] = A[row1][n] - A[row2][n];
		A[row1][n] -= A[row2][n];
	}
}

// column operations

//multiplies entire row by unit
void type1col(int A[][M], int len, int col, int unit) {
	if (unit != 1 && unit != -1) {
		printf("A type 1 column operation did not execute because the input was not a unit.");
		return;
	}
	else {
		int i;
		for(i = 0; i < len; ++i) {
			A[i][col] *= unit;
		}	
	}
}

// adds mult*column addFrom to column addTo
void type2col(int A[][M], int len, int addTo, int addFrom, int mult) {
	int i;
	for(i = 0; i < len; ++i) {
		A[i][addTo] += mult*A[i][addFrom];
	}
}

//interchanges col1 and col2
void type3col(int A[][M], int len, int col1, int col2) {
	int i;
	for(i = 0; i < len; ++i) {
		A[i][col1] += A[i][col2];
		A[i][col2] = A[i][col1] - A[i][col2];
		A[i][col1] -= A[i][col2];
	}
}

//takes transpose of a size N*N matrix
void transposeN(int A[][N]) {
	int i, j;
	for(i = 0; i < N; ++i) {
		for(j = i+1; j < N; ++j) {
			A[i][j] += A[j][i];
			A[j][i] = A[i][j] - A[j][i];
			A[i][j] -= A[j][i];
		}
	}
}

void transposeM(int A[][M]) {
	int i, j;
	for(i = 0; i < M; ++i) {
		for(j = i+1; j < M; ++j) {
			A[i][j] += A[j][i];
			A[j][i] = A[i][j] - A[j][i];
			A[i][j] -= A[j][i];
		}
	}
}


void print2ArrayM(int A[][M], int len) {
	int m, n;
	printf("\n");
	for(n = 0; n < len; ++n) {
		for(m = 0; m < M; ++m) {
			if (m == 0 && n == 0) { printf("((%i, ", A[n][m]); }
			else if (m == 0) { printf("(%i, ", A[n][m]); }
			else if (m == M - 1 && n == len - 1) { printf("%i))", A[n][m]); }
			else if (m == M - 1) { printf("%i),", A[n][m]); }
			else { printf("%i, ", A[n][m]); }
		}
		printf("\n");
	}
	printf("\n");
}

void print2ArrayN(int A[][N], int len) {
	int m, n;
	printf("\n");
	for(n = 0; n < len; ++n) {
		for(m = 0; m < N; ++m) {
			if (m == 0 && n == 0) { printf("((%i, ", A[n][m]); }
			else if (m == 0) { printf("(%i, ", A[n][m]); }
			else if (m == N - 1 && n == len - 1) { printf("%i))", A[n][m]); }
			else if (m == N - 1) { printf("%i),", A[n][m]); }
			else { printf("%i, ", A[n][m]); }
		}
		printf("\n");
	}
	printf("\n");
}

void initializeA(int A[][M]) {
	int n, m;
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			A[n][m]=(rand()%100);
		}
	}
}


//initializes P to identity
void initializeP(int P[][N]) {
	int i, j;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < N; ++j) {
			if (i == j) {
				P[i][j] = 1;
			}
			else {
				P[i][j] = 0;
			}
		}
	}
}

// initializes Q to identity
void initializeQ(int Q[][M]) {
	int i, j;
	for(i = 0; i < M; ++i) {
		for(j = 0; j < M; ++j) {
			if (i == j) {
				Q[i][j] = 1;
			}
			else {
				Q[i][j] = 0;
			}
		}
	}
}

//returns absolute value of an integer
int absolu(int x) {
	if (x < 0) {
		x = -1*x;
	}
	return x;
}

// returns 1 if x is in A, 0 if x is not in A
int contains(int A[], int len, int x) {
	int i;
	int ret = 0;
	for(i = 0; i < len; ++i) {
		if (A[i] == x) {
			ret = 1;
			break;
		}
	}
	return ret;
}

void printArray(int A[], int len) {
	int i;
	for(i = 0; i < len; ++i) {
		if(i == 0) { printf(" (%i, ", A[i]); }
		else if(i == len - 1) { printf("%i)\n", A[i]); }
		else { printf("%i, ", A[i]); }
	}
}

int done(int A[], int B[], int lenA, int lenB) {
	int i;
	int ret = 1;
	for(i = 0; i < lenA; ++i) {
		if(contains(A, lenA, i) == 0) {
			ret = 0;
			break;
		}
	}
	if (ret == 1) {
		for(i = 0; i < lenB; ++i) {
			if(contains(B, lenB, i) == 0) {
				ret = 0;
				break;
			}
		}
	}
	return ret;
}

void updateFinishedRows(int A[][M], int finishedRows[]) {
	int tempCount, m, n;
	for (n = 0; n < N; ++n) {
		tempCount = 0;
		for(m = 0; m < M; ++m) {
			if(A[n][m] != 0) { ++tempCount; }			
		}
		if (tempCount < 2) { finishedRows[n] = n; }
		else { finishedRows[n] = -1; }
	}
}

void updateFinishedColumns(int A[][M], int finishedColumns[]) {
	int tempCount, m, n;
	for(m = 0; m < M; ++m) {
		tempCount = 0;
		for(n = 0; n < N; ++n) {
			if(A[n][m] != 0) { ++tempCount; }
		}
		if(tempCount < 2) { finishedColumns[m] = m; }
		else { finishedColumns[m] = -1; }
	}
}

void findLeastEntry(int A[][M], int finishedRows[], int finishedColumns[], int *tempN, int *tempM, int *finished) {
	int m, n;
	int boole = 0;
	*finished = 0;
	int tempMin = -1;
	
	// printf("finishedRows:");
	// printArray(finishedRows, N);
	// printf("finishedColumns: ");
	// printArray(finishedColumns, M);

	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if(!(contains(finishedRows, N, n) == 1 && contains(finishedColumns, M, m) == 1)) {
				if(A[n][m] != 0) {
					*tempN = n;
					*tempM = m;
					tempMin = absolu(A[n][m]);
					boole = 1;
					break;
				}
				if (boole == 1) {
					break;
				}
			}
		}
	}
	if (tempMin == -1) {
		printf("There is no valid least entry\n");
		*finished = 1;
		return;
	}
	for(n = 0; n < N; ++n) {
		for (m = 0; m < M; ++m) {
			if(contains(finishedRows, N, n) == 0 || contains(finishedColumns, M, m) == 0) {
				if (A[n][m] != 0) {
					if(absolu(A[n][m]) < tempMin) {
						*tempN = n;
						*tempM = m;
						tempMin = absolu(A[n][m]);
					}
				}
			}
		}
	}
	if (*tempM == -1 || *tempN == -1) { *finished = 1; }
}

int dividesRowAndCol(int A[][M], int tempN, int tempM) {
	int m, n;
	int ret = 1;
	if (A[tempN][tempM] == 0) {
		printf("You have called this with A[i][j] = 0\n");
		return 0;
	}
	if (tempN == -1 || tempM == -1) {
		printf("You have called this with tempN = -1 or tempM = -1\n");
		return 0;
	}
	for(n = 0; n < N; ++n) {
		if (n != tempN) {
			if (A[n][tempM] % A[tempN][tempM] != 0) {
				ret = 0;
			}
		}
	}

	for(m = 0; m < M; ++m) {
		if (m != tempM) {
			if (A[tempN][m] % A[tempN][tempM] != 0) {
				ret = 0;
			}
		}
	}
	return ret;
}

// returns the 'q' in a = bq + r
int eucDiv(int a, int b) {
	if (a == 0) {
		return 0;
	}
	if (b == 0) {
		printf("You are calling a division by 0\n");
		return a;
	}
	int temp = a % b;
	if (temp < 0) { temp += absolu(b); }
	return (a - temp) / b;
}

void findLeastEntry2(int A[][M], int *tempN, int *tempM, int counter) {
	int m, n;
	int tempMin = -1;
	int boole = 0;

	for(n = counter; n < N; ++n) {
		for(m = counter; m < M; ++m) {
			if(A[n][m] != 0) {
				*tempN = n;
				*tempM = m;
				tempMin = absolu(A[n][m]);
				boole = 1;
				break;
			}
			if (boole == 1) {
				break;
			}
		}
	}
	if (tempMin == -1) {
		printf("This should not execute.\n");
		return;
	}

	for(n = counter; n < N; ++n) {
		for (m = counter; m < M; ++m) {
			if (A[n][m] != 0) {
				if(A[n][m] < tempMin) {
					*tempN = n;
					*tempM = m;
					tempMin = A[n][m];
				}
			}
		}
	}

	if(*tempN == -1 || *tempM == -1) { printf("tempN or tempN is equal to -1\n");}
}

void rowOperations1(int A[][M], int P[][N], int Pinv[][N], int n, int unit) {
	type1rowM(A, n, unit);
	type1rowN(P, n, unit);
	type1rowN(Pinv, n, unit);
}

void rowOperations2(int A[][M], int P[][N], int Pinv[][N], int n, int tempN, int tempMult) {
	type2rowM(A, n, tempN, -1*tempMult);
	type2rowN(P, n, tempN, -1*tempMult);
	type2rowN(Pinv, tempN, n, tempMult);
}

void columnOperations2(int A[][M], int Q[][M], int Qinv[][M], int m, int tempM, int tempMult) {
	type2col(A, N, m, tempM, -1*tempMult);
	type2col(Q, M, tempM, m, -1*tempMult);
	type2col(Qinv, M, m, tempM, tempMult);
}

void rowOperations3(int A[][M], int P[][N], int Pinv[][N], int row1, int row2) {
	if (row1 != row2) {
		// row2 < min(N, M) && row1 < min(N, M) && 
		type3rowM(A, row1, row2);
		type3rowN(P, row1, row2);
		type3rowN(Pinv, row1, row2);
	}
}

void columnOperations3(int A[][M], int Q[][M], int Qinv[][M], int col1, int col2) {
	if (col1 != col2) {
		//col2 < min(N, M) && col1 < min(N, M) &&
		type3col(A, N, col1, col2);
		type3col(Q, M, col1, col2);
		type3col(Qinv, M, col1, col2);
	}
}