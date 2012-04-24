#include <stdio.h>

/*
@author: Christian Drappi
*/

/*	TO DO:
	- diagonalize matrix of polynomials
*/


int M, N, rank, consistent;
int maxDegree;

void leastEntryAlgo(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);

void initializeSize(int *, char (*));
void changeMaxDegree(int *);

void initializeA(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);
void initializeP(float (*)[N][maxDegree+1]);
void initializeQ(float (*)[M][maxDegree+1]);
void initializeb(float (*));

// void matxvecN(float (*)[N][maxDegree+1], float (*), float (*));
// void matxvecM(float (*)[M][maxDegree+1], float (*), float (*));
// void matxvecNM(float (*)[M][maxDegree+1], float (*), float (*));

// void matNNxmatNN(float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1]);
// void matMMxmatMM(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);

// void matNNxmatNM(float (*)[N][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);
// void matNMxmatMM(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1]);

// void calcy(float (*)[N][maxDegree+1], float (*), float (*));

void type1rowM(float (*)[M][maxDegree+1], int, float);
void type2rowM(float (*)[M][maxDegree+1], int, int, float (*), int);
void type3rowM(float (*)[M][maxDegree+1], int, int);

void type1rowN(float (*)[N][maxDegree+1], int, float);
void type2rowN(float (*)[N][maxDegree+1], int, int, float (*), int);
void type3rowN(float (*)[N][maxDegree+1], int, int);

void type1col(float (*)[M][maxDegree+1], int, int, float);
void type2col(float (*)[M][maxDegree+1], int, int, int, float(*), int);
void type3col(float (*)[M][maxDegree+1], int, int, int);

void rowOperations1(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int, float);

void columnOperations2(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], int, int, float (*));
void rowOperations2(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int, int, float(*));

void rowOperations3(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int, int);
void columnOperations3(float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], float (*)[M][maxDegree+1], int, int);

void transposeN(float (*)[N][maxDegree+1]);
void transposeM(float (*)[M][maxDegree+1]);

void print2ArrayM(float (*)[M][maxDegree+1], int);
void print2ArrayN(float (*)[N][maxDegree+1], int);
void printArray(float (*)[maxDegree+1], int);

int contains(int(*), int, int);
int done(int (*), int (*), int, int);
// void makeAllDiagsPositive(float (*)[M][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1]);

void updateFinishedRows(float (*)[M][maxDegree+1], int (*));
void updateFinishedColumns(float (*)[M][maxDegree+1], int (*));

int dividesRowAndCol(float (*)[M][maxDegree+1], int, int);
void eucDiv(float (*), float (*), float (*));
int min(int, int);

// int checkSmith(float (*)[M][maxDegree+1]);
int getRank(float (*)[M][maxDegree+1]);

// void smithTransform(float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int(*)[M][maxDegree+1], int(*)[M][maxDegree+1], int);
// void orderDiagonals(float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], float (*)[N][maxDegree+1], int(*)[M][maxDegree+1], int(*)[M][maxDegree+1]);

void findLeastEntry(float (*)[M][maxDegree+1], int (*), int (*), int *, int *, int *);
void findLeastEntry2(float (*)[M][maxDegree+1], int *, int *, int);

int equalsZero(float (*));
void setEquals(float (*), float (*));

void add(float (*), float (*));
void scale(float (*), float);
void polyTimesXn(float (*), int);
void multiply(float (*), float (*));

main() {
	// initializeSize(&N, "row");
	// initializeSize(&M, "column");
	
	// changeMaxDegree(&maxDegree);
	N = 3;
	M = 3;
	maxDegree = 20;

	float A[N][M][maxDegree+1]; // matrix to diagonalize
	float P[N][N][maxDegree+1];
	float Pinv[N][N][maxDegree+1];
	float Q[M][M][maxDegree+1];
	float Qinv[M][M][maxDegree+1];
	float b[N][maxDegree+1];
	float x[M][maxDegree+1];
	float y[M][maxDegree+1];
	// int Pb[N][maxDegree+1];
	// int testB[N][maxDegree+1];
	float Acopy[N][M][maxDegree+1];
	float PAtest[N][M][maxDegree+1];
	float diagTest[M][N][maxDegree+1];
	float PPinvTest[N][N][maxDegree+1];
	float QQinvTest[M][M][maxDegree+1];

	A[0][0][0] = 4;
	A[0][0][1] = 2;
	A[0][0][2] = 3;
	A[0][1][0] = 8;
	A[0][1][1] = 4;
	A[0][1][2] = 1;
	A[0][2][0] = 9;
	A[0][2][1] = 0;
	A[0][2][2] = 4;
	A[1][0][0] = 3;
	A[1][0][1] = 5;
	A[1][0][2] = 6;
	A[1][1][0] = 7;
	A[1][1][1] = 8;
	A[1][1][2] = 9;
	A[1][2][0] = 3;
	A[1][2][1] = 2;
	A[1][2][2] = 0;
	A[2][0][0] = 4;
	A[2][0][1] = 6;
	A[2][0][2] = 8;
	A[2][1][0] = 9;
	A[2][1][1] = 3;
	A[2][1][2] = 7;
	A[2][2][0] = 3;
	A[2][2][1] = 8;
	A[2][2][2] = 0;

	//initializeA(A, Acopy);
	// initializeb(b);
	initializeP(P); // initialized to identity
	initializeP(Pinv);
	initializeQ(Q); // intiialized to identity
	initializeQ(Qinv);

	printf("\nA: ");
	print2ArrayM(A, N);

	leastEntryAlgo(A, P, Pinv, Q, Qinv);
	// matxvecN(P, b, Pb);

	printf("diag: ");
	print2ArrayM(A, N);

	printf("\nP: ");
	print2ArrayN(P, N);

	printf("Q: ");
	print2ArrayM(Q, M);

	// matNNxmatNM(P, Acopy, PAtest);
	// matNMxmatMM(PAtest, Q, diagTest);
	// printf("diagTest: ");
	// print2ArrayM(diagTest, N);


	// printf("Pb = ");
	// printArray(Pb, N);

	// calcy(A, Pb, y);
	// printf("y = ");
	// printArray(y, M);

	// matNNxmatNN(P, Pinv, PPinvTest);
	// printf("PPinvTest: ");
	// print2ArrayM(PPinvTest, N);

	// matNNxmatNN(Q, Qinv, QQinvTest);
	// printf("QQinvTest: ");
	// print2ArrayM(QQinvTest, M);

	// if (consistent == 1) {
	// 	matxvecM(Q, y, x);

	// 	printf("Solution: x = ");
	// 	printArray(x, M);

	// 	matxvecNM(Acopy, x, testB);
		
	// 	printf("Testing solution: b = ");
	// 	printArray(b, N);
	// 	printf("Testing solution: testB = ");
	// 	printArray(testB, N);
	// }

	
	// printf("diag: ");
	// print2ArrayM(A, N);

	// printf("P: ");
	// print2ArrayN(P, N);

	// printf("Q: ");
	// print2ArrayM(Q, M);

	// printf("Pinv: ");
	// print2ArrayN(Pinv, N);
	
	// printf("Qinv: ");
	// print2ArrayM(Qinv, M);
}

// adds 
void add(float addTo[], float addFrom[]) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		addTo[i] += addFrom[i];
	}
}

// takes a vector a to b*a
void scale(float a[], float b) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		a[i] *= b;
	}	
}

void polyTimesXn(float a[], int power) {
	
	if (degree(a) + power > maxDegree) {
		printf("This would create a polynomial with degree greater than the maximum degree\n");
		return;
	}
	int i;
	if (power >= 0) {
		for (i = maxDegree - 1; i >= 0; --i) {
			if (i >= power) {
				a[i] = a[i-power];
			}
			else {
				a[i] == 0;
			}
		}
	}
	else {
		for (i = 0; i <= maxDegree; ++i) {
			if (i <= maxDegree - power) {
				a[i] = a[i+power];
			}
			else {
				a[i] == 0;
			}
		}
	}
}

// multiplies 2 polynomials, assuming that deg(a*b) <= maxDegree
void multiply(float multTo[], float multFrom[]) {
	int i, j, k;
	for (i = maxDegree; i >= 0; --i) {
		for (j = i; j >= 0; --j) {
			multTo[i] += multTo[j] * multFrom[i-j];
		}
	}
}

int equalsZero(float poly[]) {
	int ret = 1;
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		if (poly[i] < 0.001) {
			ret = 0;
			break;
		}
	}
	return ret;
}

void initializeSize(int *size, char type[]) {
	int num;
	printf("Please enter the amount of %ss in your matrix: ", type);
	scanf("%d", &num);
	*size = num;
}

void changeMaxDegree(int *maxDeg) {
	int num;
	printf("Enter the maximum degree in your matrix of polynomials.\n");
	scanf("%d", &num);
	while (num < 0) {
		printf("Try again:\n");
		scanf("%d", &num);
	}
	*maxDeg = N*M*num;
}

void initializeA(float A[][M][maxDegree+1], float Acopy[][M][maxDegree+1]) {
	int i, n, m;
	float temp;
	printf("Initialize the matrix \"A\".\n");
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			for(i = 0; i <= maxDegree/(N*M); ++i) {
				printf("The coefficient of x^%d in A[%d][%d] = :", i, n, m);
				scanf("%f", &temp);
				A[n][m][i] = temp;
				Acopy[n][m][i] = temp;
			}
		}
	}
}

// void initializeb(int b[][maxDegree+1]) {
// 	printf("Initialize the vector \"b\".\n");
// 	int n, temp;
// 	for(n = 0; n < N; ++n) {
// 		printf("b[%d] = ", n);
// 		scanf("%d", &temp);
// 		b[n] = temp;
// 	}
// }

// void calcy(float A[][M][maxDegree+1], int vec[][maxDegree+1], int y[][maxDegree+1]) {
// 	int m, n = 0;
// 	consistent = 1;
// 	for(m = 0; m < M; ++m) {
// 		if (m >= rank) {
// 			y[m] = 0;
// 		}
// 		else {
// 			if (vec[m] % A[m][m] == 0) {
// 				y[m] = vec[m] / A[m][m];
// 			}
// 			else {
// 				printf("This system is not consistent\n");
// 				consistent = 0;
// 				break;
// 			}
// 		}
// 	}
// }

void leastEntryAlgo(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1]) {
	int n, m; // used in for loops over rows, columns respectively
	int p; // used in for loops over polynomial array
	int tempN, tempM, tempMin; //used to store temporary locations in A
	int finished = 0; // bool to decide if loop is over
	int divides, allEntriesPos, smith;  // 1 is true, 0 is false
	int tempRowEntry, tempColEntry, tempMult, tempEntry;
	int diag = 0;
	float q[maxDegree];
	
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

	while (finished == 0) {

		// finds least non-zero entry, stores it in tempN, tempM
		tempM = -1;
		tempN = -1;
		findLeastEntry(A, finishedRows, finishedColumns, &tempN, &tempM, &finished);
		//printf("Least entry: A[%i][%i] = %i\n", tempN, tempM, A[tempN][tempM][maxDegree+1]);

		if(finished == 1) { break; }
		
		//checks if entry divides its entire row and column
		divides = dividesRowAndCol(A, tempN, tempM);

		/* if entry divides its entire row and column,
		we annihilate the entire row and column with euclidean division */
		if (divides == 1) {
			for(n = 0; n < N; ++n) {
				//we do not operate on the row itself
				if (n != tempN) {
					eucDiv(A[n][tempM], A[tempN][tempM], q);

				/* if tempMult is 0, this operation is useless
				   so we do not execute it.*/
					if (tempMult != 0) {
						rowOperations2(A, P, Pinv, n, tempN, q);
					}
				}
			}
			for(m = 0; m < M; ++m) {
				// we do not operate on the column itself
				if (m != tempM) {
					eucDiv(A[tempN][m], A[tempN][tempM], q);
					//printf("\nMultiplier: %i\n", tempMult);
					/* if tempMult is 0, this operation is useless
				   	so we do not execute it. */
					if (tempMult != 0) {
						columnOperations2(A, Q, Qinv, m, tempM, q);
					}
				}
			}
			printf("\nA: ");
			print2ArrayM(A, N);
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
				eucDiv(A[tempRowEntry][tempM], A[tempN][tempM], q);
				rowOperations2(A, P, Pinv, tempRowEntry, tempN, q);
				
			}
			else {
				for(m = 0; m < M; ++m) {
					if(A[tempN][m] != 0 && m != tempM) {
						tempColEntry = m;
						break;
					}
				}
				if (tempColEntry == -1) { break; }
				eucDiv(A[tempN][tempColEntry], A[tempN][tempM], q);
				columnOperations2(A, Q, Qinv, tempColEntry, tempM, q);
			}
			printf("\nA: ");
			print2ArrayM(A, N);
		}

		// updates rows and colums that are finished
		updateFinishedRows(A, finishedRows);
		updateFinishedColumns(A, finishedColumns);
		finished = done(finishedRows, finishedColumns, N, M);

		//printf("finished: %i\n", finished);
		//print2ArrayM(A, N);

	}

	//rank = getRank(A);

	/* if any entries are negative at this point,
	   multiply whole row by -1 */
	// makeAllDiagsPositive(A, P, Pinv);

	// /* perform type 3 operations to order the 
	//    non-zero elements onto the diagonals
	// */
	// orderDiagonals(A, P, Pinv, Q, Qinv);

	// smith = checkSmith(A);
	// while(smith != -1) {
	// 	/* write something that takes a_n, a_n+1 and transforms
	// 	   it into gcd(a_n, a_n+1) and lcm(a_n, a_n+1) */
	// 	smithTransform(A, P, Pinv, Q, Qinv, smith);
	// 	// this call ensures all entries are ordered
	// 	orderDiagonals(A, P, Pinv, Q, Qinv);
	// 	smith = checkSmith(A);
	// }
	// print2ArrayM(A, N);
	// // once again, ensure all diagonal elements are positive
	// makeAllDiagsPositive(A, P, Pinv);

	// computation trick: transpose Q and Pinv so they are correct
	transposeM(Qinv);
	transposeN(Pinv);
}

//initializes P to identity
void initializeP(float P[][N][maxDegree+1]) {
	int i, j;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < N; ++j) {
			if (i == j) {
				P[i][j][0] = 1;
			}
			// else {
			// 	P[i][j] = 0;
			// }
		}
	}
}

// initializes Q to identity
void initializeQ(float Q[][M][maxDegree+1]) {
	int i, j;
	for(i = 0; i < M; ++i) {
		for(j = 0; j < M; ++j) {
			if (i == j) {
				Q[i][j][0] = 1;
			}
			// else {
			// 	Q[i][j] = 0;
			// }
		}
	}
}

//returns absolute value of an integer
int absolu(int x) {
	if (x < 0) {
		x = -1.0*x;
	}
	return x;
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
void type1rowM(float A[][M][maxDegree+1], int row, float unit) {
	int m;
	for(m = 0; m < M; ++m) {
		scale(A[row][m], unit);
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
void type2rowM(float A[][M][maxDegree+1], int addTo, int addFrom, float q[], int neg) {
	int m;
	float temp[maxDegree];
	for(m = 0; m < M; ++m) {
		setEquals(temp, q);
		if (neg == 1) { scale(temp, -1.0); }
		multiply(temp, A[addFrom][m]);
		add(A[addTo][m], temp);
	}
}


// interchanges row1 and row2 for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type3rowM(float A[][M][maxDegree+1], int row1, int row2) {
	int m, j;
	for(m = 0; m < M; ++m) {
		for(j = 0; j <= maxDegree; ++j) {
			A[row1][m][j] += A[row2][m][j];
			A[row2][m][j] = A[row1][m][j] - A[row2][m][j];
			A[row1][m][j] -= A[row2][m][j];
		}
	}
}

// transpose of this operation is itself
// inverse of this operation is itself
void type1rowN(float A[][N][maxDegree+1], int row, float unit) {
	int n;
	for(n = 0; n < N; ++n) {
		scale(A[row][n],unit);
	}	
}

// adds mult*row addFrom to row addTo for column size M
//transpose of this operation is flipping addTo and addFrom
//inverse of this operation is changing mult to -1*mult
void type2rowN(float A[][N][maxDegree+1], int addTo, int addFrom, float q[], int neg) {
	int n;
	float temp[maxDegree];
	for(n = 0; n < N; ++n) {
		setEquals(temp, q);
		if (neg == 1) { scale(temp, -1.0); }
		multiply(temp, A[addFrom][n]);
		add(A[addTo][n], temp);
	}
}


// interchanges row1 and row2 for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type3rowN(float A[][N][maxDegree+1], int row1, int row2) {
	int n, k;
	for(n = 0; n < N; ++n) {
		for(k = 0; k <= maxDegree; ++k) {
			A[row1][n][k] += A[row2][n][k];
			A[row2][n][k] = A[row1][n][k] - A[row2][n][k];
			A[row1][n][k] -= A[row2][n][k];
		}
	}
}

// column operations

//multiplies entire row by unit
void type1col(float A[][M][maxDegree+1], int len, int col, float unit) {
	int i;
	for(i = 0; i < len; ++i) {
		scale(A[i][col], unit);	
	}
}

// adds mult*column addFrom to column addTo
void type2col(float A[][M][maxDegree+1], int len, int addTo, int addFrom, float q[], int neg) {
	int i;
	float temp[maxDegree];
	for(i = 0; i < len; ++i) {
		setEquals(temp, q);
		if (neg == 1) { scale(temp, -1.0); }
		multiply(temp, A[i][addFrom]);
		add(A[i][addTo], temp);
	}
}

//interchanges col1 and col2
void type3col(float A[][M][maxDegree+1], int len, int col1, int col2) {
	int i, j;
	for(i = 0; i < len; ++i) {
		for(j = 0; j <= maxDegree; ++j) {
			A[i][col1][j] += A[i][col2][j];
			A[i][col2][j] = A[i][col1][j] - A[i][col2][j];
			A[i][col1][j] -= A[i][col2][j];
		}
	}
}

//takes transpose of a size N*N matrix
void transposeN(float A[][N][maxDegree+1]) {
	int i, j, k;
	for(i = 0; i < N; ++i) {
		for(j = i+1; j < N; ++j) {
			for(k = 0; k <= maxDegree; ++k) {
				A[i][j][k] += A[j][i][k];
				A[j][i][k] = A[i][j][k] - A[j][i][k];
				A[i][j][k] -= A[j][i][k];
			}
		}
	}
}

void transposeM(float A[][M][maxDegree+1]) {
	int i, j, k;
	for(i = 0; i < M; ++i) {
		for(j = i+1; j < M; ++j) {
			for(k = 0; k <= maxDegree; ++k) {
				A[i][j][k] += A[j][i][k];
				A[j][i][k] = A[i][j][k] - A[j][i][k];
				A[i][j][k] -= A[j][i][k];
			}
		}
	}
}

void printPoly(float poly[]) {
	int i;
	for(i = 0; i <= maxDegree; ++i) {
		if (i == 0) {
			printf("%f", poly[i]);
			if (degree(poly) > 0) {
				printf(" + ");
			}
		}
		else if (i == degree(poly)) {
			printf("%10.2e*x^%d", poly[i], i);
		}
		else {
			if(poly[i] != 0.0) {
				printf("%10.1f*x^%d + ", poly[i], i);
			}
		}
	}
}

void print2ArrayM(float A[][M][maxDegree+1], int len) {
	int m, n;
	printf("\n");
	for(n = 0; n < len; ++n) {
		for(m = 0; m < M; ++m) {
			if (m == 0 && n == 0) { 
				printf("(("); 
				printPoly(A[n][m]);
				printf(", ");
			}
			else if (m == 0) { 
				printf("(");
				printPoly(A[n][m]);
				printf(", "); 
			}
			else if (m == M - 1 && n == len - 1) {
				//printf(",");
				printPoly(A[n][m]);
				printf("))");
			}
			else if (m == M - 1) {
				printPoly(A[n][m]);
				printf(")");
			}
			else {
				printPoly(A[n][m]);
				printf(", ");
			}
		}
		printf("\n");
	}
	printf("\n");
}

void print2ArrayN(float A[][N][maxDegree+1], int len) {
	int m, n;
	printf("\n");
	for(n = 0; n < len; ++n) {
		for(m = 0; m < N; ++m) {
			if (m == 0 && n == 0) { 
				printf("(("); 
				printPoly(A[n][m]);
				printf(", ");
			}
			else if (m == 0) { 
				printf("(");
				printPoly(A[n][m]);
				printf(", "); 
			}
			else if (m == M - 1 && n == len - 1) {
				//printf(",");
				printPoly(A[n][m]);
				printf("))");
			}
			else if (m == M - 1) {
				printPoly(A[n][m]);
				printf(")");
			}
			else {
				printPoly(A[n][m]);
				printf(", ");
			}
		}
		printf("\n");
	}
	printf("\n");
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

void printArray(float A[][maxDegree+1], int len) {
	int i;
	for(i = 0; i < len; ++i) {
		if (i == 0) { 
			printf("(");
			printPoly(A[i]);
			printf(", "); 
		}
		else if (i == len - 1) {
			printPoly(A[i]);
			printf("))");
		}
		else {
			printPoly(A[i]);
			printf(", ");
		}
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

void updateFinishedRows(float A[][M][maxDegree+1], int finishedRows[]) {
	int tempCount, m, n;
	for (n = 0; n < N; ++n) {
		tempCount = 0;
		for(m = 0; m < M; ++m) {
			if(equalsZero(A[n][m]) == 0) { ++tempCount; }			
		}
		if (tempCount < 2) { finishedRows[n] = n; }
		else { finishedRows[n] = -1; }
	}
}

void updateFinishedColumns(float A[][M][maxDegree+1], int finishedColumns[]) {
	int tempCount, m, n;
	for(m = 0; m < M; ++m) {
		tempCount = 0;
		for(n = 0; n < N; ++n) {
			if(equalsZero(A[n][m]) == 0) { ++tempCount; }
		}
		if(tempCount < 2) { finishedColumns[m] = m; }
		else { finishedColumns[m] = -1; }
	}
}

void findLeastEntry(float A[][M][maxDegree+1], int finishedRows[], int finishedColumns[], int *tempN, int *tempM, int *finished) {
	int m, n;
	int boole = 0;
	*finished = 0;
	int tempMin = -1; // stores a minimum degree

	// find one element that does not equal zero
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if(!(contains(finishedRows, N, n) == 1 && contains(finishedColumns, M, m) == 1)) {
				if(equalsZero(A[n][m]) == 0) {
					*tempN = n;
					*tempM = m;
					tempMin = absolu(degree(A[n][m]));
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
				if (equalsZero(A[n][m]) == 0) {
					if(absolu(degree(A[n][m])) < tempMin) {
						*tempN = n;
						*tempM = m;
						tempMin = absolu(degree(A[n][m]));
					}
				}
			}
		}
	}
	if (*tempM == -1 || *tempN == -1) { *finished = 1; }
}

int dividesRowAndCol(float A[][M][maxDegree+1], int tempN, int tempM) {
	int m, n;
	int ret = 1;
	if (equalsZero(A[tempN][tempM]) == 1) {
		printf("You have called this with A[i][j] = the zero vector\n");
		return 0;
	}
	if (tempN == -1 || tempM == -1) {
		printf("You have called this with tempN = -1 or tempM = -1\n");
		return 0;
	}
	for(n = 0; n < N; ++n) {
		if (n != tempN) {
			if (dividesPoly(A[n][tempM], A[tempN][tempM]) == 0) {
				ret = 0;
			}
		}
	}

	for(m = 0; m < M; ++m) {
		if (m != tempM) {
			if (dividesPoly(A[tempN][m], A[tempN][tempM]) == 0) {
				ret = 0;
			}
		}
	}
	return ret;
}

int dividesPoly(float a[], float b[]) {
	float tempA[maxDegree];
	int i;
	int degA = degree(a);
	int degB = degree(b);
	float tempScale;
	setEquals(tempA, a);

	int d = degA - degB;
	if (d < 0) { 
		printf("This is a nonsensical call of dividesPoly\n");
		return;
	}

	for (i = d; i <= 0; --i) {
		tempScale = -1.0*tempA[degA + i - d]/b[degB + i - d];
		scale(b, tempScale);
		polyTimesXn(b, i);
		
		add(tempA, b);

		//change b back to normal
		polyTimesXn(b, -1*i);
		scale(b, 1.0/tempScale);
		// perform these operations, at end of loop, check equalsZero(tempA)   
	}

	return equalsZero(tempA);
}

// sets q equal to the polynomial s.t. a = bq + r
void eucDiv(float a[], float b[], float q[]) {
	if (equalsZero(a) == 1) {
		printf("You are calling eucDiv with a equal to the zero vector\n");
		return;
	}
	if (equalsZero(b) == 1) {
		printf("You are calling a division by 0\n");
		return;
	}

	float tempA[maxDegree];
	int i;
	int degA = degree(a);
	int degB = degree(b);
	float tempScale;
	for (i = 0; i <= maxDegree; ++i) {
		tempA[i] = a[i];
	}

	int d = degA - degB;
	if (d < 0) { 
		printf("This is a nonsensical call of eucDiv\n");
		return;
	}

	for (i = d; i <= 0; --i) {
		tempScale = -1.0*tempA[degA + i - d]/b[degB + i - d];
		scale(b, tempScale);
		polyTimesXn(b, i);
		
		q[i] = -1.0*tempScale;
		add(tempA, b);

		//change b back to normal
		polyTimesXn(b, -1*i);
		scale(b, 1.0/tempScale);
		// perform these operations, at end of loop, check equalsZero(tempA)   
	}
}

void setZero(float x[]) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		x[i] = 0.0;
	}
}

void setEquals(float x[], float y[]) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		x[i] = y[i];
	}
}

int degree(float x[]) {
	int i;
	int ret = 0;
	for (i = 0; i <= maxDegree; ++i) {
		if (x[i] != 0.0) {
			ret = i;
		}
	}
	return ret;
}

void findLeastEntry2(float A[][M][maxDegree+1], int *tempN, int *tempM, int counter) {
	int m, n;
	int tempMin = -1;
	int boole = 0;
	// if (counter >= rank) {
	// 	return;
	// }
	// else {
		for(n = counter; n < N; ++n) {
			for(m = counter; m < M; ++m) {
				if(degree(A[n][m]) != 0) {
					*tempN = n;
					*tempM = m;
					tempMin = degree(A[n][m]);
					boole = 1;
					break;
				}
				if (boole == 1) {
					break;
				}
			}
		}
		if (tempMin == -1) {
			printf("This should not execute Error in findLeastEntry2.\n");
			return;
		}

		for(n = counter; n < N; ++n) {
			for (m = counter; m < M; ++m) {
				if (degree(A[n][m]) != 0) {
					if(degree(A[n][m]) < tempMin) {
						*tempN = n;
						*tempM = m;
						tempMin = degree(A[n][m]);
					}
				}
			}
		}

		if(*tempN == -1 || *tempM == -1) { printf("tempN or tempN is equal to -1\n"); }
	// }
}

void rowOperations1(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], int n, float unit) {
	type1rowM(A, n, unit);
	type1rowN(P, n, unit);
	type1rowN(Pinv, n, unit);
}

void rowOperations2(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], int n, int tempN, float q[]) {
	type2rowM(A, n, tempN, q, 0);
	type2rowN(P, n, tempN, q, 0);
	type2rowN(Pinv, tempN, n, q, 1);

	// to handle the floating point errors
	int i, j, k;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < M; ++j) {
			for(k = 0; k <= maxDegree; ++k) {
				if (A[i][j][k] < .001){
					A[i][j][k] = 0;
				}
			}
		}
	}
}

void rowOperations3(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], int row1, int row2) {
	if (row1 != row2) {
		// row2 < min(N, M) && row1 < min(N, M) && 
		type3rowM(A, row1, row2);
		type3rowN(P, row1, row2);
		type3rowN(Pinv, row1, row2);
	}
}

void columnOperations2(float A[][M][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1], int m, int tempM, float q[]) {
	type2col(A, N, m, tempM, q, 0);
	type2col(Q, M, m, tempM, q, 0);
	type2col(Qinv, M, tempM, m, q, 1);
	// to handle the floating point errors
	int i, j, k;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < M; ++j) {
			for(k = 0; k <= maxDegree; ++k) {
				if (A[i][j][k] < .001) {
					A[i][j][k] = 0;
				}
			}
		}
	}
}

void columnOperations3(float A[][M][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1], int col1, int col2) {
	if (col1 != col2) {
		//col2 < min(N, M) && col1 < min(N, M) &&
		type3col(A, N, col1, col2);
		type3col(Q, M, col1, col2);
		type3col(Qinv, M, col1, col2);
	}
}

// void smithTransform(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1], int pos) {
// 	int tempMult;
// 	float q[maxDegree];
// 	columnOperations2(A, Q, Qinv, pos, pos + 1, -1);

// 	while (A[pos+1][pos] != 0) {
// 		eucDiv(A[pos+1][pos], A[pos][pos], q);
		
// 		// runs Euclidean division
// 		rowOperations2(A, P, Pinv, pos+1, pos, q);
		
// 		// makes A[pos][pos] entry larger than A[pos][pos+1]
// 		if (A[pos+1][pos] != 0) { 
// 			rowOperations3(A, P, Pinv, pos, pos+1);
// 		}
// 	}
// 	eucDiv(A[pos][pos+1], A[pos][pos], q);
// 	columnOperations2(A, Q, Qinv, pos+1, pos, q);
// 	if (equalsZero(A[pos][pos+1]) == 0) {
// 		printf("smithTransform failed. This should never happen");
// 	}
// }

// void makeAllDiagsPositive(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1]) {
// 	int n, m;
// 	for(n = 0; n < N; ++n) {
// 		for(m = 0; m < M; ++m) {
// 			if(A[n][m] < 0) {
// 				rowOperations1(A, P, Pinv, n, -1);
// 			}
// 		}
// 	}
// }

void orderDiagonals(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1]) {
	int counter = 0;
	int tempM, tempN;
	while(counter < rank) {
		tempM = -1;
		tempN = -1;
		findLeastEntry2(A, &tempN, &tempM, counter);
		if (tempN != counter || tempM != counter) {
			rowOperations3(A, P, Pinv, tempN, counter);
			columnOperations3(A, Q, Qinv, tempM, counter);
		}
		++counter;
	}
}

//returns first n where A[n+1][n+1] % A[n][n] != 0. 0 if in smith form
// int checkSmith(float A[][M][maxDegree+1]) {
// 	int smith = -1;
// 	int n;
// 	for(n = 0; n < rank - 1; ++n) {
// 		if (A[n+1][n+1] % A[n][n] != 0) {
// 			smith = n;
// 			break;
// 		}
// 	}
// 	return smith;
// }

// calculates the rank of a diagonal matrix
int getRank(float A[][M][maxDegree+1]) {
	int rank = 0;
	int n, m;
	for (n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if (equalsZero(A[n][m]) != 1) {
				++rank;
			}
		}
	}
	return rank;
}

// void matxvecN(float P[][N][maxDegree+1], int b[][maxDegree+1], int vec[][maxDegree+1]) {
// 	int i, j, temp;
// 	for(i = 0; i < N; ++i) {
// 		temp = 0;
// 		for(j = 0; j < N; ++j) {
// 			temp += P[i][j]*b[j];
// 		}
// 		vec[i] = temp;

// 	}
// }

// // void matxvecM(float Q[][M][maxDegree+1], int y[][maxDegree+1], int vec[][maxDegree+1]) {
// // 	int i, j, temp;
// // 	for(i = 0; i < M; ++i) {
// // 		temp = 0;
// // 		for(j = 0; j < M; ++j) {
// // 			temp += Q[i][j]*y[j];
// // 		}
// // 		vec[i] = temp;
// // 	}
// // }

// // void matxvecNM(float A[][M][maxDegree+1], int x[][maxDegree+1], int vec[][maxDegree+1]) {
// // 	int n, m, temp;
// // 	for(n = 0; n < N; ++n) {
// // 		temp = 0;
// // 		for(m = 0; m < M; ++m) {
// // 			temp += A[n][m]*x[m];
// // 		}
// // 		vec[n] = temp;
// // 	}
// // }

// // // void matNNxmatNM(float P[][N][maxDegree+1], float A[][M][maxDegree+1], int PA[][M][maxDegree+1]) {
// // // 	int n, m, n1, temp;
// // // 	for(n = 0; n < N; ++n) {
// // // 		for(m = 0; m < M; ++m) {
// // // 			temp = 0;
// // // 			for(n1 = 0; n1 < N; ++n1) {
// // // 				temp += P[n][n1]*A[n1][m];
// // // 			}
// // // 			PA[n][m] = temp;
// // // 		}
// // // 	}
// // // }

// // // // void matNMxmatMM(float A[][M][maxDegree+1], float Q[][M][maxDegree+1], float AQ[][M][maxDegree+1]) {
// // // // 	int n, m, m1, temp;
// // // // 	for(n = 0; n < N; ++n) {
// // // // 		for(m = 0; m < M; ++m) {
// // // // 			temp = 0;
// // // // 			for(m1 = 0; m1 < M; ++m1) {
// // // // 				temp += A[n][m1]*Q[m1][m];
// // // // 			}
// // // // 			AQ[n][m] = temp;
// // // // 		}
// // // // 	}
// // // // }

// // // // void matNNxmatNN(float A[][N][maxDegree+1], float B[][N][maxDegree+1], float AB[][N][maxDegree+1]) {
// // // // 	int n1, n2, n, temp;
// // // // 	for(n1 = 0; n1 < N; ++n1) {
// // // // 		for(n2 = 0; n2 < N; ++n2) {
// // // // 			temp = 0;
// // // // 			for(n = 0; n < N; ++n) {
// // // // 				temp += A[n1][n]*B[n][n2];
// // // // 			}
// // // // 			AB[n1][n2] = temp;
// // // // 		}
// // // // 	}
// // // // }


// // // // void matMMxmatMM(float A[][M][maxDegree+1], float B[][M][maxDegree+1], float AB[][M][maxDegree+1]) {
// // // // 	int m1, m2, m, temp;
// // // // 	for(m1 = 0; m1 < M; ++m1) {
// // // // 		for(m2 = 0; m2 < M; ++m2) {
// // // // 			temp = 0;
// // // // 			for(m = 0; m < M; ++m) {
// // // // 				temp += A[m1][m]*B[m][m2];
// // // // 			}
// // // // 			AB[m1][m2] = temp;
// // // // 		}
// // // // 	}
// // // // }