#include <stdio.h>
#include "matrixOperations.h"
#include "polynomialOperations.h"

void initializeA(float A[][M][maxDegree+1], float Acopy[][M][maxDegree+1]) {
	int i, n, m;
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			setZero(A[n][m]);
			setZero(Acopy[n][m]);
		}
	}

	// A[0][0][0] = 1;
	// A[0][0][1] = 0;
	// A[0][0][2] = 1;
	// A[0][1][0] = 0;
	// A[0][1][1] = -4;
	// A[0][1][2] = 3;
	// A[1][0][0] = 0;
	// A[1][0][1] = 8;
	// A[1][0][2] = 0;
	// A[1][1][0] = 0;
	// A[1][1][1] = 0;
	// A[1][1][2] = -1;

	float temp;
	float negative;
	// printf("Initialize the matrix \"A\".\n");
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			for(i = 0; i <= maxEntry; ++i) {
				// printf("Enter the coefficient of x^%d in A[%d][%d]: ", i, n, m);
				// scanf("%f", &temp);

				//after this, negative will be either +1 or -1
				negative = (float) (rand()%2);
				negative = negative - 0.5;
				negative *= 2;

				temp = (float) (rand()%10);
				temp *= negative;

				A[n][m][i] = temp;
				Acopy[n][m][i] = temp;
			}
		}
	}
	
}

void init(float A[][M][maxDegree+1], float Acopy[][M][maxDegree+1], float P[][N][maxDegree+1], 
		  float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1]) {
	initializeA(A, Acopy);
	initializeP(P); // initialized to identity
	initializeP(Pinv);
	initializeQ(Q); // intiialized to identity
	initializeQ(Qinv);
}

void printAll(float A[N][M][maxDegree+1], float Acopy[N][M][maxDegree+1], 
			  float P[N][N][maxDegree+1], float Pinv[N][N][maxDegree+1], 
			  float Q[M][M][maxDegree+1], float Qinv[M][M][maxDegree+1], 
			  float PAtest[N][M][maxDegree+1], float diagTest[M][N][maxDegree+1], 
			  float PPinvTest[N][N][maxDegree+1], float QQinvTest[M][M][maxDegree+1]) {

	printf("diag: ");
	print2ArrayM(A, N);

	// printf("\nP: ");
	// print2ArrayN(P, N);

	// printf("Q: ");
	// print2ArrayM(Q, M);

	// matNNxmatNM(P, Acopy, PAtest);
	// matNMxmatMM(PAtest, Q, diagTest);
	// printf("diagTest: ");
	// print2ArrayM(diagTest, N);

	// matxvecN(P, b, Pb);

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

	// printf("P: ");
	// print2ArrayN(P, N);

	// printf("Q: ");
	// print2ArrayM(Q, M);

	// printf("Pinv: ");
	// print2ArrayN(Pinv, N);
	
	// printf("Qinv: ");
	// print2ArrayM(Qinv, M);
	
}

//initializes P to identity
void initializeP(float P[][N][maxDegree+1]) {
	int i, j;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < N; ++j) {
			setZero(P[i][j]);
			if (i == j) {
				P[i][j][0] = 1;
			}
		}
	}
}

// initializes Q to identity
void initializeQ(float Q[][M][maxDegree+1]) {
	int i, j;
	for(i = 0; i < M; ++i) {
		for(j = 0; j < M; ++j) {
			setZero(Q[i][j]);
			if (i == j) {
				Q[i][j][0] = 1;
			}
		}
	}
}

// calculates the rank of a matrix where there is only one entry per row and one entry per column
int getRank(float A[][M][maxDegree+1]) {
	int rank = 0;
	int n, m;
	for (n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if (equalsZero(A[n][m]) == 0) {
				++rank;
				break;
			}
		}
	}
	return rank;
}

void rowOperations1(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], int n, float unit) {
	type1rowM(A, n, unit);
	type1rowN(P, n, unit);
	type1rowN(Pinv, n, unit);
}

void rowOperations2(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], int n, int tempN, float q[]) {
	type2rowM(A, n, tempN, q, 1);
	type2rowN(P, n, tempN, q, 1);
	type2rowN(Pinv, tempN, n, q, 0);

	int i, j;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < M; ++j) {
			clearZeroes(A[i][j]);
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
	type2col(A, N, m, tempM, q, 1);
	type2col(Q, M, m, tempM, q, 1);
	type2col(Qinv, M, tempM, m, q, 0);
	
	// to handle the floating point errors
	int i, j;
	for(i = 0; i < N; ++i) {
		for(j = 0; j < M; ++j) {
			clearZeroes(A[i][j]);
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
	float temp[maxDegree+1];
	for(m = 0; m < M; ++m) {
		setZero(temp);
		//printf("a = q, b = A[%d][%d], c = temp", addFrom, m);
		multiply(q, A[addFrom][m], temp);
		if (neg == 1) {
			subtract(A[addTo][m], temp);
		}
		else {
			add(A[addTo][m], temp);
		}
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
	float temp[maxDegree+1];
	for(n = 0; n < N; ++n) {
		//printf("a = q, b = A[%d][%d], c = temp", addFrom, n);
		multiply(q, A[addFrom][n], temp);
		if (neg == 1) {
			subtract(A[addTo][n], temp);
		}
		else {
			add(A[addTo][n], temp);
		}
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
	float temp[maxDegree+1];
	for(i = 0; i < len; ++i) {
		//printf("a = q, b = A[%d][%d], c = temp", i, addFrom);
		multiply(q, A[i][addFrom], temp);
		if (neg == 1) {
			subtract(A[i][addTo], temp);
		}
		else {
			add(A[i][addTo], temp);
		}
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

void print2ArrayM(float A[][M][maxDegree+1], int len) {
	int m, n;
	printf("\n");
	for(n = 0; n < len; ++n) {
		for(m = 0; m < M; ++m) {
			if (m == 0 && n == 0) { 
				printf("(("); 
				printPoly(A[n][m], 0);
				printf(", ");
			}
			else if (m == 0) { 
				printf("(");
				printPoly(A[n][m], 0);
				printf(", "); 
			}
			else if (m == M - 1 && n == len - 1) {
				//printf(",");
				printPoly(A[n][m], 0);
				printf("))");
			}
			else if (m == M - 1) {
				printPoly(A[n][m], 0);
				printf(")");
			}
			else {
				printPoly(A[n][m], 0);
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
				printPoly(A[n][m], 0);
				printf(", ");
			}
			else if (m == 0) { 
				printf("(");
				printPoly(A[n][m], 0);
				printf(", "); 
			}
			else if (m == M - 1 && n == len - 1) {
				//printf(",");
				printPoly(A[n][m], 0);
				printf("))");
			}
			else if (m == M - 1) {
				printPoly(A[n][m], 0);
				printf(")");
			}
			else {
				printPoly(A[n][m], 0);
				printf(", ");
			}
		}
		printf("\n");
	}
	printf("\n");
}

void printArray(float A[][maxDegree+1], int len) {
	int i;
	for(i = 0; i < len; ++i) {
		if (i == 0) { 
			printf("(");
			printPoly(A[i], 0);
			printf(", "); 
		}
		else if (i == len - 1) {
			printPoly(A[i], 0);
			printf("))");
		}
		else {
			printPoly(A[i], 0);
			printf(", ");
		}
	}
	printf("\n");
}

void makeAllMonic(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1]) {
	int i;
	int tempDeg;
	float tempUnit;
	for(i = 0; i < rank; ++i) {
		tempDeg = degree(A[i][i]);
		tempUnit = A[i][i][tempDeg];
		rowOperations1(A, P, Pinv, i, 1/tempUnit);
	}
}