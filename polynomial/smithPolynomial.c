#include <stdio.h>
#include "smithPolynomial.h"

/*
@author: Christian Drappi
*/

/*	TO DO:
	- given a matrix of differential operators, 
	  find the constraints for a vector of 
	  arbitrary functions that this matrix acts on
*/


/*	State as of Apr 28, 2:08am:
	- make all polynomials monic
	- diagonalizes
*/


main() {
	initializeSize(&N, "row");
	initializeSize(&M, "column");
	changeMaxDegree(&maxDegree, &maxEntry);

	float A[N][M][maxDegree+1]; // matrix to diagonalize
	float Acopy[N][M][maxDegree+1];
	float P[N][N][maxDegree+1];
	float Pinv[N][N][maxDegree+1];
	float Q[M][M][maxDegree+1];
	float Qinv[M][M][maxDegree+1];
	float PAtest[N][M][maxDegree+1];
	float diagTest[M][N][maxDegree+1];
	float PPinvTest[N][N][maxDegree+1];
	float QQinvTest[M][M][maxDegree+1];
	float b[N][maxDegree+1];

	
	init(A, Acopy, P, Pinv, Q, Qinv, PAtest, diagTest, PPinvTest, QQinvTest, b);
	// initializeA(A, Acopy);
	// initializeP(P); // initialized to identity
	// initializeP(Pinv);
	// initializeQ(Q); // intiialized to identity
	// initializeQ(Qinv);

	printf("\nA: ");
	print2ArrayM(A, N);

	leastEntryAlgo(A, P, Pinv, Q, Qinv);

	// printf("diag: ");
	// print2ArrayM(A, N);

	printAll(A, Acopy, P, Pinv, Q, Qinv, PAtest, diagTest, PPinvTest, QQinvTest, b);
}

void leastEntryAlgo(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1]) {
	int n, m; // used in for loops over rows, columns respectively
	int p; // used in for loops over polynomial array
	int tempN, tempM, tempMin; //used to store temporary locations in A
	int finished = 0; // bool to decide if loop is over
	int tempRowEntry, tempColEntry, tempMult, tempEntry;
	int diag = 0;
	float q[maxDegree+1];
	
	// for finishedRows and finishedColumns,
	// set entry i to -1 if row/col i is not finished, or
	// set entry i to i if row/col i is finished
	int finishedRows[N];
	int finishedColumns[M];


	int maxN;
	int maxM;
	int maxD;
	int tempD;
	//initializes finishedRows and finishedColumns
	for(n = 0; n < N; ++n) {
		finishedRows[n] = -1;
	}
	for(m = 0; m < M; ++m) {
		finishedColumns[m] = -1;
	}
	// int c = 0;
	while (finished == 0) {

		// finds least non-zero entry, stores it in tempN, tempM
		tempM = -1;
		tempN = -1;
		findLeastEntry(A, finishedRows, finishedColumns, &tempN, &tempM, &finished);
		// printf("\nLeast entry: A[%i][%i] = ", tempN, tempM);
		// printPoly(A[tempN][tempM], 1);

		if(finished == 1) { break; }
		
		tempRowEntry = -1;
		tempColEntry = -1;
		if(contains(finishedColumns, M, tempM) == 0) {
			for(n = 0; n < N; ++n) {
				if(equalsZero(A[n][tempM]) == 0 && n != tempN) {
					tempRowEntry = n;
					break;
				}
			}
			if (tempRowEntry == -1) { break; }
			
			// printf("attempts eucdiv on rows: A[%d][%d] = ", tempRowEntry, tempM);
			// printPoly(A[tempRowEntry][tempM], 1);
			// printf("\nA: ");
			// print2ArrayM(A, N);
			
			eucDiv(A[tempRowEntry][tempM], A[tempN][tempM], q);
			// printf("q_final: ");
			// printPoly(q, 1);
			// printf("operate on rows, A[%d][%d] = ", tempRowEntry, tempM);
			// printPoly(A[tempRowEntry][tempM], 1);

			maxD = -1;
			for(n = 0; n < N; ++n) {
				for(m = 0; m < M; ++m) {
					tempD = degree(A[n][m]);
					if(tempD > maxD) {
						maxD = tempD;
						maxN = n;
						maxM = m;
					}
				}
			}
			// printf("Max degree pre operation: %d, A[%d][%d] = ", maxD, maxN, maxM);
			// printPoly(A[maxN][maxM], 1);

			rowOperations2(A, P, Pinv, tempRowEntry, tempN, q);

			// printf("after eucdiv on rows: A[%d][%d] = ", tempRowEntry, tempM);
			// printPoly(A[tempRowEntry][tempM], 1);

			// printf("\nA: ");
			// print2ArrayM(A, N);
		}
		else {
			for(m = 0; m < M; ++m) {
				if(equalsZero(A[tempN][m]) == 0 && m != tempM) {
					tempColEntry = m;
					break;
				}
			}
			if (tempColEntry == -1) { break; }
			
			// printf("attempts eucdiv on cols: A[%d][%d] = ", tempN, tempColEntry);
			// printPoly(A[tempN][tempColEntry], 1);
			
			eucDiv(A[tempN][tempColEntry], A[tempN][tempM], q);
			// printf("q_final: ");
			// printPoly(q, 1);
			// printf("operate on cols A[%d][%d] = ", tempN, tempColEntry);
			// printPoly(A[tempN][tempColEntry], 1);

			maxD = -1;
			for(n = 0; n < N; ++n) {
				for(m = 0; m < M; ++m) {
					tempD = degree(A[n][m]);
					if(tempD > maxD) {
						maxD = tempD;
						maxN = n;
						maxM = m;
					}
				}
			}
			// printf("Max degree pre operation: %d, A[%d][%d] = ", maxD, maxN, maxM);
			// printPoly(A[maxN][maxM], 1);

			columnOperations2(A, Q, Qinv, tempColEntry, tempM, q);

			
			// printf("after eucdiv on cols: A[%d][%d] = ", tempN, tempColEntry);
			// printPoly(A[tempN][tempColEntry], 1);

			// printf("\nA: ");
			// print2ArrayM(A, N);
		}
		// printf("\nA: ");
		// print2ArrayM(A, N);


		maxD = -1;
		for(n = 0; n < N; ++n) {
			for(m = 0; m < M; ++m) {
				tempD = degree(A[n][m]);
				if(tempD > maxD) {
					maxD = tempD;
					maxN = n;
					maxM = m;
				}
			}
		}
		// printf("Max degree post operation: %d, A[%d][%d] = ", maxD, maxN, maxM);
		// printPoly(A[maxN][maxM], 1);

		// updates rows and colums that are finished
		updateFinishedRows(A, finishedRows);
		updateFinishedColumns(A, finishedColumns);
		finished = done(finishedRows, finishedColumns, N, M);

		//printf("finished: %i\n", finished);
		//print2ArrayM(A, N);
		// ++c;

	}

	rank = getRank(A);

	// /* perform type 3 operations to order the 
	//    non-zero elements onto the diagonals
	// */
	orderDiagonals(A, P, Pinv, Q, Qinv);
	
	/* make all polynomials monic */
	makeAllMonic(A, P, Pinv);


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
	//makeAllMonic(A, P, Pinv);

	// computation trick: transpose Q and Pinv so they are correct
	transposeM(Qinv);
	transposeN(Pinv);
}

void initializeSize(int *size, char type[]) {
	int num;
	printf("Please enter the amount of %ss in your matrix: ", type);
	scanf("%d", &num);
	*size = num;
}

void changeMaxDegree(int *maxDeg, int *maxEnt) {
	int num;
	printf("Enter the maximum degree in your matrix of polynomials: ");
	scanf("%d", &num);
	while (num < 0) {
		printf("Try again:\n");
		scanf("%d", &num);
	}
	*maxEnt = num;
	*maxDeg = N*M*(num + 1);
}

// void matxvecDiff(float P[][N][maxDegree+1], float b[][maxDegree+1], float Pb[][][maxDegree+1]) {
// 	int i, j;
// 	float temp[maxDegree+1];
// 	float temp2[maxDegree+1];
// 	for(i = 0; i < N; ++i) {
// 		for(j = 0; j < N; ++j) {
// 			add(Pb[i][j], P[i][j]);
// 		}
// 	}
// }

void matNNxmatNM(float P[][N][maxDegree+1], float A[][M][maxDegree+1], float PA[][M][maxDegree+1]) {
	int n, m, n1;
	float temp[maxDegree+1];
	float temp2[maxDegree+1];
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			setZero(temp);
			for(n1 = 0; n1 < N; ++n1) {
				multiply(P[n][n1], A[n1][m], temp2);
				add(temp, temp2);
			}
			setEquals(PA[n][m], temp);
		}
	}
}

void matNMxmatMM(float A[][M][maxDegree+1], float Q[][M][maxDegree+1], float AQ[][M][maxDegree+1]) {
	int n, m, m1;
	float temp[maxDegree+1];
	float temp2[maxDegree+1];
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			setZero(temp);
			for(m1 = 0; m1 < M; ++m1) {
				multiply(A[n][m1], Q[m1][m], temp2);
				add(temp, temp2);
			}
			setEquals(AQ[n][m], temp);
		}
	}
}

void matNNxmatNN(float A[][N][maxDegree+1], float B[][N][maxDegree+1], float AB[][N][maxDegree+1]) {
	int n1, n2, n;
	float temp[maxDegree+1];
	float temp2[maxDegree+1];
	for(n1 = 0; n1 < N; ++n1) {
		for(n2 = 0; n2 < N; ++n2) {
			setZero(temp);
			for(n = 0; n < N; ++n) {
				multiply(A[n1][n], B[n][n2], temp2);
				add(temp, temp2);
			}
			setEquals(AB[n1][n2], temp);
		}
	}
}


void matMMxmatMM(float A[][M][maxDegree+1], float B[][M][maxDegree+1], float AB[][M][maxDegree+1]) {
	int m1, m2, m;
	float temp[maxDegree+1];
	float temp2[maxDegree+1];
	for(m1 = 0; m1 < M; ++m1) {
		for(m2 = 0; m2 < M; ++m2) {
			setZero(temp);
			for(m = 0; m < M; ++m) {
				multiply(A[m1][m], B[m][m2], temp2);
				add(temp, temp2);
			}
			setEquals(AB[m1][m2], temp);
		}
	}
}


/****************************MATRIX_OPERATIONS_C***********************************/

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

				temp = (float) (rand()%5);
				temp *= negative;

				A[n][m][i] = temp;
				Acopy[n][m][i] = temp;
			}
		}
	}
	
}

void init(float A[][M][maxDegree+1], float Acopy[][M][maxDegree+1], float P[][N][maxDegree+1], 
		  float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1], 
		  float PAtest[N][M][maxDegree+1], float diagTest[M][N][maxDegree+1], 
		  float PPinvTest[N][N][maxDegree+1], float QQinvTest[M][M][maxDegree+1], 
		  float b[N][maxDegree+1]) {

	initializeA(A, Acopy);
	initializeP(P); // initialized to identity
	initializeP(Pinv);
	initializeQ(Q); // intiialized to identity
	initializeQ(Qinv);

	int m, n;
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			setZero(diagTest[n][m]);
			setZero(PAtest[n][m]);
		}
	}


	initializeP(PPinvTest);
	initializeQ(QQinvTest);
	
	for(n = 0; n < N; ++n) {
		setZero(b[n]);
	}

}

void printAll(float A[N][M][maxDegree+1], float Acopy[N][M][maxDegree+1], 
			  float P[N][N][maxDegree+1], float Pinv[N][N][maxDegree+1], 
			  float Q[M][M][maxDegree+1], float Qinv[M][M][maxDegree+1], 
			  float PAtest[N][M][maxDegree+1], float diagTest[M][N][maxDegree+1], 
			  float PPinvTest[N][N][maxDegree+1], float QQinvTest[M][M][maxDegree+1], 
			  float b[N][maxDegree+1]) {

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

	matNNxmatNM(P, Acopy, PAtest);
	matNMxmatMM(PAtest, Q, diagTest);
	printf("diagTest: ");
	print2ArrayM(diagTest, N);

	printf("Pb = \n");
	printDiff(P, 0);
	printf("\n");

	matNNxmatNN(P, Pinv, PPinvTest);
	int n1, n2;
	for(n1 = 0; n1 < N; ++n1) {
		for(n2 = 0; n2 < N; ++n2) {
			clearZeroes2(PPinvTest[n1][n2]);
		}
	}
	printf("PPinvTest: ");
	print2ArrayM(PPinvTest, N);

	matNNxmatNN(Q, Qinv, QQinvTest);
	int m1, m2;
	for(m1 = 0; m1 < M; ++m1) {
		for(m2 = 0; m2 < M; ++m2) {
			clearZeroes2(QQinvTest[m1][m2]);
		}
	}
	printf("QQinvTest: ");
	print2ArrayM(QQinvTest, M);

	if (N > M) {
		printf("Given this matrix, the following are conditions that must be satisfied so this system of differential equations is consistent:\n");
		printDiff(P, M);
	}
	
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

	// int i, j;
	// for(i = 0; i < N; ++i) {
	// 	for(j = 0; j < M; ++j) {
	// 		clearZeroes(A[i][j]);
	// 	}
	// }
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
	// int i, j;
	// for(i = 0; i < N; ++i) {
	// 	for(j = 0; j < M; ++j) {
	// 		clearZeroes(A[i][j]);
	// 	}
	// }
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

void printDiff(float P[][N][maxDegree+1], int start) {
	int n;
	for(n = start; n < N; ++n) {
		printDiffRow(P[n]);
		if (start == M) {
			printf(" = 0");
		}
		printf("\n");
	}
}

void printDiffRow(float p[][maxDegree+1]) {
	int n;
	printf("(");
	for(n = 0; n < N; ++n) {
		if (equalsZero(p[n]) == 0) {	
			printf("(");
			printPoly(p[n], 0);
			printf(")");
			printf("* f_%d", n);
			if (n + 1 < N) {
				printf(" + ");
			}
		}
	}
	printf(")");
}

void makeAllMonic(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1]) {
	int i, j;
	int tempDeg;
	float tempUnit;
	for(i = 0; i < rank; ++i) {
		for(j = 0; j <= maxDegree; ++j) {
			if(j <= maxEntry*min(N, M)) {
				if (A[i][i][j] < precision && A[i][i][j] > -precision) {
					A[i][i][j] = 0.0;
				}
			}
			else {
				if (A[i][i][j] < precision2 && A[i][i][j] > -precision2) {
					A[i][i][j] = 0.0;
				}
			}
		}
		tempDeg = degree(A[i][i]);
		tempUnit = A[i][i][tempDeg];
		rowOperations1(A, P, Pinv, i, 1/tempUnit);
	}
}


/****************************POLYNOMIAL_OPERATIONS_C***********************************/

void clearZeroes(float x[]) {
	int i;
	for(i = 0; i <= maxDegree; ++i) {
		if (x[i] < precision && x[i] > -precision) {
			x[i] = 0.0;
		}
	}
}

void clearZeroes2(float x[]) {
	int i;
	for(i = 0; i <= maxDegree; ++i) {
		if (x[i] < precision2 && x[i] > -precision2 && i != 0) {
			x[i] = 0.0;
		}
	}
}

void add(float addTo[], float addFrom[]) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		addTo[i] += addFrom[i];
	}
	clearZeroes(addTo);
}

void subtract(float subTo[], float subFrom[]) {
	int i, j;
	for(i = 0; i <= maxDegree; ++i) {
		subTo[i] -= subFrom[i];
	}

	clearZeroes(subTo);
	// printf("tempA: ");
	// printPoly(subTo);

	// printf("b scaled: ");
	// printPoly(subFrom);
}

// multiplies 2 polynomials, overwriting their product as the first argument
// this assumes that deg(a*b) <= maxDegree
void multiply(float a[], float b[], float c[]) {
	int i, j;
	setZero(c);
	// int degA = degree(a);
	// int degB = degree(b);
	// int d = degA + degB;
	for (i = maxDegree; i >= 0; --i) {
		for (j = i; j >= 0; --j) {
			c[i] += a[j] * b[i-j];
			//printf("j = %i, i = %i, a[j] = %f, b[i-j] = %f\n", j, i, a[j], b[i-j]);
		}
	}
	//clearZeroes(c);
}

int degree(float x[]) {
	int i;
	int ret = 0;
	for (i = maxDegree; i >= 0; --i) {
		if (x[i] > precision || x[i] < -1*precision) {
			ret = i;
			break;
		}
	}
	return ret;
}

// takes a vector a to b*a
void scale(float a[], float b) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		a[i] *= b;
	}
}

int equalsZero(float poly[]) {
	int ret = 1;
	int i;
	clearZeroes(poly);
	for (i = 0; i <= maxDegree; ++i) {
		if (poly[i] != 0) {
			ret = 0;
			break;
		}
	}
	return ret;
}

void setZero(float x[]) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		x[i] = 0.0;
	}
}

//sets arg1 equal to arg2
void setEquals(float x[], float y[]) {
	int i;
	for (i = 0; i <= maxDegree; ++i) {
		x[i] = y[i];
	}
}

void polyTimesXn(float a[], int power) {
	
	if (degree(a) + power > maxDegree) {
		printf("deg: %d, power: %d\n", degree(a), power);
		printf("This would create a polynomial with degree greater than the maximum degree\n");
		return;
	}
	int i;
	if (power == 0) {
		return;
	}
	else if (power > 0) {
		for (i = maxDegree; i >= 0; --i) {
			if (i >= power) {
				a[i] = a[i-power];
			}
			else {
				a[i] = 0;
			}
		}
	}
	else {
		for (i = 0; i <= maxDegree; ++i) {
			if (i - power <= maxDegree) {
				a[i] = a[i - power];
			}
			else {
				a[i] = 0;
			}
		}
	}
}

void printPoly(float poly[], int line) {
	int i;
	for(i = 0; i <= maxDegree; ++i) {
		if (i == 0) {
			if (degree(poly) > 0) {
				if (poly[i] != 0) {
					printf("%1.3f", poly[i]);
					printf(" + ");
				}
			}
			else {
				printf("%1.3f", poly[i]);
			}
		}
		else if (i == degree(poly)) {
			printf("%1.3f*D^%d", poly[i], i);
		}
		else {
			if(poly[i] > precision || poly[i] < -precision) {
				printf("%1.3f*D^%d + ", poly[i], i);
			}
		}
	}
	if (line == 1) {
		printf("\n");
	}
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
			if(contains(finishedRows, N, n) == 0 || contains(finishedColumns, M, m) == 0) {
				if(equalsZero(A[n][m]) == 0) {
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
					if(degree(A[n][m]) < tempMin) {
						*tempN = n;
						*tempM = m;
						tempMin = degree(A[n][m]);
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
	int i;
	// creates tempA as a clone of a
	float tempA[maxDegree+1];
	setEquals(tempA, a);
	float tempB[maxDegree+1];
	setEquals(tempB, b);	
	float q[maxDegree+1];
	setZero(q);

	int degA = degree(a);
	int degB = degree(b);
	int d = degA - degB;
	
	float tempQ;
	
	// algorithm for euclidean division
	for (i = d; i >= 0; --i) {
		tempQ = tempA[degA + i - d] / tempB[degB];
		q[i] = tempQ;
		if (tempQ != 0) {
			scale(tempB, tempQ);
			polyTimesXn(tempB, i);
			subtract(tempA, tempB);
			polyTimesXn(tempB, -i);
			scale(tempB, 1.0/tempQ);
		}
	}
	multiply(b, q, tempA);
	subtract(tempA, a);
	if(equalsZero(tempA) == 1) {
		return 1;
	}
	else {
		return 0;
	}
}

//the input for c will be tempA from eucdiv
// void ldMult(float b[], float q[], float c[]) {
// 	int i, j;
// 	for (i = maxDegree; i >= 0; --i) {
// 		for (j = i; j >= 0; --j) {
// 			c[i] -= b[j] * q[i-j];
// 		}
// 	}

// 	//eliminate inevitable floating point precision related errors
// 	clearZeroes(c);

// 	// printf("b: ");
// 	// printPoly(b);
// 	// printf("q: ");
// 	// printPoly(q);
// 	// printf("tempA: ");
// 	// printPoly(c);
// }

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

	// printf("a in eucdiv call: ");
	// printPoly(a);

	// printf("b in eucdiv call: ");
	// printPoly(b);
	int i;
	// printf("\narg1 of eucdiv: ");
	// printPoly(a);
	// printf("\narg2 of eucdiv: ");
	// printPoly(b);
	int degA = degree(a);
	int degB = degree(b);
	//printf("\ndeg arg1: %d, deg arg2: %d\n", degA, degB);

	int d = degA - degB;
	if (d < 0) { 
		printf("This is a nonsensical call of eucDiv\n");
		return;
	}

	float tempQ;
	setZero(q);

	// creates tempA as a clone of a
	float tempA[maxDegree+1];
	setEquals(tempA, a);
	float tempB[maxDegree+1];
	setEquals(tempB, b);
	// float tempPoly[maxDegree+1];

	// algorithm for euclidean division
	for (i = d; i >= 0; --i) {
		tempQ = tempA[degA + i - d] / tempB[degB];
		q[i] = tempQ;
		if (tempQ != 0) {
			scale(tempB, tempQ);
			polyTimesXn(tempB, i);
			subtract(tempA, tempB);
			polyTimesXn(tempB, -i);
			scale(tempB, 1.0/tempQ);
		}
	}
}

void orderHelper(float A[][M][maxDegree+1], int *tempN, int *tempM, int counter) {
	int m, n;
	int tempMin = -1;
	int boole = 0;

	for(n = counter; n < N; ++n) {
		for(m = counter; m < M; ++m) {
			if(equalsZero(A[n][m]) == 0) {
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
		printf("This should not execute Error in orderHelper.\n");
		return;
	}

	for(n = counter; n < N; ++n) {
		for (m = counter; m < M; ++m) {
			if (equalsZero(A[n][m]) == 0) {
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


void orderDiagonals(float A[][M][maxDegree+1], float P[][N][maxDegree+1], float Pinv[][N][maxDegree+1], float Q[][M][maxDegree+1], float Qinv[][M][maxDegree+1]) {
	int counter = 0;
	int tempM, tempN;
	while(counter < rank) {
		tempM = -1;
		tempN = -1;
		orderHelper(A, &tempN, &tempM, counter);
		if (tempN != counter || tempM != counter) {
			rowOperations3(A, P, Pinv, tempN, counter);
			columnOperations3(A, Q, Qinv, tempM, counter);
		}
		++counter;
	}
}

int min(int a, int b) {
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}