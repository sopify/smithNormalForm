#include <stdio.h>
#include <stdlib.h>
#include "smithInteger.h"

/*
@author: Christian Drappi

This is my first time writing anything in C, so excuse me if I break conventions or
have organized my code poorly.

*/

main(int argc, char* argv[]) {
	// M = 3; //initialize M, # columns
	// N = 3; //initialize N, # rows
	initializeSize(&N, "row");
	initializeSize(&M, "column");

	int A[N][M]; // matrix to diagonalize
	int P[N][N];
	int Pinv[N][N];
	int Q[M][M];
	int Qinv[M][M];
	int b[N];
	int x[M];
	int y[M];
	int Pb[N];
	int testB[N];
	int Acopy[N][M];
	int PAtest[N][N];
	int diagTest[M][N];
	int PPinvTest[N][N];
	int QQinvTest[M][M];

	init (A, Acopy, P, Pinv, Q, Qinv, b);

	printf("\nA: ");
	print2ArrayM(A, N);

	leastEntryAlgo(A, P, Pinv, Q, Qinv);

	int diags[rank];
	int elementaryDivisors[maxED];
	
	printAll(A, Acopy, P, Pinv, Q, Qinv, PAtest, diagTest, PPinvTest, QQinvTest, b, testB, Pb, elementaryDivisors, y, x, diags);
}

void init(int A[][M], int Acopy[][M], int P[][N], int Pinv[][N], int Q[][M], int Qinv[][M], int b[]) {
	initializeA(A, Acopy);
	initializeb(b);
	initializeP(P); // initialized to identity
	initializeP(Pinv);
	initializeQ(Q); // intiialized to identity
	initializeQ(Qinv);
}

void printAll(int A[][M], int Acopy[][M], int P[][N], int Pinv[][N], 
			  int Q[][M], int Qinv[][M], int PAtest[][M], int diagTest[][N], 
			  int PPinvTest[][N], int QQinvTest[][M], int b[], int testB [], 
			  int Pb[], int elementaryDivisors[], int y[], int x[], int diags[]) {

	initializeZero(elementaryDivisors, maxED);
	
	getDiags(A, diags);
	getElementaryDivisors(diags, elementaryDivisors);
	qsort(elementaryDivisors,maxED,sizeof(int),comp);

	// print diag, P, Q
	printf("diag: ");
	print2ArrayM(A, N);

	printf("\nP: ");
	print2ArrayN(P, N);

	printf("Q: ");
	print2ArrayM(Q, M);

	// if you want y or Pb, you can uncomment these
	// printf("Pb = ");
	// printArray(Pb, N);

	// printf("y = ");
	// printArray(y, M);

	// test whether the matrix was correctly diagonalized
	matNNxmatNM(P, Acopy, PAtest);
	matNMxmatMM(PAtest, Q, diagTest);
	matNNxmatNN(P, Pinv, PPinvTest);
	matNNxmatNN(Q, Qinv, QQinvTest);

	// if(reveal == 1) {
	// 	printf("diagTest: ");
	// 	print2ArrayM(diagTest, N);

	// 	printf("PPinvTest: ");
	// 	print2ArrayM(PPinvTest, N);

	// 	printf("QQinvTest: ");
	// 	print2ArrayM(QQinvTest, M);
	// }
	if (choiceB == 1) {
		matxvecN(P, b, Pb);
		calcy(A, Pb, y);
		if (consistent == 1) {
			matxvecM(Q, y, x);

			printf("Solution: x = ");
			printArray(x, M);

			matxvecNM(Acopy, x, testB);
			
			printf("Testing solution: b = ");
			printArray(b, N);
			printf("Testing solution: testB = ");
			printArray(testB, N);
		}
	}

	printEDs(elementaryDivisors, maxED);
}

void leastEntryAlgo(int A[][M], int P[][N], int Pinv[][N], int Q[][M], int Qinv[][M]) {
	int n, m; // used in for loops
	int tempN, tempM, tempMin; //used to store temporary locations in A
	int finished = 0; // bool to decide if loop is over
	int divides, allEntriesPos, smith;  // 1 is true, 0 is false
	int tempRowEntry, tempColEntry, tempMult, tempEntry;
	int diag = 0;
	
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

		//printf("finished: %i\n", finished);
		//print2ArrayM(A, N);

	}

	rank = getRank(A);

	/* if any entries are negative at this point,
	   multiply whole row by -1 */
	makeAllDiagsPositive(A, P, Pinv);

	/* perform type 3 operations to order the 
	   non-zero elements onto the diagonals
	*/
	orderDiagonals(A, P, Pinv, Q, Qinv);

	smith = checkSmith(A);
	while(smith != -1) {
		/* write something that takes a_n, a_n+1 and transforms
		   it into gcd(a_n, a_n+1) and lcm(a_n, a_n+1) */
		smithTransform(A, P, Pinv, Q, Qinv, smith);
		// this call ensures all entries are ordered
		orderDiagonals(A, P, Pinv, Q, Qinv);
		smith = checkSmith(A);
	}
	// once again, ensure all diagonal elements are positive
	makeAllDiagsPositive(A, P, Pinv);

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

void initializeZero(int arr[], int len) {
	int i;
	for(i = 0; i < len; ++i) {
		arr[i] = 0;
	}
}

void initializeA(int A[][M], int Acopy[][M]) {
	int n, m, temp, negative;
	int choiceA1;
	printf("If you would like to manually enter a matrix of integers, press 0.\nIf you would like have one generated randomly, press 1: ");
	scanf("%d", &choiceA1);
	while (choiceA1 != 0 && choiceA1 != 1) {
		printf("\ntry again: ");
		scanf("%d", &choiceA1);
	}

	choiceA = choiceA1;

	if (choiceA == 1) {
		for(n = 0; n < N; ++n) {
			for(m = 0; m < M; ++m) {
				// this makes negative -1 or +1
				negative = (rand()%2) + 1;
				negative *= 2;
				negative -= 3;

				// temp will take on an integer on random in [-10, 10]
				temp = (rand()%10);
				temp *= negative;

				A[n][m] = temp;
				Acopy[n][m] = temp;
			}
		}
	}

	else {
		printf("Please initialize the matrix \"A\".\n");
		for(n = 0; n < N; ++n) {
			for(m = 0; m < M; ++m) {
				printf("A[%d][%d] = ", n, m);
				scanf("%d", &temp);
				A[n][m] = temp;
				Acopy[n][m] = temp;
			}
		}
	}
}

// puts all the diagonal elements into an array
void getDiags(int A[][M], int diags[]) {
	int i;
	for(i = 0; i < rank; ++i) {
		diags[i] = A[i][i];
	}
}

// gets the elementary divisors using Prof. Zhou's algo
void getElementaryDivisors(int diags[], int EDs[]) {
	int i, j;
	int ret;
	for (i = 0; i < rank; ++i) {
		primeFactorize(diags[i], EDs);
	}
}

// prompts the user to enter a "b"
void initializeb(int b[]) {
	int n;
	if (choiceA == 0) {
		int choiceB1;
		printf("Would you like to enter a \"b\" to solve Ax = b?\nIf so, press 1. If not, press 0: ");
		scanf("%d", &choiceB1);
		
		while (choiceB1 != 0 && choiceB1 != 1) {
			printf("\ntry again: ");
			scanf("%d", &choiceB1);
		}
		choiceB = choiceB1;
	}
	if (choiceB == 1) {
		printf("Please initialize the vector \"b\".\n");
		int n, temp;
		for(n = 0; n < N; ++n) {
			printf("b[%d] = ", n);
			scanf("%d", &temp);
			b[n] = temp;
		}
	}
	else {
		for(n = 0; n < N; ++n) {
			b[n] = 0;
		}
	}
}

// computes y to solve for x
void calcy(int A[][M], int vec[], int y[]) {
	int m, n = 0;
	consistent = 1;
	for(m = 0; m < M; ++m) {
		if (m >= rank) {
			y[m] = 0;
		}
		else {
			if (vec[m] % A[m][m] == 0) {
				y[m] = vec[m] / A[m][m];
			}
			else {
				printf("This system is not consistent\n");
				consistent = 0;
				break;
			}
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
		// this swaps entries without creating a temp
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
		// this swaps entries without creating a temp
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
		// this swaps entries without creating a temp
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
			// this swaps entries without creating a temp
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
			// this swaps entries without creating a temp
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

// compares objects for sorting
int comp(const void *pa, const void *pb) {   
	int a = *(const int*)pa;
   	int b = *(const int*)pb;
	if (a==b) {
    	return 0;
	}
  	else if (a < b) {
    	return -1;
    }
    else {
    	return 1;
	}
}


void printArray(int A[], int len) {
	int i;
	if (len == 1) {
		printf("(%d)\n", A[0]);
	}
	else {
		for(i = 0; i < len; ++i) {
			if(i == 0) { printf(" (%i, ", A[i]); }
			else if(i == len - 1) { printf("%i)\n", A[i]); }
			else { printf("%i, ", A[i]); }
		}
	}
}

void printEDs(int A[], int len) {
	int i, s;
	// s is assigned to the number of elementary divisors
	for(i = 0; i < len; ++i) {
		if (A[i] != 0) {
			s = i;
			break;
		}
	}
	printf("The set of elementary divisors is: \n");
	for(i = s; i < len; ++i) {
		if(i == s) { printf(" (%i, ", A[i]); }
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

// returns 1 if the entry divides its entire row and col; 0 otherwise
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
	if (counter >= rank) {
		return;
	}
	else {
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
			printf("This should not execute Error in findLeastEntry2.\n");
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
}

// combines all type 1 row operations into one function
void rowOperations1(int A[][M], int P[][N], int Pinv[][N], int n, int unit) {
	type1rowM(A, n, unit);
	type1rowN(P, n, unit);
	type1rowN(Pinv, n, unit);
}

// combines all type 2 row operations into one function
void rowOperations2(int A[][M], int P[][N], int Pinv[][N], int n, int tempN, int tempMult) {
	type2rowM(A, n, tempN, -1*tempMult);
	type2rowN(P, n, tempN, -1*tempMult);
	type2rowN(Pinv, tempN, n, tempMult);
}

// combines all type 3 row operations into one function
void rowOperations3(int A[][M], int P[][N], int Pinv[][N], int row1, int row2) {
	if (row1 != row2) {
		// row2 < min(N, M) && row1 < min(N, M) && 
		type3rowM(A, row1, row2);
		type3rowN(P, row1, row2);
		type3rowN(Pinv, row1, row2);
	}
}

// combines all type 2 column operations into one function
void columnOperations2(int A[][M], int Q[][M], int Qinv[][M], int m, int tempM, int tempMult) {
	type2col(A, N, m, tempM, -1*tempMult);
	type2col(Q, M, m, tempM, -1*tempMult);
	type2col(Qinv, M, tempM, m, tempMult);
}

// combines all type 2 column operations into one function
void columnOperations3(int A[][M], int Q[][M], int Qinv[][M], int col1, int col2) {
	if (col1 != col2) {
		//col2 < min(N, M) && col1 < min(N, M) &&
		type3col(A, N, col1, col2);
		type3col(Q, M, col1, col2);
		type3col(Qinv, M, col1, col2);
	}
}

void smithTransform(int A[][M], int P[][N], int Pinv[][N], int Q[][M], int Qinv[][M], int pos) {
	int tempMult;
	columnOperations2(A, Q, Qinv, pos, pos + 1, -1);

	while (A[pos+1][pos] != 0) {
		tempMult = eucDiv(A[pos+1][pos], A[pos][pos]);
		
		// runs Euclidean division
		rowOperations2(A, P, Pinv, pos+1, pos, tempMult);
		
		// makes A[pos][pos] entry larger than A[pos][pos+1]
		if (A[pos+1][pos] != 0) { 
			rowOperations3(A, P, Pinv, pos, pos+1);
		}
	}
	tempMult = eucDiv(A[pos][pos+1], A[pos][pos]);
	columnOperations2(A, Q, Qinv, pos+1, pos, tempMult);
	if (A[pos][pos+1] != 0) {
		printf("smithTransform failed. This should never happen");
	}
}

void makeAllDiagsPositive(int A[][M], int P[][N], int Pinv[][N]) {
	int n, m;
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if(A[n][m] < 0) {
				rowOperations1(A, P, Pinv, n, -1);
			}
		}
	}
}

void orderDiagonals(int A[][M], int P[][N], int Pinv[][N], int Q[][M], int Qinv[][M]) {
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
int checkSmith(int A[][M]) {
	int smith = -1;
	int n;
	for(n = 0; n < rank - 1; ++n) {
		if (A[n+1][n+1] % A[n][n] != 0) {
			smith = n;
			break;
		}
	}
	return smith;
}

// calculates the rank of a diagonal matrix
int getRank(int A[][M]) {
	int rank = 0;
	int n, m;
	for (n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			if (A[n][m] != 0) {
				++rank;
			}
		}
	}
	return rank;
}

void matxvecN(int P[][N], int b[], int vec[]) {
	int i, j, temp;
	for(i = 0; i < N; ++i) {
		temp = 0;
		for(j = 0; j < N; ++j) {
			temp += P[i][j]*b[j];
		}
		vec[i] = temp;

	}
}

void matxvecM(int Q[][M], int y[], int vec[]) {
	int i, j, temp;
	for(i = 0; i < M; ++i) {
		temp = 0;
		for(j = 0; j < M; ++j) {
			temp += Q[i][j]*y[j];
		}
		vec[i] = temp;
	}
}

void matxvecNM(int A[][M], int x[], int vec[]) {
	int n, m, temp;
	for(n = 0; n < N; ++n) {
		temp = 0;
		for(m = 0; m < M; ++m) {
			temp += A[n][m]*x[m];
		}
		vec[n] = temp;
	}
}

void matNNxmatNM(int P[][N], int A[][M], int PA[][M]) {
	int n, m, n1, temp;
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			temp = 0;
			for(n1 = 0; n1 < N; ++n1) {
				temp += P[n][n1]*A[n1][m];
			}
			PA[n][m] = temp;
		}
	}
}

void matNMxmatMM(int A[][M], int Q[][M], int AQ[][M]) {
	int n, m, m1, temp;
	for(n = 0; n < N; ++n) {
		for(m = 0; m < M; ++m) {
			temp = 0;
			for(m1 = 0; m1 < M; ++m1) {
				temp += A[n][m1]*Q[m1][m];
			}
			AQ[n][m] = temp;
		}
	}
}

void matNNxmatNN(int A[][N], int B[][N], int AB[][N]) {
	int n1, n2, n, temp;
	for(n1 = 0; n1 < N; ++n1) {
		for(n2 = 0; n2 < N; ++n2) {
			temp = 0;
			for(n = 0; n < N; ++n) {
				temp += A[n1][n]*B[n][n2];
			}
			AB[n1][n2] = temp;
		}
	}
}


void matMMxmatMM(int A[][M], int B[][M], int AB[][M]) {
	int m1, m2, m, temp;
	for(m1 = 0; m1 < M; ++m1) {
		for(m2 = 0; m2 < M; ++m2) {
			temp = 0;
			for(m = 0; m < M; ++m) {
				temp += A[m1][m]*B[m][m2];
			}
			AB[m1][m2] = temp;
		}
	}
}

// slight customization for Prof. Zhou's prime factorization algorithm
void primeFactorize(int Num, int prime[]) {
	int i, j, k;
	int index;
	int start;
	int primePower[maxED];
	initializeZero(primePower, maxED);
	int tempEntry;
	for(j = 0; j < maxED; ++j) {
		if(prime[j] == 0) {
			index = j;
			start = j;
			break;
		}
	}

	//check if 2 is a factor
	if(Num > 1 && Num % 2 == 0) {
		prime[index] = 2;
	}

	//take out all factors of 2
	while (Num > 1 && Num % 2 == 0) {
		primePower[index] += 1;
		Num = Num / 2;
	}

	while (Num > 1 && i > 0 && 4*i*i < Num) {
		if (Num % (2*i + 1) == 0) {
			index++;
			prime[index] = 2*i + 1;
		}
		while (Num > 1 && Num % (2*i + 1) == 0) {
			primePower[index]++;
			Num = Num / (2*i + 1);
		}
		++i;
	}

	// it should actually never hit here...
	if(Num > 1) {
		index++;
		primePower[index] = 1;
		prime[index] = Num;
	}

	for(j = start; j < index + 1; ++j) {
		tempEntry = prime[j];
		for(k = 1; k < primePower[j]; ++k) {
			prime[j] *= tempEntry;
		}
	}
}