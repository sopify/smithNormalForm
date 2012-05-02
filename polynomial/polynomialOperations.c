#include <stdio.h>
#include "matrixOperations.h"
#include "polynomialOperations.h"

void clearZeroes(float x[]) {
	int i;
	for(i = 0; i <= maxDegree; ++i) {
		if (x[i] < precision && x[i] > -precision) {
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
	clearZeroes(c);
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
					printf("%1.2f", poly[i]);
					printf(" + ");
				}
			}
			else {
				printf("%1.2f", poly[i]);
			}
		}
		else if (i == degree(poly)) {
			printf("%1.2f*x^%d", poly[i], i);
		}
		else {
			if(poly[i] > precision || poly[i] < -precision) {
				printf("%1.2f*x^%d + ", poly[i], i);
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
	clearZeroes(tempA);
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