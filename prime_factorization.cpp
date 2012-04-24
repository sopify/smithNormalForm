// module theory.cpp primary factorization
//

// #include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
// #include <conio.h>
#include <math.h>
#include <string.h>
#include <time.h>
using namespace std;




int main()
{
	while(1)
	{
		cout << "Input an integer to decompose:" << endl;

		unsigned int PCount = 0;
		unsigned long int prime[50] = {0};
		unsigned int primepower[50] = {0};

		unsigned long long int Num;// input 502400; the number to be decomposed =2^7 * 5^2 * 157^1
		cin >> Num;

		unsigned long int i = 1;

		if(Num > 1 && Num % 2 == 0)//determine if 2 is a factor
			prime[PCount] = 2;


		while(Num >1 && Num % 2 == 0)
			{
				primepower[0]++;
				Num = Num/2;
			}

		
		while(Num > 1 && i>0 && 4*i*i<Num)
		{		
			if(Num % (2*i + 1) == 0)
			{
				PCount++; //a new prime factor found
				prime[PCount] = 2*i+1;

			}

			while(Num > 1 && Num % (2*i + 1) == 0)
			{
				primepower[PCount]++;
				Num = Num / (2*i + 1);
			}
			i++;

		}

		if(Num > 1)
		{
			PCount++;
			primepower[PCount] = 1;
			prime[PCount] = Num;
		}
	

			//print results
		if(primepower[0])
		printf("%d^%d\n",2,primepower[0]);
		for(i=1; i  <= PCount; i++)
			if(primepower[i])
				printf("%lu^%u\n",prime[i], primepower[i]);
	}
	return 0;
}