// SimpleLU.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

void LU_Decomposition(double * A, double * L, double * U, int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			U[N * i + j] = A[N * i + j];
		}
	}
	for (int i = 0; i < N; i++)
	{
		L[i * N + i] = 1;
		for (int k = i + 1; k < N; k++)
		{
			double mu = U[k * N + i] / U[N * i + i];
			for (int j = i; j < N; j++)
			{
				U[k * N + j] -= mu * U[i * N + j];
			}
			L[k * N + i] = mu;
			L[i * N + k] = 0;
		}
	}
	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			U[i * N + j] = 0;
		}
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	const int N = 10;
	double* A = new double[N*N]; 
	double* L = new double[N*N]; 
  double* U = new double[N*N]; 
  time_t begin, end;

	//Генерация тестовой матрицы A
	for(int i=0;i<N*N;i++)
	{
		A[i] = rand()%100 + 1;
	}
  begin = clock();
	//Разложение
	LU_Decomposition(A,L,U,N);
  end = clock();
  
	//Вывод результатов
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			cout.width(4);
			cout<<A[i * N + j]<<" ";
		}
    cout<<endl;
  }
  cout<<endl;
	cout.precision(3);
	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{
			cout.width(6);
			cout<<L[i * N + j]<<" ";
		}
    cout<<endl;
  }
  cout<<endl;
	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{
			cout.width(6);
			cout<<U[i * N + j]<<" ";
		}
		cout<<endl;
	}

  
  cout<<" Total time = "<< (end - begin)/1000.0<<"sec."<<endl;


	return 0;
}

