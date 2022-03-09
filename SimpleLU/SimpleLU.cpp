// SimpleLU.cpp : Defines the entry point for the console application.
//

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

void LU_DecompositionP(double* A, double* L, double* U, int *p, int N)
{
	for (int i = 0; i < N; i++)
	{
		p[i] = i;
		for (int j = 0; j < N; j++)
		{
			U[N * i + j] = A[N * i + j];
		}
	}
	for (int i = 0; i < N; i++)
	{
		//поиск позиции максимального элемента в столбце
		int pos = i;
		for (int k = i + 1; k < N; k++)
		{
			if (fabs(U[p[k] * N + i]) > fabs(U[p[pos] * N + i]))
				pos = k;
		}
		//перестановка
		int tmp = p[i];
		p[i] = p[pos];
		p[pos] = tmp;

		L[p[i] * N + i] = 1;
		for (int k = i + 1; k < N; k++)
		{
			double mu = U[p[k] * N + i] / U[p[i] *N + i];
			for (int j = i; j < N; j++)
			{
				U[p[k] * N + j] -= mu * U[p[i] * N + j];
			}
			L[p[k] * N + i] = mu;
			L[p[i] * N + k] = 0;
		}
	}
	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			U[p[i] * N + j] = 0;
		}
	}
}



void PrintResult(double * A, double * L, double * U, int N)
{
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
	cout << endl;
}

void PrintResultP(double* A, double* L, double* U,int *p, int N)
{
	//Вывод результатов
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(4);
			cout << A[i * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl << "Permutation = ";
	for (int i = 0; i < N; i++)
	{
		cout << p[i] << " ";
	}
	cout << endl<<endl;
	
	cout.precision(3);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(6);
			cout << L[p[i] * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(6);
			cout << U[p[i] * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


int main(int argc, char* argv[])
{
	const int N = 20;
	double* A = new double[N*N]; 
	double* L = new double[N*N]; 
	double* U = new double[N*N];

	int* p = new int[N];

	time_t begin, end;

	//Генерация тестовой матрицы A
	for(int i=0;i<N*N;i++)
	{
		A[i] = rand()%100 + 1;
	}
	begin = clock();
	//Разложение
	//LU_Decomposition(A, L, U, N);
	LU_DecompositionP(A,L,U,p,N);
	end = clock();

//	PrintResult(A, L, U, N);

	PrintResultP(A,L,U,p,N);
    
	cout<<" Total time = "<< (end - begin)/1000.0<<" sec."<<endl;

	delete[] A;
	delete[] L;
	delete[] U;
	delete[] p;

	return 0;
}

