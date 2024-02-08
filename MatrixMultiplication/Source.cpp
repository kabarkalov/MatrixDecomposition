#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

void matmul(double* A, double* B, double* C, int N)
{
	int i, j, k;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) 
			for (k = 0; k < N; k++) 
				C[i*N + j] += A[i*N + k] * B[k*N + j];
}

void matmulblock(double* A, double* B, double* C, int Size, int BlockSize = 50) {
	int GridSize = int(Size / double(BlockSize));
	for (int n = 0; n < GridSize; n++)
		for (int m = 0; m < GridSize; m++)
			for (int s = 0; s < GridSize; s++)
				for (int i = n * BlockSize; i < (n + 1) * BlockSize ; i++)
					for (int j = m * BlockSize; j < (m + 1) * BlockSize ; j++)
						for (int k = s * BlockSize; k < (s + 1) * BlockSize ; k++)
							C[i * Size + j] += A[i * Size + k] * B[k * Size + j];
}

void matmulblockp(double* A, double* B, double* C, int Size, int BlockSize = 50) {
	int GridSize = int(Size / double(BlockSize));
	#pragma omp parallel for
	for (int n = 0; n < GridSize; n++)
		for (int m = 0; m < GridSize; m++)
			for (int s = 0; s < GridSize; s++)
				for (int i = n * BlockSize; i < (n + 1) * BlockSize; i++)
					for (int j = m * BlockSize; j < (m + 1) * BlockSize; j++)
						for (int k = s * BlockSize; k < (s + 1) * BlockSize; k++)
							C[i * Size + j] += A[i * Size + k] * B[k * Size + j];
}


void PrintResult(double* A, double* L, double* U, int N)
{
	//Вывод результатов
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(6);
			cout << A[i * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(3);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(6);
			cout << L[i * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(6);
			cout << U[i * N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


int main(int argc, char* argv[])
{
	const int N = 2000;
	double* A = new double[N * N];
	double* B = new double[N * N];
	double* C = new double[N * N];
	time_t begin, end;

	//Генерация тестовых матриц
	for (int i = 0; i < N * N; i++)
	{
		A[i] = rand() % 10 + 1;
		B[i] = rand() % 10 + 1;
		C[i] = 0;
	}

	begin = clock();
	//Произведение
	matmul(A, B, C, N);
	//matmulblock(A, B, C, N);
	end = clock();
//	PrintResult(A, B, C, N);
	cout << " Total time = " << (end - begin) / 1000.0 << " sec." << endl;

	delete[] A;
	delete[] B;
	delete[] C;

	return 0;
}
