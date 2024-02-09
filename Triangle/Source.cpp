#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

void solve(double* A, double* x, double* b, int N)
{
	int i, j, k;
	for (i = N - 1; i >= 0; i--)
	{
		x[i] = 0;
		for (j = i + 1; j < N; j++)
		{
			x[i] -= A[i * N + j] * x[j];
		}
		x[i] += b[i];
		x[i] /= A[i * N + i];
	}
}

void PrintResult(double* A, double* x, double* b, int N)
{
	cout.precision(3);
	//Вывод результатов
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			cout.width(4);
			cout << 0 << " ";
		}
		for (int j = i; j < N; j++)
		{
			cout.width(4);
			cout << A[i * N + j] << " ";
		}
		cout.width(4);
		cout << x[i] << "   " << b[i];
		cout << endl;
	}
	cout << endl;
}

double CheckResult(double* x, int N)
{
	double err = 0;
	for (int i = 0; i < N; i++)
	{
		if (fabs(x[i] - 0.1) > err)
		{
			err = fabs(x[i] - 0.1);
		}
	}
	return err;
}


int main(int argc, char* argv[])
{
	const int N = 100;
	double* A = new double[N * N];
	double* x = new double[N];
	double* b = new double[N];
	time_t begin, end;

	//Генерация тестовой матрицы
	for (int i = 0; i < N; i++)
	{
		b[i] = 0;
		for (int j = i; j < N; j++)
		{
			A[i*N+j] = rand() % 10 + 1;
			b[i] += A[i * N + j];
		}
		b[i] *= 0.1;
	}

	begin = clock();
	//Произведение
	solve(A, x, b, N);
	end = clock();
//	PrintResult(A, x, b, N);
	cout << "Error = " << scientific <<CheckResult(x, N)<<endl;
	cout << "Total time = " << (end - begin) / 1000.0 << " sec." << endl;

	delete[] A;
	delete[] x;
	delete[] b;

	return 0;
}
