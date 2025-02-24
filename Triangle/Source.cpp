#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

template <typename T>
void solve(T* A, T* x, T* b, int N)
{
	//Решение треугольной системы
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
template <typename T>
void PrintResult(T* A, T* x, T* b, int N)
{
	cout << " A x = b"<<endl;
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
		cout <<" | " << x[i] << " | " << b[i];
		cout << endl;
	}
	cout << endl;
}

template <typename T>
T CheckResult(T* x, int N)
{
	T err = 0;
	for (int i = 0; i < N; i++)
	{
		if (fabs(x[i] - (T)0.1) > err)
		{
			err = fabs(x[i] - (T)0.1);
		}
	}
	return err;
}


int main(int argc, char* argv[])
{
	const int N = 50;
	double* A = new double[N * N];
	double* x = new double[N];
	double* b = new double[N];
	time_t begin, end;

	//Генерация тестовой задачи
	for (int i = 0; i < N; i++)
	{
		b[i] = 0;
		for (int j = i; j < N; j++)
		{
			A[i*N+j] = rand() % 100 + 1.0;
			b[i] += A[i * N + j];
		}
		b[i] *= 0.1;
	}

	begin = clock();
	//Решение системы
	solve(A, x, b, N);
	end = clock();
//	PrintResult(A, x, b, N);
	cout << "Error = " << scientific << CheckResult(x, N) <<endl;
	cout << "Total time = " << (end - begin) / 1000.0 << " sec." << endl;

	delete[] A;
	delete[] x;
	delete[] b;

	return 0;
}
