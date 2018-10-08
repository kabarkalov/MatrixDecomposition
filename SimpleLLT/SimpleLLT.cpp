// SimpleLLT.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <time.h>

void Cholesky_Decomposition(double * A, double * L, int N)
{
	memset(L, 0, sizeof(double)*N*N);
	for (int i = 0; i < N; i++)
	{
		L[i*N + i] = A[i*N + i];
		for (int k = 0; k < i; k++)
		{
			L[i*N + i] -= L[i*N + k] * L[i*N + k];
		}
		L[i*N + i] = sqrt(L[i*N + i]);
		
		for (int j = i + 1; j < N; j++)
		{			
			L[j*N + i] = A[j*N + i];
			for (int k = 0; k < i; k++)
			{
				L[j*N + i] -= L[i*N + k] * L[j*N + k];
			}
			L[j*N + i] /= L[i*N + i];
		}
	}
}

void PrintResult(double * A, double * L, int N)
{
	//����� �����������
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			printf("%3.0lf ",A[i * N + j]);
		}
		printf("\n");
  }
	printf("\n");

	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{

      if(j>i)
        printf("   0  ");
      else 
        printf("%2.3lf ",L[i * N + j]);
		}
    printf("\n");
  }
  printf("\n");
	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{
      if(i>j)
        printf("   0  ");
      else
			  printf("%2.3lf ",L[j * N + i]);
		}
		printf("\n");
	}
}

int main(int argc, char* argv[])
{
  const int N = 10;
	double A[N*N], L[N*N];
  time_t begin, end;
	//��������� �������� ������������ ������� A > 0
	for(int i=0;i<N;i++)
	{
  	for(int j=i+1;j<N;j++)
    {
		  A[i*N+j] = rand()%100 + 1;
      A[j*N+i] = A[i*N+j];
    }
	}
	for(int i=0;i<N;i++)
	{
    double sum = 0;
    A[i*N+i]=0;
  	for(int j=0;j<N;j++)
    {
		  sum += A[i*N+j];
    }
    A[i*N+i]=sum;
	}

  begin = clock();
	//����������
	Cholesky_Decomposition(A,L,N);
  end = clock();

  PrintResult(A,L,N);

  printf(" Total time = %f sec.\n", (end - begin)/1000.0);

	return 0;
}

