// SimpleQR.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

 //Умножение матрицы A на матрицу Хаусхолдера (правый нижний квадрат начиная с номера k)
 void RowHouse(double *A, double *v, int N, int k)
 {
   double dp = 0;
   //Скалярное произведение (v,v)
   for(int i=k;i<N;i++)
   {
     dp += v[i]*v[i];
   }
   //Коэффициент бета
   double beta = -2/dp;
   //Матрично-векторное произведение A^t*v
   double *w = new double[N];
   for(int j=k;j<N;j++)
   {
     w[j] = 0;
     for(int i=k;i<N;i++)
     {
       w[j] += A[i*N + j]*v[i];
     }
     w[j] *= beta;
   }
   //Обновление матрицы
   for(int i=k;i<N;i++)
   {
     for(int j=k;j<N;j++)
     {
       A[i*N + j] += v[i]*w[j]; 
     }
   }
   delete[] w;
 }

//QR-разложение
void QR_Decomposition(double * A, double * Q, double * R, int N)
{
  double * t=new double[N];
  int j,i,k;
  double norm,d,sum,alpha,beta,tau;

  for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			R[N * i + j] = A[N * i + j];
		}
	}

  time_t begin, end;
  begin = clock();
  for (j = 0; j < N; j++)
  {
    /* Вычисление преобразований Хаусхолдера, чтобы привести j-ый столбец
       матрицы к множителю j-ого единичного вектора. */

    norm = 0.;
    for(i = j+1; i < N; i++)
    {
      d = fabs(R[i*N + j]);
      if( d > norm)
      {
        norm = d;
      }
    }
    sum = 0.;
    for(i = j+1; i < N; i++)
    {
      d = R[i*N + j]/norm;
      sum += d*d;
    }
    norm = sqrt(sum)*norm;


    alpha = R[j*N+j];
    beta = -(alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, norm) ;
    if(beta == 0.)
    {
      cout<<" matrix is singular"<<endl;
      return;
    }
    t[j] = tau = (beta-alpha)/beta;
    d   = (alpha-beta);
    for(i = j+1; i < N; i++)
    {
      R[i*N + j] /= d;
    }
    R[j*N + j] = beta;
    for(k = j+1; k < N; k++)
    {
      sum = R[j*N + k];
      for(i = j+1; i < N; i++)
        sum += R[i*N + k]*R[i*N + j];
      R[j*N + k] -= tau*sum;
      for(i = j+1; i < N; i++)
        R[i*N+k] -= tau*sum*R[i*N+j];
    }
   }
  end = clock();

  cout<<"QR time = "<< (end - begin)/1000.0<<"sec."<<endl;

//  double * Q_k=new double[N*N];
//  double * M=new double[N*N];
  double * v=new double[N];

  //Быстрое формирование Q обратным произведением
  for(i = 0; i < N*N; i++)
  {
      Q[i] = 0.0;
  }
  for(i = 0; i < N; i++)
  {
      Q[i*N + i] = 1;
  }

  for(j = N-1; j >=0; j--)
  {
    v[j] = 1;
    for(i=j+1;i<N;i++)
    {
      v[i] = R[i*N + j];
    }
    RowHouse(Q,v,N,j);
    //Печать текущей матрицы Q
   // cout<<"Q"<<k<<endl;
	  //for(i=0;i<N;i++)
   // {
		 // for(k=0;k<N;k++)
		 // {
			//  cout.width(6);
			//  cout<<Q[i * N + k]<<" ";
		 // }
		 // cout<<endl;
	  //}
   // cout<<endl;

  }

/*

  //Формирование Q обратным произведением
  for(k = N-1; k >=0; k--)
  {
    double tk = t[k];    
    
    Q_k[k*N + k] = 1-tk;
    for(i = k+1; i < N; i++)
    {
      Q_k[k*N + i] = Q_k[i*N + k] = -R[i*N + k]*t[k];
    }
    for(j = k+1; j < N; j++)
    {
      for(i = k+1; i < j; i++)
        Q_k[i*N+j] = Q_k[j*N + i] = -R[i*N+k]*R[j*N + k]*t[k];
      Q_k[j*N + j] = 1 - R[i*N + k]*R[j*N + k]*t[k];
    }

    //Печать текущей матрицы Q
    cout<<"Q"<<k<<endl;
	  for(i=0;i<N;i++)
    {
		  for(j=0;j<N;j++)
		  {
			  cout.width(6);
			  cout<<Q[i * N + j]<<" ";
		  }
		  cout<<endl;
	  }
    cout<<endl;
    

    if(k==N-1)
    {
      for(i = 0; i < N*N; i++)
      {
          Q[i] = 0.0;
          M[i] = 0.0;
      }
      for(i = 0; i < N-1; i++)
      {
          Q[i*N + i] = 1;
          M[i*N + i] = 1;
      }
      Q[(N-1)*N + N-1] =  Q_k[(N-1)*N + N-1];      
    }else
    { //M = Q_k*Q
      for(i=k;i<N;i++)
      {
        for(j=k;j<N;j++)
        {
          M[i*N+j] = 0;
          for(int s=k;s<N;s++)
          {
            M[i*N+j] +=  Q_k[i*N + s] * Q[s*N + j];
          }          
        }
      }
      //Q = M
      for(i=k;i<N;i++)
      {
        for(j=k;j<N;j++)
        {          
          Q[i*N+j] =  M[i*N + j];
        }
      }
    }
 
  }
 */ 
  /*
  //формирование Q прямым произведением
  for(k = 0; k < N; k++)
  {
    double tk = t[k];
    for(j = 0; j < k; j++)
    {
      for(i = 0; i < k; i++)
      {
        Q_k[i*N+j] = Q_k[j*N+i] = 0.;
      }
      Q_k[j*N+j] = 1.;
      for(i = k; i < N; i++)
      {
        Q_k[i*N + j] = Q_k[j*N + i] = 0.;
      }
    }
    Q_k[k*N + k] = 1-tk;
    for(i = k+1; i < N; i++)
    {
      Q_k[k*N + i] = Q_k[i*N + k] = -R[i*N + k]*t[k];
    }
    for(j = k+1; j < N; j++)
    {
      for(i = k+1; i < j; i++)
        Q_k[i*N+j] = Q_k[j*N + i] = -R[i*N+k]*R[j*N + k]*t[k];
      Q_k[j*N + j] = 1 - R[i*N + k]*R[j*N + k]*t[k];
    }
    
   // //Печать Qk
   // cout<<"Q"<<k<<endl;
	  //for(i=0;i<N;i++)
   // {
		 // for(j=0;j<N;j++)
		 // {
			//  cout.width(6);
			//  cout<<Q_k[i * N + j]<<" ";
		 // }
		 // cout<<endl;
	  //}
   // cout<<endl;
   // 

    if(k==0)
    {
      for(i=0;i<N*N;i++)
      {
        Q[i] = Q_k[i];
      }
    }else
    { //M = Q*Q_k
      for(i=0;i<N;i++)
      {
        for(j=0;j<N;j++)
        {
          M[i*N+j] = 0;
          for(int s=0;s<N;s++)
          {
            M[i*N+j] += Q[i*N + s] * Q_k[s*N + j];
          }
        }
      }
      //Q = M
      for(i=0;i<N*N;i++)
      {
        Q[i] = M[i];
      }
    }
  }
  */


//  delete[] Q_k;
//  delete[] M;
  delete[] v;
  
  //Обнуление нижнего треугольника R
  for(i=1;i<N;i++)
  {
    for(j=0;j<i;j++)
    {
      R[i*N + j] = 0;
    }
  }
 }

 void PrintResult(double *A, double *Q, double *R, int N)
 {
	double* T = new double[N*N]; 
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
			cout.width(7);
			cout<<Q[i * N + j]<<" ";
		}
    cout<<endl;
  }
  cout<<endl;
  
	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{
			cout.width(7);
			cout<<R[i * N + j]<<" ";
		}
		cout<<endl;
	}
  cout<<endl;

  //проверка ортогональности
  //T = QT*Q
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      T[i*N+j] = 0;
      for(int k=0;k<N;k++)
      {
        T[i*N+j] += Q[i*N + k] * Q[j*N + k];
      }
    }
  }

  cout<<endl;
 	cout.precision(3); 
	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{
			cout.width(7);
			cout<<fixed<<T[i * N + j]<<" ";
		}
		cout<<endl;
	}

  //проверка корректности
  //A = Q*R
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      A[i*N+j] = 0;
      for(int k=0;k<N;k++)
      {
        A[i*N+j] += Q[i*N + k] * R[k*N + j];
      }
    }
  }

  cout<<endl;
 	cout.precision(3); 
	for(int i=0;i<N;i++)
  {
		for(int j=0;j<N;j++)
		{
			cout.width(7);
			cout<<fixed<<A[i * N + j]<<" ";
		}
		cout<<endl;
	}
  delete[] T;
 }

 

int _tmain(int argc, _TCHAR* argv[])
{
  const int N = 10;
	double* A = new double[N*N]; 
	double* Q = new double[N*N]; 
  double* R = new double[N*N]; 

  time_t begin, end;
	//Генерация тестовой матрицы A
	for(int i=0;i<N*N;i++)
	{
		A[i] = rand()%200 + 1;
	}
  //Разложение
  begin = clock();
	QR_Decomposition(A,Q,R,N);
  end = clock();

  //Печать результатов
  PrintResult(A,Q,R,N);

  cout<<" Total time = "<< (end - begin)/1000.0<<"sec."<<endl;
  delete[] A;
  delete[] Q;
  delete[] R;
	return 0;
}

