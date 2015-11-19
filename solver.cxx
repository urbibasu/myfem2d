#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <omp.h>
using namespace std;
#include "solver.hpp"

void test_system
( int n, double **A, double *b )
{
    int i,j;
    //fprintf(stderr,"In %s: (0)\n",__func__);
    for(i=0; i<n; i++)
    {
        //fprintf(stderr,"In %s: (1) i=%d/%d\n",__func__,i,n);
        b[i] = 2.0*n;
#pragma omp parallel
        for(j=0; j<n; j++)
        {
            //fprintf(stderr,"In %s: (2) j=%d/%d\n",__func__,j,n);
            if (j==i)
                A[i][j] = 2.0;
            else if (j==i-1)
                A[i][j]= -1;
            else if (j==i+1)
                A[i][j]= -1;
            else
                A[i][j]= 0;
        }
        
        A[i][i] = n + 1.0;
    }
}

//------------------------------------------------
void jacobi(int n, double **A, double *b, double epsilon, int maxit, double *x )
{
    double y[n+1];
    int i, j, k;
    double diag, sum;
    printf("jacobi\n");
    for (k = 1; k <= maxit; k++)
    {
        cp_vector(n, y, x);
        for (i = 0; i < n; i++)
        {
            sum = b[i];
            diag = A[i][i];
            if (fabs(diag) < epsilon)
            {
                printf("diagonal element too small\n");
                return;
            }
            for (j = 0; j < n; j++)
                if (j != i)
                    sum -= A[i][j]*y[j];
            x[i] = sum/diag;
            //printf("%d:\tx[%d] = %2.4f\n", k, i, x[i]);
        }
        if (epsilon_comp(n, sub_vector(n, x, y)))
        {
            printf("k = %d\n", k);
            for (i = 0; i < n; i++)
                printf("x[%d] = %2.4f\n", i, x[i]);
            return;
        }
    }
    printf("maximum iterations reached\n");
    return;
}

//--------------------------------------------
void cp_vector(int size, double* y, const double* x)
{
    int i;
    
    for (i = 1; i <= size; i++)
        y[i] = x[i];
}

//---------------------------------------------
double* sub_vector(int size, const double* x, const double* y)
{
    int i;
    double* d;
    
    d = new double[size+1];
    for (i = 1; i <= size; i++)
        d[i] = x[i] - y[i];
    return d;
}
//--------------------------------------------
bool epsilon_comp(int size, const double* x)
{
    const double epsilon = 0.00005;
    int i, count;
    
    count = 0;
    for (i = 1; i <= size; i++)
        if (x[i] < epsilon)
            count++;
    if (count == size)
        return true;
    else
        return false;
}

//------------gauss
void run_gauss_seidel_method
( int n, double **A, double *b,
 double epsilon, int maxit,
 int *numit, double *x )
{
    printf("gauss\n");
    double *dx = (double*) malloc(n*sizeof(double));
    int i,j,k;
    
    //fprintf(stderr,"maxit=%d n=%d\n",maxit,n);
    for(k=0; k<maxit; k++)
    {
        double sum = 0.0;
        for(i=0; i<n; i++)
        {
            dx[i] = b[i];
            for(j=0; j<n; j++)
                dx[i] -= A[i][j]*x[j];
            dx[i] /= A[i][i]; x[i] += dx[i];
            sum += ( (dx[i] >= 0.0) ? dx[i] : -dx[i]);
            
        }
        if(sum <= epsilon) break;
        printf("x[%d] = %2.4f\n", k, sum);
    }
    *numit = k+1; free(dx);
}

//---------------------------------------COnJUGATE GRADIEnT
void cg(int n, double **A, double *b, double epsilon, int maxit,double *x)
{
    double dot[n];
    double r[n],p0[n];
    double A_p0[n],beta,rk1;
    // zero the dot array
    for(unsigned int i=0; i < n; i++)
    {
        dot[i]=0.0;
    }
    
    // compute A * x0
    for(unsigned int i = 0; i < n; i++)
    {
        for(unsigned int j = 0; j < n; j++)
        {
            dot[i] += A[i][j] * x[j];
        }
        //		cout<<dot[i]<<" ";
        //		std::cout << "\n";
    }
    // compute residual
    for(unsigned int i = 0; i < n; i++)
    {
        r[i] = b[i] - dot[i];
        //		std::cerr << r[i]<<" ";
    }
    //	cout<<"\n";
    for(int i=0; i<n; i++)
    {
        p0[i]=r[i];
        //    	cout<<p0[i]<<" ";
    }
    for(unsigned int step=0; step<1000000; step++)
    {
        
        //		double r_tran[1][2]={0.0,0.0};
        //		double r_tran_A[1][2]={0.0,0.0};
        double upper1=0.0;
        double upper2=0.0;
        double lower1=0.0;
        double lower2=0.0;
        double sum=0.0;
        for(unsigned int i = 0; i < n; i++)
        {
            upper1 += p0[i] * r[i];
            //			cout<<upper1<<"\n";
        }
        //
        
        // r^T*A*r
        for(unsigned int i = 0; i < n; i++)
        {
            for(unsigned int j = 0; j < n; j++)
            {
                lower1 +=  p0[i] * A[i][j] * p0[j];
            }
        }
        //		cout<<lower1<<"\n";
        // new residual
        if(lower1 != 0)
        {
            rk1=upper1/lower1;
            //    	  cout<<rk1<<"\n";
        }
        else
        {
            std::cerr << "R^T*A*r is ZERO!! Aborting!!" << "\n";
            //        	exit(-1);
        }
        
        // update the solution
        for(unsigned int i = 0; i < n; i++)
        {
            x[i] += rk1 * p0[i];
            //    	  cout<<x[i]<<" ";
        }
        //        cout<<"\n";
        for(unsigned int i=0; i < n; i++)
        {
            A_p0[i]=0.0;
        }
        //update the residual
        for(unsigned int i = 0; i < n; i++)
        {
            for(unsigned int j = 0; j < n; j++)
            {
                A_p0[i] += A[i][j] * p0[j];
            }
        }
        for(int i=0;i<n;i++)
        {
            r[i]=r[i]-rk1*A_p0[i];
            //        	cout<<r[i]<<" ";
        }
        
        for(unsigned int i = 0; i < n; i++)
        {
            upper2 += p0[i] * r[i];
        }
        //		cout<<upper2<<"\n";
        for(unsigned int i = 0; i < n; i++)
        {
            for(unsigned int j = 0; j < n; j++)
            {
                lower2 +=  p0[i] * A[i][j] * p0[j];
            }
        }
        //		cout << lower2;
        if(lower2 != 0)
        {
            beta=upper2/lower2;
        }
        //update p0;
        for(int i=0;i<n;i++)
        {
            p0[i]=r[i]-beta*p0[i];
        }
        // computing the norm of the residual.
        for(int i=0;i<n;i++)
        {
            sum = sum+pow(r[i],2);
        }
        cout<<sqrt(sum)<<"\n";
        if( sqrt(sum) <= epsilon)
        {
            std::cout <<"n step to converge(CG method):  " << step <<"\n" ;
            std::cout <<"the epsilon is : " << epsilon << "\n" ;
            break;
        }
    }
}





