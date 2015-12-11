#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <omp.h>
using namespace std;
#include "solver.hpp"


void assembly (int n, int nelem)
{
    int e,i,j;
    double **A; // STIFFNESS MATRIX
    A = (double**) malloc(3*sizeof(double*));
    for(i=0; i<3; i++)
        A[i] = (double*) malloc(3*sizeof(double));
    for(i=0; i< 3; i++)
    {
        for(j=0; j< 3; j++)
        {
            A[0][0] = 2; A[0][1] = -1; A[0][2] = -1;
            A[1][0] = -1;A[1][1] = 1;  A[1][2] = 0;
            A[2][0] = -1; A[2][1] = 0; A[2][2] = 1;
        }
    }
    double **Z; // MASS MATRIX
    Z = (double**) malloc(3*sizeof(double*));
    for(i=0; i<3; i++)
        Z[i] = (double*) malloc(3*sizeof(double));
    for(i=0; i< 3; i++)
    {
        for(j=0; j< 3; j++)
        {
            Z[0][0] = 2; Z[0][1] = 1; Z[0][2] = 1;
            Z[1][0] = 1;Z[1][1] = 2;  Z[1][2] = 1;
            Z[2][0] = 1; Z[2][1] = 1; Z[2][2] = 2;
        }
    }
    double B[n][1]={};    //GLOBAL LOAD VECTOR
    double K[n][n]={};
    double M[n][n]={}; // declaring GLOBAL MASS and STIFFNESS MATRIX
    
    double **k; //element stiffness matrix
    k = (double**) malloc(nelem*sizeof(double*));
    for(e=0; e< nelem; i++)
        k[e] = (double*) malloc(nelem*sizeof(double));
    double **m; //element mass matrix
    m = (double**) malloc(nelem*sizeof(double*));
    for(e=0; e< nelem; i++)
        m[e] = (double*) malloc(nelem*sizeof(double));
    
    
    // LOOP OVER ALL ELEMENTS to create K,M, B //
    for (e=0;e<nelem;e++)
    {
        /*   SHAPE FUNCTIONS:   phi1 = 1-x-y;   phi2=x;   phi3=y; */
        ///element stiffness matirx= A *2area;
        //printf("%2.0f", k[e]);
        double loc_node[3][2] = {};
        
        /* the following loop uses the *var.coord connectivity array to get the coordinates of the nodes of an element*/
        for (i=0;i< 3;i++)
        {
            for (j=0;j<2;j++)
                const double *loc_node[i][j]= (*var.coord)[i];
        }
        
        /* the following loop loops over all nodes in an element and identifies if its on the boundary or not and changes the values in the element stiffness matrix;*/
        
        for (i=0;i<3;i++)
        {
            // looping over all nodes, i is the no of node
            if (loc_node[i][0] == 0.0 || 500.e03) || (loc_node[i][1] == 0.0 || -100.e03)
                for (j=0;j<3;j++)
                {
                    if (j==i)
                        k[i][j]==1;
                    else
                        k[i][j]==0;
                }
        }
        
        m[e]= Z * 1/12 * (*var.volume)[e]; //all element mass matrices
        k[e]= A * (*var.volume)[e];        //all element stiffness matrices
        
        
        /* creating the FORCE array; he heat source is specified in a region of {(10,-2),(11,-2),(10,-4),(11,-4)};
         THE values of the force vector at the heat source is 1 and other places is zero*/
        
        double *force;  // (3*1) array
        force = (double*) malloc(3*sizeof(double));
        for (i=0;i<3;i++)
        {
            // looping over all nodes, i is the no of node
            if (loc_node[i][0] == 10.0 || 11.0) || (loc_node[i][1] == -2.0 || -4.0)
                f[i]==1;
            else
                f[i]==0;
        }
        /* creating the load vector for each element*/
        double load_Vector[3][1] = {};
        load_Vector = (double*) malloc(3*sizeof(double));
        for (i=0;i<3;i++)
            load_Vector[i]= f[i] *1/3 * (*var.volume)[e];
        
        
        // GLOBAL mass and stiffness matrix, load vector
        const int *conn = (*var.connectivity)[e];
        K(conn(e,i),conn(e,j)) +=k;
        M(conn(e,i),conn(e,j)) +=m;
        B(conn(e,i)) +=load_Vector;
    }
    
    
    /* creating an equation of the form Ax=B to solve the linear system*/
    int lambda =5.5; //heat capacity
    
    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            new_M[i][j]= lambda * M[i][j];      //creating matrix K
        }
    }
    
    //double BB[n][1]={};  //TIMESTEP*B
    int timestep =0.0001;
    BB= timestep * B;
    double LL[n][n]={}; // new_M -timestep
    LL= new_M -timestep;
    //double L[n][n]={};    // (LL-K)
    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            L[i][j]= LL[i][j] * K [i][j];
        }
    }
}



//------------------------ creating vector new_B (matrix B of the linear system)------//


void Loadvector( int n,double **L, double *BB, double *node_temp )
{
    int i,j;
    //double new_B[n][1]={};
    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
        {
            new_B= (L[i][j] * node_temp[i]) + BB;  //final B vector of Ax=B //
            //new_B= (L[i][j] * node_temp[i]);
            
        }
       // new_B= new_B + BB;
    }


/*void test_system
( int n, double **A, double *b )
{
    int i,j;
    //fprintf(stderr,"In %s: (0)\n",__func__);
    for(i=0; i<n; i++)
    {
        //fprintf(stderr,"In %s: (1) i=%d/%d\n",__func__,i,n);
        if (i==0)
            b[i] = 0.5;
        else if (i==n-1)
            b[i] = 0.5;
        else
            b[i] = 0.0;
        //std::cout << b[i] << std::endl;
    }
#pragma omp parallel
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            //fprintf(stderr,"In %s: (2) j=%d/%d\n",__func__,j,n);
            if (j==i)
                A[i][j] = 2.0;
            else if (j==i-1)
                A[i][j]= -1.0;
            else if (j==i+1)
                A[i][j]= -1.0;
            else
                A[i][j]= 0.0;
            // std::cout << A[i][j] << std::endl;
        }
    }
}*/
//------------------------------------------------
void jacobi(int n, double **A, double *b, int maxit, int eps, double *x )
{
    double y[n];
    const double delta= 0.00000001;
    int i, j, k;
    double diag, sum;
    //printf("jacobi\n");
    
    for (k = 1; k <= maxit; k++)
    {
        cp_vector(n, y, x);
        for (i = 0; i < n; i++)
        {
            sum = b[i];
            diag = A[i][i];
            if (fabs(diag) < delta)
            {
                printf("diagonal element too small\n");
                return;
            }
            for (j = 0; j < n; j++){
                if (j != i)
                    sum -= A[i][j]*y[j];}
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
    
    for (i = 0; i < size; i++)
        y[i] = x[i];
}

//---------------------------------------------
double* sub_vector(int size, const double* x, const double* y)
{
    int i;
    double* d;
    
    d = new double[size+1];
    for (i = 0; i < size; i++)
        d[i] = x[i] - y[i];
    return d;
}
//--------------------------------------------
bool epsilon_comp(int size, const double* x)
{
    const double epsilon = 0.00001;
    int i, count;
    
    count = 0;
    for (i = 0; i < size; i++)
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
 double eps, int maxit,
 int *numit, double *x )
{
    //printf("gauss\n");
    double *dx = (double*) malloc(n*sizeof(double));
    int i,j,k;
    
    //fprintf(stderr,"maxit=%d n=%d\n",maxit,n);
    for(k=0; k<maxit; k++)
    {
        double residual = 0.0;
        double res_sum = 0.0;
        for(i=0; i<n; i++)
        {
            dx[i] = b[i];
            for(j=0; j<n; j++)
                dx[i] -= A[i][j]*x[j];
            dx[i] /= A[i][i];
            x[i] += dx[i];
            res_sum += dx[i]*dx[i]; //( (dx[i] >= 0.0) ? dx[i] : -dx[i]);
        }
        residual = sqrt(res_sum);
        node_temp[i]=x[i];
        Loadvector(n, L,BB, node_temp);
        //printf("Residual at iteration %4d : %.3e\n",k,residual);
        if(residual <= eps) break;
        //std::cout <<sum << std::endl;
        //printf("x[%d] = %2.4f\n", k, sum);
    }
    *numit = k+1;
    //std::cout <<numit << std::endl;
    
    std::cerr << "converged soluion:\n";
    for(i=0; i<n; i++)
        std::cerr << x[i] <<"\n";
    
    free(dx);
}

//---------------------------------------COnJUGATE GRADIEnT
void cg(int n, double **A, double *b, double eps, int maxit,double *x)
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
            // cout<<x[i]<<" ";
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
        //cout<<sqrt(sum)<<"\n";
        if( sqrt(sum) <= eps)
        {
            std::cout <<"n step to converge(CG method):  " << step <<"\n" ;
            for(unsigned int i=0; i<n; i++)
                std::cerr << x[i] <<"\n";
            //std::cout <<"the eps is : " << eps << "\n" ;
            break;
        }
    }
}