#include <iostream>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "fields.hpp"
#include "solver.hpp"


void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.temperature = new double_vec(n);

    var.shpdx = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);
}



void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{
    int n,nelem,i;
    printf("no of nodes : ");
    scanf("%d",&n);
    printf("no of elements : ");
    scanf("%d",&nelem);
    /* Declaring all the global variables*/
    double **new_M;
    double **L;
    double *new_B;
    new_B = (double*) malloc(n*sizeof(double));
    new_M = (double**) malloc(n*sizeof(double*));
    for(i=0; i<n; i++)
        new_M[i] = (double*) malloc(n*sizeof(double));
    L = (double**) malloc(n*sizeof(double*));
    for(i=0; i<n; i++)
        L[i] = (double*) malloc(n*sizeof(double));
    double *BB;
    BB = (double*) malloc(n*sizeof(double));
    double *node_temp;
    node_temp = (double*) malloc(n*sizeof(double));
    for (i=0;i<n;i++)
    {
        node_temp[i] = 32; //initial conditions at all the nodes is 32F//
    }
    
    /* calling the assembly function to create the K matrix and B vector*/
    assembly(n,nelem,L,BB,new_M);
    //Loadvector(n,L,BB,node_temp,new_B);
    
    double eps = 1.0e-5;
    int maxit = 10000; //2*n*n;
    int cnt = 0;
    
    //jacobi(n,A,b,maxit,eps,x);
    run_gauss_seidel_method(n,new_M,new_B,eps,maxit,&cnt,node_temp,L,BB);
    printf("computed %d iterations\n",cnt);
    //cg(n,A,b,eps,maxit,x);
    
    // Free arrays.
#if 0
    for (i = 0; i < n; i++) {
        free (new_M[i]);
    }
#endif
    free(new_M);
    free(new_B);
    free(node_temp);
    
    return;
}