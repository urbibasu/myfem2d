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
    int n,i;
        printf("give the dimension : ");
        scanf("%d",&n);
    double *b;
    b = (double*) malloc(n*sizeof(double));
    double **A;
    A = (double**) malloc(n*sizeof(double*));
    for(i=0; i<n; i++)
        A[i] = (double*) malloc(n*sizeof(double));
    test_system(n,A,b);
    double *x;
    x = (double*) malloc(n*sizeof(double));
    /* we start at an array of all zeroes */
    for(i=0; i<n; i++) x[i] = 0.0;
    double eps = 1.0e-5;
    int maxit = 2*n*n;
    int cnt = 0;
    //jacobi(n,A,b,eps,maxit,x);
    //run_gauss_seidel_method(n,A,b,eps,maxit,&cnt,x);
    cg(n,A,b,eps,maxit,x);
    printf("computed %d iterations\n",cnt);
    double sum = 0.0;
    for(i=0; i<n; i++) /* compute the error */
    {
        double d = x[i] - 1.0;
        sum += (d >= 0.0) ? d : -d;
    }
    printf("error : %.3e\n",sum);
    //}
    for (i = 0; i < n; i++);
    free (A[i]);
    free (A);
    free(b);
    free(x);
    return;
}

