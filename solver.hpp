#ifndef solver_hpp
#define solver_hpp

void test_system( int n, double **A, double *b );
void assembly(int n, int nelem);
void Loadvector(int n, double **L, double **BB, double *node_temp);
void jacobi(int n, double **A, double *b, int maxit, int eps, double *x);
void cp_vector(int, double* , const double* );
bool epsilon_comp(int, const double* );
double* sub_vector(int, const double* , const double* );
void run_gauss_seidel_method( int n, double **A, double *b, double eps, int maxit, int *numit, double *x );
void cg( int n, double **A, double *b, double eps, int maxit,double *x);

#endif
