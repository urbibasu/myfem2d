#ifndef solver_hpp
#define solver_hpp

void test_system( int n, double **A, double *b );
void assembly(int n, int nelem,double **L, double *BB, double **new_M);
void Loadvector(int n, double **L, double *BB, double *node_temp,double *new_B);
void jacobi(int n, double **A, double *b, int maxit, int eps, double *x);
void cp_vector(int, double* , const double* );
bool epsilon_comp(int, const double* );
double* sub_vector(int, const double* , const double* );
void run_gauss_seidel_method( int n, double **A, double *b, double eps, int maxit, int *numit, double *x, double **L, double *BB);
void cg( int n, double **A, double *b, double eps, int maxit,double *x);

#endif
