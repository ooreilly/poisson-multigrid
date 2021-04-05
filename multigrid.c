#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define EPSILON 1e-4
#define START_L 10 //1,2,3,...,9,10
double *multiGrid(const int,double*,double*);
double *multigridVCycle(int,double*,double*,double*);
double *L(const int,const double*,double*);
double *L_sub(const int,const double*,const double*,double*);
double *S(const int,const int,double*,const double*);
double *RB(const int,const int,double*,const double*);
double *init_f(int,double,double*);
double *R(const int,double*);
double *P_sub(const int,const int,double*,const double*);
double *sub(const int,const double*,const double*,double*);
double norm(const int,const double,const double*);
void printMatrix(int,double*);
void checkSolution(double*,int);

int main() {
        int l = START_L;
        int m = 1 << l;  // m = 2^l
        double *u = (double *)calloc((m + 1) * (m + 1), sizeof(double));
        double h = (double)1 / m;  // Step length
        double *f = (double *)calloc((m + 1) * (m + 1),
                                     sizeof(double));  // Matrices to Lu = f
        f = init_f(m, h, f);
        // Intialize values for forcing function
        // Solve for u
        u = multiGrid(l, u, f);
        // printMatrix(m,u);
        if (START_L < 10) checkSolution(u, m);
        free(u);
        free(f);
        return 0;
}
double *multiGrid(const int l, double *u, double *f) {
        // Initialize data
        int m = 1 << l;  // m = 2^l
        double h = (double)1 / m;
        double r;  // Residual
        double *v = (double *)calloc((m + 1) * (m + 1), sizeof(double));
        r = norm(m, h, f);
        double norm_;
        // Norm of initial residual
        // Iterate until solution converges
        do {
                u = multigridVCycle(l, u, f, v);
                v = L_sub(m, u, f, v);
                norm_ = norm(m, h, v);
        } while (norm_ > r * EPSILON);
        free(v);
        return u;
}
/*
Recursive multigrid function
*/
double *multigridVCycle(int l, double *u, double *f, double *v) {
        int m = 1 << l;            // m = 2^l
        int m_lMinus1 = m >> 1;    // m^(l-1)
        double h = (double)1 / m;  // Step length
        double *v_lMinus1 =
            (double *)calloc((m_lMinus1 + 1) * (m_lMinus1 + 1), sizeof(double));
        double *d = (double *)calloc((m + 1) * (m + 1), sizeof(double));
        // Solve the linear equation and return the solution
        if (l == 1) {
                u[4] = 0.25 * f[4];
                return u;
        }
        // Perform a Gauss-Seidel sweep, smoothing out the high frequencies
        u = S(l, m, u, f);
        if (l < 3) u = S(l, m, u, f);
        d = L_sub(
            m, u, f,
            d);  // Compute residual, d = Lu -f //d = sub(m,L(m, u, v),f,d);
        d = R(m_lMinus1, d);  // Restrict d to coarser grid, d^(l-1) = Rd^l
        // Perform a new multigrid step on a coarser grid, return v^(l-1)
        v_lMinus1 = multigridVCycle(l - 1, v_lMinus1, d, v);
        // Prolongate and subract from solution, u^l = u^l - Pv^(l-1)
        u = P_sub(m, m_lMinus1, u, v_lMinus1);
        // Perfom a new Gauss-Seidel sweep
        u = S(l, m, u, f);
        if (l != 10) u = S(l, m, u, f);
        // Clean up
        free(v_lMinus1);
        free(d);
        return u;
}

/*
Compute the L1 Matrix Norm
*/
double norm(const int m, const double h, const double *matrix) {
        int i, j;
        double sum = 0;
        double dA = h * h;  // Area element
        for (i = 0; i < (m + 1) * (m + 1); i++) {
                sum += dA * fabs(matrix[i]);
        }
        return sum;
}

/*
Matrix subtraction (u-v)
in: m - size of matrix (m+1)x(m+1)
u - matrix to operate on
v - matrix to subtract
d - matrix to store result in
out: d
*/
double *sub(const int m, const double *u, const double *v, double *d) {
        int i, j;
        // The boundaries are set to zero
        for (i = 0; i < (m + 1) * (m + 1); i++) {
                d[i] = u[i] - v[i];
        }
        return d;
}

/*
Forcing function Lu = f
in: m - size of matrix (m+1)x(m+1)
h - step length
f - matrix to store result in
out: f
*/
double *init_f(int m, double h, double *f) {
        int i, j;
        for (i = 1; i < m; i++) {
                for (j = 1; j < m; j++) {
                        // h is from discretization of operator L
                        f[i * (m + 1) + j] = h * h * sin(PI * h * i) *
                                             sin(PI * h * j) *
                                             sin(sqrt(2) * PI * h * i) *
                                             sin(sqrt(3) * PI * h * j);
                }
        }
        return f;
}

/*
Difference operator L
in: m - size of matrix (m+1)x(m+1)
u - matrix to operate on
v - matrix to store result in , v = Lu
out: v
*/
double *L(const int m, const double *u, double *v) {
        int i, j;
        for (i = 0; i < m; i++) {
                for (j = 0; j < m; j++) {
                        v[i * (m + 1) + j] =
                            4 * u[i * (m + 1) + j] - u[(i - 1) * (m + 1) + j] -
                            u[(i + 1) * (m + 1) + j] - u[i * (m + 1) + j - 1] -
                            u[i * (m + 1) + j + 1];
                        v[j] = 0;
                        if (i == m - 1) v[i * m + j] = 0;
                }
                v[i * m + j] = 0;
                v[i * m] = 0;
        }
        return v;
}
/*
Difference and subraction operator L_sub
in: m - size of matrix (m+1)x(m+1)
u - matrix that L operate on
w - matrix to subtact on
v - matrix to store result in , v = Lu - w
out: v
*/
double *L_sub(const int m, const double *u, const double *w, double *v) {
        int i, j;
        for (i = 1; i < m; i++) {
                for (j = 1; j < m; j++) {
                        v[i * (m + 1) + j] =
                            4 * u[i * (m + 1) + j] - u[(i - 1) * (m + 1) + j] -
                            u[(i + 1) * (m + 1) + j] - u[i * (m + 1) + j - 1] -
                            u[i * (m + 1) + j + 1] - w[i * (m + 1) + j];
                }
        }
        return v;
}
/*
Gauss-Seidel Sweep - Removes high frequency errors
in: m - size of matrix (m+1)x(m+1)
u - matrix to operate on
f - forcing matrix
out: u
*/
double *S(const int l, const int m, double *u, const double *f) {
        int i, j, n;
        for (i = 1; i < m; i++) {
                n = i * m;
                for (j = 1; j < m; j++) {
                        u[n + i + j] =
                            0.25 * (u[n - m + i - 1 + j] +
                                    u[n + m + i + 1 + j] + u[n + i + j - 1] +
                                    u[n + i + j + 1] + f[n + i + j]);
                }
        }
        return u;
}

/*
Restriction operator - Restricts a matrix to a coarser grid (d^(l-1) = Rd^l)
in: m_lMinus1 - m^(l-1)
d - matrix to operate on
out: d - d^(l-1)
*/
double *R(const int m_lMinus1, double *d) {
        int i, j;
        for (i = 0; i < (m_lMinus1 + 1); i++) {
                for (j = 0; j < (m_lMinus1 + 1); j++) {
                        d[i * (m_lMinus1 + 1) + j] =
                            4 * d[2 * i * (2 * m_lMinus1 + 1) + 2 * j];
                }
        }
        return d;
}
/*
Prolongation and subraction operator - Prolongs a matrix to a finer grid
(v^l = Pv^(l-1) and subtrct that value on u^l, u^l = u^l - v^l)
in: m - size of matrix (m+1)x(m+1)
m_lMinus1 - m^(l-1)
u - matrix to use the P operator on to store result in after subtraction
v_lMinus1 - matrix to operate on
out: u - v^l
*/
double *P_sub(const int m, const int m_lMinus1, double *u,
              const double *v_lMinus1) {
        int i, j;
        for (i = 0; i < m_lMinus1; i++) {
                for (j = 0; j < m_lMinus1; j++) {
                        u[2 * i * (2 * m_lMinus1 + 1) + 2 * j] -=
                            v_lMinus1[i * (m_lMinus1 + 1) + j];
                        u[(2 * i + 1) * (2 * m_lMinus1 + 1) + 2 * j] -=
                            0.5 * (v_lMinus1[i * (m_lMinus1 + 1) + j] +
                                   v_lMinus1[(i + 1) * (m_lMinus1 + 1) + j]);
                        u[(2 * i) * (2 * m_lMinus1 + 1) + 2 * j + 1] -=
                            0.5 * (v_lMinus1[i * (m_lMinus1 + 1) + j] +
                                   v_lMinus1[i * (m_lMinus1 + 1) + j + 1]);
                        u[(2 * i + 1) * (2 * m_lMinus1 + 1) + 2 * j + 1] -=
                            0.25 *
                            (v_lMinus1[i * (m_lMinus1 + 1) + j] +
                             v_lMinus1[i * (m_lMinus1 + 1) + j + 1] +
                             v_lMinus1[(i + 1) * (m_lMinus1 + 1) + j] +
                             v_lMinus1[(i + 1) * (m_lMinus1 + 1) + j + 1]);
                }
        }
        return u;
}

/*
Prints a matrix
in: m - size of matrix (m+1)x(m+1)
u - matrix to print
*/
void printMatrix(int m, double *u) {
        int i, j;
        for (i = 0; i < (m + 1); i++) {
                for (j = 0; j < (m + 1); j++)
                        printf(" %e ", u[i * (m + 1) + j]);
                printf("\n");
        }
}
/*
Checks that the solutions is valid
in: m - size of matrix (m+1)x(m+1)
u - matrix to check
*/
void checkSolution(double *u, int m) {
        double h = (double)1 / m;
        int i, j, N = m + 1, K = m + 1;
        double temp;
        double *grid = (double *)calloc((m + 1) * (m + 1), sizeof(double));
        FILE *myfile;
        char filnamn[] = "mg_?.dat";
        filnamn[3] = 48 + START_L;
        printf("Reading gridfunction of size %d x %d from file %s \n", N, K,
               filnamn);
        myfile = fopen(filnamn, "r");
        if (myfile == NULL) {
                printf("File not found!\n");
        }
        for (i = 0; i < N * K; i++) {
                fscanf(myfile, "%le\n", &grid[i]);
        }
        if (norm(m, h, sub(m, u, grid, grid)) < EPSILON)
                printf("Solution OK!\n");
        else
                printf("Solution not OK!\n");
        fclose(myfile);
}
