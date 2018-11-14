#ifndef RHFUNCU

#define RHFUNCU

#define newer
/*-----------------------------------------------------------------------
     Some Fundamental Functions for RHC 
                     T. Ohtsuka  '97/10/30~'97/10/31 (rhfunc.c)
								 '00/01/27 (rhfuncu.c)
-----------------------------------------------------------------------*/

/*---- Matrix ----*/

/*---- a[m][n] -> b[m][n] ----*/
void mmov( int m, int n, double *a[], double *b[] );

/*---- a[m][n] + b[m][n] -> c[m][n] ----*/
void madd( int m, int n, double *a[], double *b[], double *c[] );

/*---- a[m][n] - b[m][n] -> c[m][n] ----*/
void msub( int m, int n, double *a[], double *b[], double *c[] );

/*---- k * a[m][n] -> b[m][n] ----*/
void mmulsc( int m, int n, double *a[], double k, double *b[] );

/*---- a[m][n] / k -> b[m][n] ----*/
void mdivsc( int m, int n, double *a[], double k, double *b[] );

/*---- Vector ----*/

/*---- a[m] -> b[m] ----*/
void vmov( int m, double *a, double *b);

/*---- a[m] + b[m] -> c[m] ----*/
void vadd( int m, double *a, double *b, double *c );

/*---- a[m] - b[m] -> c[m] ----*/
void vsub( int m, double *a, double *b, double *c );

/*---- k * a[m] -> b[m] ----*/
void vmulsc( int m, double *a, double k, double *b );

/*---- a[m] / k -> b[m] ----*/
void vdivsc( int m, double *a, double k, double *b );

/*---- Inner Product of a[m] and b[m] ----*/
double	mvinner( int m, double *a, double *b );

/*---- k Givens Roations used in GMRES ----*/
void givrot(int k, double *c, double *s, double *v);

/*--------------------------------------------------------
 Linear Equation Subroutine
	GMRES (Generalized Minimum Residual) 
	Cf. C.T.Kelley: Iterative Methods for Linear and Nonlinear Equations
						T.Ohtsuka  '00/01/24
	axfunc(int n, double *x, double *ax): a function that gives A*x
	n: dim x
	b: the right-hand-side of the equation, A*x = b
	x: initial guess, to be replaced with the solution
	kmax: number of iteration
	err: residual norms, |b - A*x_i| (i=1,...,n+1)
----------------------------------------------------------*/
//void nfgmres(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int pkmax, double *err);
//void nfgmres(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int kmax, double *err);
//void nfgmres(void (*axfunc)(int , double *, double *, int), int n, double *b, double *x, int pkmax, double *err, int leg_number);
/*------------------------------------------------------------------------
			ODE.C: nfrkginpex() and nfadamsinp() for NLSF/SCM
                  	T.Ohtsuka  	 '90/09/26
								~'91/07/06 (for UNIX (Sun))
								 '92/10/29 ( adams() )
								 '93/06/25 ode.c
								 '99/12/19 with input u[]
-------------------------------------------------------------------------*/

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Runge-Kutta_Gill Method	(Original by T.Murayama)

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Euler Method (Forward Differentce)
				T.Ohtsuka  '99/12/19
----------------------------------------------------------*/

#ifdef newer
void nfeulerinp(void (*func)(double , double *, double *, double *, int, int), double x, double y[], double u[], double h, int dim, double ans[], int dvnum, int leg_number);
#endif

#ifdef original
void nfeulerinp(void (*func)(double , double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[]);
//void nfeulerinp(void (*func)(double , double *, double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[], double nu[]);
//void nfrkginpex(void (*func)(double , double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[]);
#endif

#endif