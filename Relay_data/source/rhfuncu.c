#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>

#include "rhfuncu.h"
#include "agex.h"
#include "my_math.h"

#define TDIV 100  /* Minimal resolution of computation time is 1/TDIV */

/*-----------------------------------------------------------------------
     Some Fundamental Functions for RHC 
                     T. Ohtsuka  '97/10/30~'97/10/31 (rhfunc.c)
								 '00/01/27 (rhfuncu.c)
-----------------------------------------------------------------------*/

/*---- Matrix ----*/

/*---- a[m][n] -> b[m][n] ----*/
void mmov( int m, int n, double *a[], double *b[] )
{
     int	i, j;//, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++){
			  b[i][j] = a[i][j];
               //tmp = i*n + j ;
               //b[tmp] = a[tmp] ;
          }
}

/*---- a[m][n] + b[m][n] -> c[m][n] ----*/
void madd( int m, int n, double *a[], double *b[], double *c[] )
{
     int	i, j;//, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++) {
               c[i][j] = a[i][j] + b[i][j];
               //tmp = i*n + j ;
               //c[tmp] = a[tmp] + b[tmp] ;
          }
}

/*---- a[m][n] - b[m][n] -> c[m][n] ----*/
void msub( int m, int n, double *a[], double *b[], double *c[] )
{
     int	i, j;//, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++) {
               c[i][j] = a[i][j] - b[i][j];
               //tmp = i*n + j ;
               //c[tmp] = a[tmp] - b[tmp] ;
          }
}

/*---- k * a[m][n] -> b[m][n] ----*/
void mmulsc( int m, int n, double *a[], double k, double *b[] )
{
     int	i, j;//, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++){
               b[i][j] = k * a[i][j];
               //tmp = i*n + j ;
               //b[tmp] = k * a[tmp];
          }
}

/*---- a[m][n] / k -> b[m][n] ----*/
void mdivsc( int m, int n, double *a[], double k, double *b[] )
{
     int	i, j;//, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++){
               b[i][j] = a[i][j] / k;
               //tmp = i*n + j ;
               //b[tmp] = a[tmp] / k;
          }
}

/*---- Vector ----*/

/*---- a[m] -> b[m] ----*/
void vmov( int m, double *a, double *b)
{
     int i;
     for(i = 0; i < m; i++)	b[i] = a[i];
}

/*---- a[m] + b[m] -> c[m] ----*/
void vadd( int m, double *a, double *b, double *c )
{
     int i;
	 for(i = 0; i < m; i++)	c[i] = a[i] + b[i];
}

/*---- a[m] - b[m] -> c[m] ----*/
void vsub( int m, double *a, double *b, double *c )
{
     int i;
     for(i = 0; i < m; i++)	c[i] = a[i] - b[i];
}

/*---- k * a[m] -> b[m] ----*/
void vmulsc( int m, double *a, double k, double *b )
{
     int i;
     for(i = 0; i < m; i++)	b[i] = k * a[i];

}

/*---- a[m] / k -> b[m] ----*/
void vdivsc( int m, double *a, double k, double *b )
{
     int i;
	 for(i = 0; i < m; i++)	b[i] = a[i] / k;
}

/*---- Inner Product of a[m] and b[m] ----*/
double	mvinner( int m, double *a, double *b )
{
	int	i;
	double tmp;
	tmp =0.;
	for(i=0; i<m; i++)
		tmp += a[i] * b[i];
	return tmp;
}

/*---- Linear Interpolation of x(t) ----*/

/*---- k Givens Roations used in GMRES ----*/
void givrot(int k, double *c, double *s, double *v)
{
	int i;
	double w1, w2;
	for(i=0; i<k; i++){
		w1 = c[i] * v[i] - s[i] * v[i+1];
		w2 = s[i] * v[i] + c[i] * v[i+1];
		v[i] = w1;
		v[i+1] = w2;
	}
}

#ifdef newer
/*------------------------------------------------------------------------
			ODE.C: nfrkginpex() and nfadamsinp() for NLSF/SCM
                  	T.Ohtsuka  	 '90/09/26
								~'91/07/06 (for UNIX (Sun))
								 '92/10/29 ( adams() )
								 '93/06/25 ode.c
								 '99/12/19 with input u[]
-------------------------------------------------------------------------*/

#define	DIMRK	50	/* Maximum Dimension of the state vector (for the Runge-Kutta_Gill Method) */
#define	DIMAD	50	/* Maximum Dimension of the state vector (for the Adams Method) */
#define	DIMEU	50	/* Maximum Dimension of the state vector (for the Euler Method) */

///*--------------------------------------------------------
// Simultaneous Ordinaly Diferential Equation Subroutine
//	Runge-Kutta_Gill Method	(Original by T.Murayama)
//
//---- Variation of rkg() : dx/dt is also returned. ----
//----  (for Start of the Adams method)  ----
//----------------------------------------------------------*/
//void nfrkginpex(func,x,y,u,h,dim,ans,fxy)
//void (*func)();
//double x,y[],u[],h,ans[],fxy[]; 
//int dim;
//{
//	int i;
//	double fval[DIMRK],k1[DIMRK],k2[DIMRK],k3[DIMRK],k4[DIMRK],
//			yp1[DIMRK],yp2[DIMRK],yp3[DIMRK],
//			q1[DIMRK],q2[DIMRK],q3[DIMRK],
//			c1 = 0.2928932188134528,
//			c2 = 0.1213203435596426,
//			c3 = 0.5857864376269054,
//			c4 = 1.707106781186548,
//			c5 = 4.121320343559643,
//			c6 = 3.414213562373097;
//
///*	fval = (double *)malloc( (size_t)( dim * sizeof(double) ) );
//	k1   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
//	k2   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
//	k3   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
//	k4   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
//*/
//
//	func(x,y,u,fxy);
//	for(i = 0; i < dim; i++)
//	{
//		 k1[i] = h * fxy[i];
//		yp1[i] = y[i] + 0.5 * k1[i];
//		 q1[i] = k1[i];
//	}
//	func(x + 0.5 * h,yp1,u,fval);
//	for(i = 0; i < dim; i++)
//	{
//		 k2[i] = h * fval[i];
//		yp2[i] = yp1[i] + c1 * (k2[i] - q1[i]);;
//		 q2[i] = c2 * q1[i] + c3 * k2[i];
//	}
//	func(x + 0.5 * h,yp2,u,fval);
//	for(i = 0; i < dim; i++)
//	{
//		 k3[i] = h * fval[i];
//		yp3[i] = yp2[i] + c4 * (k3[i] - q2[i]);
//		 q3[i] = -c5 * q2[i] + c6 * k3[i];
//	}
//	func(x + h,yp3,u,fval);
//	for(i = 0; i < dim; i++)
//	{
//		 k4[i] = h * fval[i];
//		ans[i] = yp3[i] + k4[i] / 6.0 - q3[i] / 3.0;
//	}
//}

///*--------------------------------------------------------
// Simultaneous Ordinaly Diferential Equation Subroutine
//	Adams Method (Predictor-Corrector Method)
//		Predictor: Adams-Bashforth
//		Corrector: Adams-Moulton
//				T.Ohtsuka  '92/10/24
//----------------------------------------------------------*/
//
//void nfadamsinp(func,x,y,u,f1,f2,f3,h,dim,ans)
//void (*func)();
//double x,y[],u[],f1[],f2[],f3[],h,ans[];  /* fi: func(x,y) at i-steps ago. */
//int dim;						   /* *Notice! fi are updated. */
//{
//	int i;
//	double y1[DIMAD], fval[DIMAD], cd;
///*	double *y1, *fval, cd;
//
//	fval = (double *)malloc( (size_t)( dim * sizeof(double) ) );
//	y1   = (double *)malloc( (size_t)( dim * sizeof(double) ) ); */
//	cd = h / 24.;
//
///*---- Adams-Bashforth Predictor ----*/
//	func(x,y,u,fval);
//	for(i=0; i<dim; i++)
//	{
//		y1[i] = 55.* fval[i] -59.* f1[i] +37.* f2[i] -9.* f3[i];
//		y1[i] *= cd;
//		y1[i] += y[i];
//		f3[i] = f2[i]; f2[i] = f1[i]; f1[i] = fval[i]; /* Shift */
//	}
//
///*---- Adams-Moulton Corrector ----*/
//	func(x+h,y1,u,fval);
//	for(i=0; i<dim; i++)
//	{
//		y1[i] = 9.* fval[i] +19.* f1[i] -5.* f2[i] + f3[i];
//		y1[i] *= cd;
//		ans[i] = y1[i] + y[i];
//	}
//
///*	free(fval);
//	free(y1); */
//}

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Euler Method (Forward Differentce)
				T.Ohtsuka  '99/12/19
----------------------------------------------------------*/
//与えられた関数を用いて前進差分のオイラー法で時間的に更新をしていく関数
void nfeulerinp(void (*func)(double , double *, double *, double *, int, int), double x, double y[], double u[], double h, int dim, double ans[], int dvnum, int leg_number)
{
	int i;
	double fval[DIMEU];

	func(x,y,u,fval,dvnum,leg_number);

	for(i=0; i<dim; i++)
		ans[i] = y[i] + h * fval[i];
}
#endif

#ifdef original
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
void nfgmres(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int pkmax, double *err)
{
	int i,j,k;
	double rho, nu; 

	static double cvec[kmax+1];
	static double svec[kmax+1];
	static double gvec[kmax+1];
	static double tmpvec[dv*DIMUC];
	static double hmat[kmax+1][kmax+1];
	static double vmat[kmax+1][dv*DIMUC];
	
	axfunc(n, x, tmpvec); 
	
	vsub(n, b, tmpvec, tmpvec); 
	rho = my_sqrt(mvinner(n, tmpvec, tmpvec));
	gvec[0] = rho; 
	for(i=1; i<pkmax+1; i++){
		gvec[i] = 0;
	}
	err[0] = rho;
	
	vdivsc(n, tmpvec, rho, vmat[0]);
	
	for(k=0; k<pkmax; k++){
		axfunc(n, vmat[k], vmat[k+1]); 
		
		for(j=0; j<=k; j++){
			hmat[k][j] = mvinner(n, vmat[j], vmat[k+1]);
			vmulsc(n, vmat[j], hmat[k][j], tmpvec);
			vsub(n, vmat[k+1], tmpvec, vmat[k+1]); 
		}
		
		hmat[k][k+1] = my_sqrt(mvinner(n, vmat[k+1], vmat[k+1]));
		
		if( hmat[k][k+1] != 0){
			vdivsc(n, vmat[k+1], hmat[k][k+1], vmat[k+1]);
		}
		
		givrot(k, cvec, svec, hmat[k]); 
		nu = my_sqrt(hmat[k][k] * hmat[k][k] + hmat[k][k+1] * hmat[k][k+1]);
		if( nu != 0 ){
			cvec[k] = hmat[k][k] / nu;
			svec[k] = - hmat[k][k+1] / nu;
			hmat[k][k] = cvec[k] * hmat[k][k] - svec[k] * hmat[k][k+1];
			hmat[k][k+1] = 0;
			givrot(1, cvec+k, svec+k, gvec+k);
		}
		
		rho = fabs(gvec[k+1]);
		err[k+1] = rho;
		
		//if(rho < 1.0e-6)
		//{
		//	//printf("rho is small:rho=%1.10lf\n",rho);
		//	break;
		//}
	}
	
	for(i=k-1; i>=0; i--){
		for(nu=gvec[i], j=i+1; j<k; j++){
			nu -= hmat[j][i] * cvec[j];
		}
		//if(nu != 0.0)	cvec[i] = nu / hmat[i][i] ; 
		cvec[i] = nu / hmat[i][i] ; 
	}
	for(i=0; i<n; i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j][i] * cvec[j];
		}
		x[i] += nu;
	}
}
/*
void nfgmres(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int pkmax, double *err)
{
	int i,j,k;
	double rho, nu; 

	double cvec[20];
	double svec[20];
	double gvec[20];
	double tmpvec[200];
	double hmat[400];
	double vmat[4000];
	
	axfunc(n, x, tmpvec); 
	
	vsub(n, b, tmpvec, tmpvec); 
	rho = my_sqrt(mvinner(n, tmpvec, tmpvec));
	gvec[0] = rho; 
	for(i=1; i<pkmax+1; i++){
		gvec[i] = 0;
	}
	err[0] = rho;
	
	vdivsc(n, tmpvec, rho, vmat);
	
	for(k=0; k<pkmax; k++){
		axfunc(n, vmat+k*200, vmat+(k+1)*200); 
		
		for(j=0; j<=k; j++){
			hmat[k*20+j] = mvinner(n, vmat+j*200, vmat+(k+1)*200);
			vmulsc(n, vmat+j*200, hmat[k*20+j], tmpvec);
			vsub(n, vmat+(k+1)*200, tmpvec, vmat+(k+1)*200); 
		}
		
		hmat[k*20+(k+1)] = my_sqrt(mvinner(n, vmat+(k+1)*200, vmat+(k+1)*200));
		
		if( hmat[k*20+(k+1)] != 0){
			vdivsc(n, vmat+(k+1)*200, hmat[k*20+(k+1)], vmat+(k+1)*200);
		}
		
		givrot(k, cvec, svec, hmat+k*20); 
		nu = my_sqrt(hmat[k*20+k] * hmat[k*20+k] + hmat[k*20+(k+1)] * hmat[k*20+(k+1)]);
		if( nu != 0 ){
			cvec[k] = hmat[k*20+k] / nu;
			svec[k] = - hmat[k*20+(k+1)] / nu;
			hmat[k*20+k] = cvec[k] * hmat[k*20+k] - svec[k] * hmat[k*20+(k+1)];
			hmat[k*20+(k+1)] = 0;
			givrot(1, cvec+k, svec+k, gvec+k);
		}
		
		rho = fabs(gvec[k+1]);
		err[k+1] = rho;
		
	}

	//for(i=0;i<k-1;i++)	printf("%f\t",hmat[i][i]);
	//printf("\n");
	
	for(i=k-1; i>=0; i--){
		for(nu=gvec[i], j=i+1; j<k; j++){
			nu -= hmat[j*20+i] * cvec[j];
		}
		//if(nu != 0.0)	cvec[i] = nu / hmat[i][i] ; 
		if(err[i] != 0.0)	cvec[i] = nu / hmat[i*20+i] ; 
		else			cvec[i] = 0.0;
		//printf("%lf\t%lf\t%lf\n",nu,hmat[i][i],cvec[i]);
	}
	for(i=0; i<n; i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j*200+i] * cvec[j];
		}
		x[i] += nu;
	}
}
*/
/*------------------------------------------------------------------------
			ODE.C: nfrkginpex() and nfadamsinp() for NLSF/SCM
                  	T.Ohtsuka  	 '90/09/26
								~'91/07/06 (for UNIX (Sun))
								 '92/10/29 ( adams() )
								 '93/06/25 ode.c
								 '99/12/19 with input u[]
-------------------------------------------------------------------------*/

#define	DIMEU	50	/* Maximum Dimension of the state vector (for the Euler Method) */
#define	DIMRK	50	/* Maximum Dimension of the state vector (for the Runge-Kutta_Gill Method) */

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Runge-Kutta_Gill Method	(Original by T.Murayama)

---- Variation of rkg() : dx/dt is also returned. ----
----  (for Start of the Adams method)  ----
----------------------------------------------------------*/
void nfrkginpex(void (*func)(double , double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[])
{
	int i;
	double fval[DIMRK],k1[DIMRK],k2[DIMRK],k3[DIMRK],k4[DIMRK],
			yp1[DIMRK],yp2[DIMRK],yp3[DIMRK],
			q1[DIMRK],q2[DIMRK],q3[DIMRK],
			c1 = 0.2928932188134528,
			c2 = 0.1213203435596426,
			c3 = 0.5857864376269054,
			c4 = 1.707106781186548,
			c5 = 4.121320343559643,
			c6 = 3.414213562373097;

	func(x,y,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k1[i] = h * fval[i];
		yp1[i] = y[i] + 0.5 * k1[i];
		 q1[i] = k1[i];
	}
	func(x + 0.5 * h,yp1,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k2[i] = h * fval[i];
		yp2[i] = yp1[i] + c1 * (k2[i] - q1[i]);;
		 q2[i] = c2 * q1[i] + c3 * k2[i];
	}
	func(x + 0.5 * h,yp2,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k3[i] = h * fval[i];
		yp3[i] = yp2[i] + c4 * (k3[i] - q2[i]);
		 q3[i] = -c5 * q2[i] + c6 * k3[i];
	}
	func(x + h,yp3,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k4[i] = h * fval[i];
		ans[i] = yp3[i] + k4[i] / 6.0 - q3[i] / 3.0;
	}
}

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Euler Method (Forward Differentce)
				T.Ohtsuka  '99/12/19
----------------------------------------------------------*/
//void nfeulerinp(void (*func)(double , double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[])
void nfeulerinp(void (*func)(double , double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[])
//void nfeulerinp(void (*func)(double , double *, double *, double *, double *), double x, double y[], double u[], double h, int dim, double ans[], double nu[])
{
	int i;
	double fval[DIMEU];

	func(x,y,u,fval);
	for(i=0; i<dim; i++)
		ans[i] = y[i] + h * fval[i];
}
#endif