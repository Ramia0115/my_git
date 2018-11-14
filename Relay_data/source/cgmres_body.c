#include<stdio.h>
#include<math.h>

#include "main.h"
#include "cgmres.h"
#include "cgmres_body.h"
#include "rhfuncu.h"
#include "myfunc.h"
#include "my_math.h"
#include "PRS.h"

/*-------------- Global Variagles Defined in Program -------------- */
int		dimeq_b;
double onemhdir_b, hdirbht_b, onemzetahdir_b;
double x0_b[DIMXB];// = {X0, Y0, theta0};//Å¶
double u0_b[DIMUCB] = {0.0, 0.0, 0.0};
double x1s_b[DIMXB], bvec_b[dv*DIMUCB], duvec_b[dv*DIMUCB], dutmp_b[dv*DIMUCB], errvec_b[kmax_b+1];
double htau_b;
double ts_b;
double tauf1_b;
double tau_b[dv+1];
double xtau_b[dv+1][DIMXB], xtau1_b[dv+1][DIMXB], ltau_b[dv+1][DIMXB];
double utau_b[dv][DIMUCB], utau1_b[dv][DIMUCB], hutau_b[dv][DIMUCB], hutau1_b[dv][DIMUCB], hutau2_b[dv][DIMUCB];
double lmd0_b[DIMXB], hu0_b[DIMUCB];

extern double u_b_opt[dv][DIMXB];

/*-------------- Initial Conditions -------------- */
void dhu0func_b(int dimu, double *du0, double *dhu)
{
	double du[DIMUCB], u[DIMUCB], hu[DIMUCB];
	vmulsc(DIMUCB, du0, hdir, du);
	vadd(DIMUCB, u0_b, du, u);
	hufunc_b(tsim0, x0_b, lmd0_b, u, hu, 0);
	vsub(DIMUCB, hu, hu0_b, dhu);
	vdivsc(DIMUCB, dhu, hdir, dhu);
}

void cgmres_initialization_b()
{
    int i, j;
	
	for(i=0;i<DIMXB;i++)		x1s_b[i]		= x0_b[i];
	for(i=0;i<dv*DIMUCB;i++)	bvec_b[i]		= 0.0;
	for(i=0;i<dv*DIMUCB;i++)	duvec_b[i]	= 0.0;
	for(i=0;i<dv*DIMUCB;i++)	dutmp_b[i]	= 0.0;
	for(i=0;i<kmax_b+1;i++)	errvec_b[i]	= 0.0;
	for(i=0;i<dv+1;i++)		tau_b[i]		= 0.0;
	for(i=0;i<dv+1;i++)	for(j=0;j<DIMXB;j++)		xtau_b[i][j]		= x0_b[j];
	for(i=0;i<dv+1;i++)	for(j=0;j<DIMXB;j++)		xtau1_b[i][j]		= x0_b[j];
	for(i=0;i<dv+1;i++)	for(j=0;j<DIMXB;j++)		ltau_b[i][j]		= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUCB;j++)	utau_b[i][j]		= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUCB;j++)	utau1_b[i][j]		= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUCB;j++)	hutau_b[i][j]		= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUCB;j++)	hutau1_b[i][j]	= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUCB;j++)	hutau2_b[i][j]	= 0.0;
}

void nfgmres_init_b(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int pkmax, double *err)
{
	int i,j,k;
	double rho, nu; 

	double cvec[DIMUCB+1];
	double svec[DIMUCB+1];
	double gvec[DIMUCB+1];
	double tmpvec[DIMUCB];
	double hmat[DIMUCB+1][DIMUCB+1];
	double vmat[DIMUCB+1][DIMUCB];
	
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
		//	break;
		//}
	}
	
	for(i=k-1; i>=0; i--){
		for(nu=gvec[i], j=i+1; j<k; j++){
			nu -= hmat[j][i] * cvec[j];
		}
		cvec[i] = nu / hmat[i][i] ; 
	}
	for(i=0; i<n; i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j][i] * cvec[j];
		}
		x[i] += nu;
	}
}

void nfgmres_b(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int pkmax, double *err)
{
	int i,j,k;
	double rho, nu; 

	static double cvec[kmax_b+1];
	static double svec[kmax_b+1];
	static double gvec[kmax_b+1];
	static double tmpvec[dv*DIMUCB];
	static double hmat[kmax_b+1][kmax_b+1];
	static double vmat[kmax_b+1][dv*DIMUCB];

	axfunc(n, x, tmpvec); 
	vsub(n, b, tmpvec, tmpvec); 
	rho = sqrt(mvinner(n, tmpvec, tmpvec));
	gvec[0] = rho; 
	for(i=1; i<pkmax+1; i++){
		gvec[i] = 0;
	}
	err[0] = rho;

	vdivsc(n, tmpvec, rho, vmat[0]);
	for(k=0; k<pkmax; k++){
		axfunc(n, vmat[k], vmat[k+1]); 

		/* Modified Gram-Schmidt */
		for(j=0; j<=k; j++){
			hmat[k][j] = mvinner(n, vmat[j], vmat[k+1]);
			vmulsc(n, vmat[j], hmat[k][j], tmpvec);
			vsub(n, vmat[k+1], tmpvec, vmat[k+1]); 
		}
		hmat[k][k+1] = sqrt(mvinner(n, vmat[k+1], vmat[k+1]));
		
		/* No Breakdown? */
		if( hmat[k][k+1] != 0){
			vdivsc(n, vmat[k+1], hmat[k][k+1], vmat[k+1]);
		}
		else{
			printf("gmress() : breakdown \n");
		}
		
		/* Givens Rotation */
		givrot(k, cvec, svec, hmat[k]); 
		nu = sqrt(hmat[k][k] * hmat[k][k] + hmat[k][k+1] * hmat[k][k+1]);
		if( nu != 0 ){
			cvec[k] = hmat[k][k] / nu;
			svec[k] = - hmat[k][k+1] / nu;
			hmat[k][k] = cvec[k] * hmat[k][k] - svec[k] * hmat[k][k+1];
			hmat[k][k+1] = 0;
			givrot(1, cvec+k, svec+k, gvec+k);
		}
		else printf("nu is zero!\n");
		
		/* Residual Update */
		rho = fabs(gvec[k+1]);
		err[k+1] = rho;
	}
	
	/* Solve hmat * y = gvec (hmat: upper triangular) */
	for(i=k-1; i>=0; i--){
		for(nu=gvec[i], j=i+1; j<k; j++){
			nu -= hmat[j][i] * cvec[j];
		}
		cvec[i] = nu / hmat[i][i] ; 
		/*for(nu=0, j=i+1; j<k; j++){
			nu += hmat[j][i] * cvec[j];
		}
		cvec[i] = (gvec[i] - nu) / hmat[i][i] ;*/
	}
	/* Ans. */
	for(i=0; i<n; i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j][i] * cvec[j];
		}
		x[i] += nu;
	}

}

void cgmres_init_b(double x[], double u[])
{
    int i;

	dimeq_b = dv * DIMUCB;
	onemhdir_b = 1 - hdir / ht;
	hdirbht_b = hdir / ht; 
    onemzetahdir_b = 1 - zeta * hdir; 

    vmov(DIMXB, x0_b, x);
	vmov(DIMUCB, u0_b, u);
    for(i=0; i<dv; i++){
        vmov(DIMUCB, u0_b, utau_b[i]);
    }
	for(i=0; i<dimeq_b; i++){
			duvec_b[i] = 0;
	}
}


/*-------------- Control Update -------------- */
//void errfunc(double t, double x[DIMXB], double u[][DIMUCB], double hu[][DIMUCB])
void errfunc_b(double t, double tauf, double x[DIMXB], double u[][DIMUCB], double hu[][DIMUCB])
{
	int i;
	double taut, linp[DIMXB+DIMUCB];
	
	htau_b = tauf / (double)dv;
	vmov(DIMXB, x, xtau_b[0]);

	for(taut = t, i=0; i < dv; taut += htau_b, i++){
		nfeulerinp(xpfunc_b,taut,xtau_b[i],u[i],htau_b,DIMXB,xtau_b[i+1], i, 0);	
		tau_b[i] = taut;
	}
	tau_b[i] = taut; 

	phix_b(taut, xtau_b[dv], ltau_b[dv]);

	for(i = dv-1; i >= 0; i--){

		vmov(DIMXB, xtau_b[i], linp);
		vmov(DIMUCB, u[i], linp+DIMXB);

		nfeulerinp(lpfunc_b,taut,ltau_b[i+1],linp,-htau_b,DIMXB,ltau_b[i], i, 0);
		taut -= htau_b; 

		hufunc_b(taut, xtau_b[i], ltau_b[i+1], u[i], hu[i], i);
	}

//	printf("hu[0] = %lf, %lf, %lf\n", hu[0][0], hu[0][1], hu[0][2]);
}

void adufunc_b(int n, double *du, double *adu)
{
	int j;
	
	vmulsc(n, du, hdir, dutmp_b);
	for(j=0;j<dv;j++)	vadd(DIMUCB, utau_b[j], &dutmp_b[j*DIMUCB], utau1_b[j]);
	errfunc_b(ts_b, tauf1_b, x1s_b, utau1_b, hutau2_b);
	for(j=0;j<dv;j++)	vsub(DIMUCB, hutau2_b[j], hutau1_b[j], &adu[j*DIMUCB]);
	vdivsc(n, adu, hdir, adu);
}

//#define HDIR_EQ_HT

double unew_b(double t, double x[], double x1[], double u[])
{
	int i,j;
	double x2[DIMXB];
	double hu2 = 0.0;
	double tauf;

	ts_b = t + hdir;

	vmulsc(DIMXB, x, onemhdir_b, x2);
	vmulsc(DIMXB, x1, hdirbht_b , x1s_b);
	vadd(DIMXB, x2, x1s_b, x1s_b);

	tauf = tf * (1.0 - exp(-alpha * t) );
	tauf1_b = tf * (1.0 - exp(-alpha * ts_b) );

	errfunc_b(t, tauf, x, utau_b, hutau_b);
	errfunc_b(ts_b, tauf1_b, x1s_b, utau_b, hutau1_b);

	for(i=0;i<dv;i++)	vmulsc(DIMUCB, hutau_b[i], onemzetahdir_b, &bvec_b[i*DIMUCB]);
	for(i=0;i<dv;i++)	vsub(DIMUCB, &bvec_b[i*DIMUCB], hutau1_b[i], &bvec_b[i*DIMUCB]);
	for(i=0;i<dv;i++)	vdivsc(DIMUCB, &bvec_b[i*DIMUCB], hdir, &bvec_b[i*DIMUCB]);

	nfgmres_b(adufunc_b, dimeq_b, bvec_b, duvec_b, kmax_b, errvec_b);

	for(i=0; i<dv; i++)
	{
		for(j=0;j<DIMUCB;j++)
		{
			utau_b[i][j] += ht * duvec_b[i*DIMUCB+j];
		}
	}

	vmov(DIMUCB, utau_b[0], u);

	//ì¸óÕéûånóÒÇÃäiî[
	for(i=0;i<dv;i++)
	{
		for(j=0;j<DIMUCB;j++)
		{
			u_b_opt[i][j] = utau_b[i][j];
		}
	}

	for(i=0; i<dv; i++)
	{
		for (j = 0; j<DIMUCB; j++)
		{
			 hu2 += hutau_b[i][j] * hutau_b[i][j];

		}
	}

	 return hu2;
}
