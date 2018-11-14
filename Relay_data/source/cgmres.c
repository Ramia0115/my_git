#include<stdio.h>
#include<math.h>

#include "main.h"
#include "agex_body.h"
#include "cgmres.h"
#include "rhfuncu.h"
#include "myfunc.h"
#include "my_math.h"
#define MyPI	3.141592653589793
#define MyPIh	1.570796326794897
#define MyPIq	0.785398163397448

/*-------------- Global Variagles Defined in Program -------------- */
int		dimeq;
double onemhdir, hdirbht, onemzetahdir;

//X方向直進
double x0[WHEEL_NUM][DIMX];// =
//{
//	{ 1+0.22806	, 1+0.28651	, MyPI/6	, 7*MyPI/12		, MyPI/2},
//	{ 1-0.22806	, 1+0.28651	, 5*MyPI/6	, 5*MyPI/12		, MyPI/2},
//	{ 1-0.22806	, 1-0.28651	, -5*MyPI/6	, -5*MyPI/12	, -MyPI/2},
//	{ 1+0.22806	, 1-0.28651	, -MyPI/6	, -7*MyPI/12	, -MyPI/2}
//};

double u0[WHEEL_NUM][DIMUC] = 
{
	{ 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0},
	{ 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0},
	{ 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0},
	{ 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0}
};

double x1s[DIMX], bvec[dv*DIMUC], duvec[dv*DIMUC], dutmp[dv*DIMUC], errvec[kmax+1];
double htau;
double ts;
double tauf1;
double tau[dv+1];
double xtau[dv+1][DIMX], xtau1[dv+1][DIMX], ltau[dv+1][DIMX];
double utau[4][dv][DIMUC], utau1[4][dv][DIMUC], hutau[dv][DIMUC], hutau1[dv][DIMUC], hutau2[dv][DIMUC],test[dv+1][15];
double lmd0[DIMX], hu0[DIMUC];

extern double xref[dv+1][DIMX], uref[dv][DIMUC];
extern double u_l_opt[WHEEL_NUM][dv][DIMUC];

#ifdef newer
/*-------------- Initial Conditions -------------- */
void dhu0func(int dimu, double *du0, double *dhu, int leg_number)
{
	double du[DIMUC], u[DIMUC], hu[DIMUC];
	vmulsc(DIMUC, du0, hdir, du);
	vadd(DIMUC, u0[leg_number], du, u);
	hufunc(tsim0, x0[leg_number], lmd0, u, hu, dv-1, leg_number);
	vsub(DIMUC, hu, hu0, dhu);
	vdivsc(DIMUC, dhu, hdir, dhu);
}
#endif

void cgmres_initialization()
{
    int i, j, k;
	
	for(i=0;i<DIMX;i++) x1s[i]	= x0[0][i];
	for(i=0;i<dv*DIMUC;i++)	bvec[i]		= 0.0;
	for(i=0;i<dv*DIMUC;i++)	duvec[i]	= 0.0;
	for(i=0;i<dv*DIMUC;i++)	dutmp[i]	= 0.0;
	for(i=0;i<kmax+1;i++)	errvec[i]	= 0.0;
	for(i=0;i<dv+1;i++)		tau[i]		= 0.0;
	for(i=0;i<dv+1;i++)	for(j=0;j<DIMX;j++)		xtau[i][j]		= x0[0][j];
	for(i=0;i<dv+1;i++)	for(j=0;j<DIMX;j++)		xtau1[i][j]		= x0[0][j];
	for(i=0;i<dv+1;i++)	for(j=0;j<DIMX;j++)		ltau[i][j]		= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUC;j++)	for(k=0;k<WHEEL_NUM;k++) utau[k][i][j]		= u0[k][j];
	for(i=0;i<dv;i++)	for(j=0;j<DIMUC;j++)	for(k=0;k<WHEEL_NUM;k++) utau1[k][i][j]		= u0[k][j];
	for(i=0;i<dv;i++)	for(j=0;j<DIMUC;j++)	hutau[i][j]		= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUC;j++)	hutau1[i][j]	= 0.0;
	for(i=0;i<dv;i++)	for(j=0;j<DIMUC;j++)	hutau2[i][j]	= 0.0;

}

//void nfgmres(void (*axfunc)(int , double *, double *), int n, double *b, double *x, int kmax, double *err)
void nfgmres(void (*axfunc)(int , double *, double *, int), int n, double *b, double *x, int pkmax, double *err, int leg_number)
{
	int i,j,k;
	double rho, nu; 

	static double cvec[kmax+1];
	static double svec[kmax+1];
	static double gvec[kmax+1];
	static double tmpvec[dv*DIMUC];
	static double hmat[kmax+1][kmax+1];
	static double vmat[kmax+1][dv*DIMUC];

	axfunc(n, x, tmpvec,leg_number); 
	vsub(n, b, tmpvec, tmpvec); 
	rho = sqrt(mvinner(n, tmpvec, tmpvec));
	gvec[0] = rho; 
	for(i=1; i<pkmax+1; i++){
		gvec[i] = 0;
	}
	err[0] = rho;

	vdivsc(n, tmpvec, rho, vmat[0]);

	for(k=0; k<pkmax; k++){
		axfunc(n, vmat[k], vmat[k+1],leg_number); 

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
	}
		/* Ans. */
	for(i=0; i<n; i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j][i] * cvec[j];
		}
		x[i] += nu;
	}

}

void nfgmres_init(void (*axfunc)(int , double *, double *, int), int n, double *b, double *x, int pkmax, double *err, int leg_number)
{
	int i,j,k;
	double rho, nu; 

	double cvec[DIMUC+1];
	double svec[DIMUC+1];
	double gvec[DIMUC+1];
	double tmpvec[DIMUC];
	double hmat[DIMUC+1][DIMUC+1];
	double vmat[DIMUC+1][DIMUC];

	axfunc(n, x, tmpvec, leg_number); 
	vsub(n, b, tmpvec, tmpvec); 
	rho = sqrt(mvinner(n, tmpvec, tmpvec));
	gvec[0] = rho; 
	for(i=1; i<pkmax+1; i++){
		gvec[i] = 0;
	}
	err[0] = rho;

	vdivsc(n, tmpvec, rho, vmat[0]);
	for(k=0; k<pkmax; k++){
		axfunc(n, vmat[k], vmat[k+1], leg_number); 

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
/*		for(nu=0, j=i+1; j<k; j++){
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

void init(double x[], double u[], int leg_number)
{
    int i, j;

	dimeq = dv * DIMUC;
	onemhdir = 1 - hdir / ht;
	hdirbht = hdir / ht; 
    onemzetahdir = 1 - zeta * hdir; 

    vmov(DIMX, x0[leg_number], x);
	vmov(DIMUC, u0[leg_number], u);

	for(i=0;i<dv;i++)
	{
		vmov(DIMUC, u0[leg_number], utau[leg_number][i]);
	}

	for(i=0; i<dimeq; i++){
			duvec[i] = 0;
	}
}

void errfunc(double t, double x[DIMX], double u[][DIMUC], double hu[][DIMUC], int leg_number)
{
     int i;
     double tauf, taut, linp[DIMX+DIMUC];

     tauf = tf * (1. - exp(-alpha * t) ); 
     htau = tauf / (double)dv;
	 vmov(DIMX, x, xtau[0]);  

	for(taut = t, i=0; i < dv; taut += htau, i++)
	{
		nfeulerinp(xpfunc,taut,xtau[i],u[i],htau,DIMX,xtau[i+1],i,leg_number);//終端方向への計算→（以前のu[]を使ってxtau[]を求める）
		tau[i] = taut; 
	}

	tau[i] = taut; 
	phix(taut, xtau[dv], ltau[dv], leg_number);//終端の計算（折り返し地点）（ltau[dv]が求まる）
	
	//初期方向へ計算←
	for(i = dv-1; i >= 0; i--)
	{
		vmov(DIMX, xtau[i], linp);
		vmov(DIMUC, u[i], linp+DIMX);
		nfeulerinp(lpfunc,taut,ltau[i+1],linp,-htau,DIMX,ltau[i],i,leg_number);//ltau[]が求まる
		taut -= htau; 
		hufunc(taut, xtau[i], ltau[i+1], u[i], hu[i], i, leg_number);//Hu[]を求める
	}
}

#ifdef newer
void adufunc(int n, double *du, double *adu, int leg_number)
{
	vmulsc(n, du, hdir, dutmp);
	vadd(n, utau[leg_number][0], dutmp, utau1[leg_number][0]);
	errfunc(ts, x1s, utau1[leg_number], hutau2, leg_number);
	vsub(n, hutau2[0], hutau1[0], adu);
	vdivsc(n, adu, hdir, adu);
}
#endif

double unew(double t, double x[], double x1[], double u[], int leg_number)//x:以前の状態,x1:現在状態
{
	int i, j,k;
	double hu2 = 0.0;
	double taut;
	double x2[DIMX];

	ts = t + hdir;

	vmulsc(DIMX, x, onemhdir, x2);//x2=x*onemhdir
	vmulsc(DIMX, x1, hdirbht , x1s);//x1s=x1*hdirbht
	vadd(DIMX, x2, x1s, x1s);//x1s=x2+x1s

	errfunc(t, x, utau[leg_number], hutau, leg_number);
	errfunc(ts, x1s, utau[leg_number], hutau1, leg_number);

	vmulsc(dimeq, hutau[0], onemzetahdir, bvec);
	vsub(dimeq, bvec, hutau1[0], bvec);
	vdivsc(dimeq, bvec, hdir, bvec);

	nfgmres(adufunc, dimeq, bvec, duvec, kmax, errvec,leg_number);

	for(i=0; i<dv; i++){
		for(j=0;j<DIMUC;j++)
		{
			utau[leg_number][i][j] += ht * duvec[i*DIMUC+j];
		}
	}
	vmov(DIMUC, utau[leg_number][0], u);

	//入力時系列群をバッファする
	for(i=0;i<dv;i++)
	{
		for(j=0;j<DIMUC;j++)
		{
			u_l_opt[leg_number][i][j] = utau[leg_number][i][j];
		}
	}

	for(i=0; i<DIMUC; i++) 
	 {
		 for(j=0; j<dv; j++)
		 {
			 hu2 += hutau[j][i] * hutau[j][i];
		 }
	 }

	 return hu2;

}
