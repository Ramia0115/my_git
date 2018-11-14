#include <math.h>
#include <stdio.h>

#include "PRS.h"
#include "agex_body.h"
#include "rhfuncu.h"
#include "cgmres.h"
#include "myfunc.h"
#include "my_math.h"

/*-------------- Global Variables Defined by User -------------- */

//直進
#ifdef STR
double  q[DIMX] = {1.0, 1.0, 0.01, 0.01, 0.01};
double  r[DIMUC-DIMC] = {0.01, 0.01, 0.01, 0.01};
double  sf[DIMX] = {3.0, 3.0, 0.01, 0.01, 0.01};
#endif

//姿勢角変化
#ifdef ROT
double  q[DIMX] = {1.0, 1.0, 0.1, 0.1, 0.1};
double  r[DIMUC-DIMC] = {0.01, 0.01, 0.01, 0.01};
double  sf[DIMX] = {3.0, 3.0, 0.1, 0.1, 0.1};
#endif

//鋸歯軌道
#ifdef SAW
double  q[5] = {4, 4, 0.1, 0.1, 0.1};
double  r[4] = {3, 0.3, 0.3, 0.3};
double  sf[5] = {6, 6, 0.1, 0.1, 0.1};
double  W[3] = {10,10,10};
double  rhoinv = 0.001;
#endif

//*重心から車輪位置までの相対的なオフセット位置
//1:x, 2:y, 3:theta1, 4:theta2, 5:theta3
double  xref[DIMX*WHEEL_NUM];
double thoff[WHEEL_NUM] = {-MyPIh, -MyPIh, MyPIh, MyPIh};
double l[4] = {0.1767766953, 0.13, 0.03677, 0.061};
double beta[WHEEL_NUM] = {MyPIq, 3*MyPIq, -3*MyPIq, -MyPIq};

extern double x_b_opt[dv+1][DIMXB];

/*------------------------------------------------------------*/


void comp_constraint( double x_b[DIMXB], double x_l[WHEEL_NUM][DIMX], double constraint[DIMC*WHEEL_NUM] )
{

	constraint[0] = x_l[0][0] - l[3]*cos(x_l[0][4]) - l[2]*cos(x_l[0][3]) - l[1]*cos(x_l[0][2]) - l[0]*cos(x_b[2]+beta[0]) - x_b[0]; //CoG-Leg1-X
	constraint[1] = x_l[0][1] - l[3]*sin(x_l[0][4]) - l[2]*sin(x_l[0][3]) - l[1]*sin(x_l[0][2]) - l[0]*sin(x_b[2]+beta[0]) - x_b[1]; //CoG-Leg1-Y
	constraint[2] = x_l[1][0] - l[3]*cos(x_l[1][4]) - l[2]*cos(x_l[1][3]) - l[1]*cos(x_l[1][2]) - l[0]*cos(x_b[2]+beta[1]) - x_b[0]; //CoG-Leg2-X
	constraint[3] = x_l[1][1] - l[3]*sin(x_l[1][4]) - l[2]*sin(x_l[1][3]) - l[1]*sin(x_l[1][2]) - l[0]*sin(x_b[2]+beta[1]) - x_b[1]; //CoG-Leg2-Y
	constraint[4] = x_l[2][0] - l[3]*cos(x_l[2][4]) - l[2]*cos(x_l[2][3]) - l[1]*cos(x_l[2][2]) - l[0]*cos(x_b[2]+beta[2]) - x_b[0]; //CoG-Leg3-X
	constraint[5] = x_l[2][1] - l[3]*sin(x_l[2][4]) - l[2]*sin(x_l[2][3]) - l[1]*sin(x_l[2][2]) - l[0]*sin(x_b[2]+beta[2]) - x_b[1]; //CoG-Leg3-Y
	constraint[6] = x_l[3][0] - l[3]*cos(x_l[3][4]) - l[2]*cos(x_l[3][3]) - l[1]*cos(x_l[3][2]) - l[0]*cos(x_b[2]+beta[3]) - x_b[0]; //CoG-Leg4-X
	constraint[7] = x_l[3][1] - l[3]*sin(x_l[3][4]) - l[2]*sin(x_l[3][3]) - l[1]*sin(x_l[3][2]) - l[0]*sin(x_b[2]+beta[3]) - x_b[1]; //CoG-Leg4-Y

}

void pfunc(double t, double p1[], int i)
{
	double CoG_ref[3];

	//重心の目標軌道を算出
	pfunc_b(t, CoG_ref);

	//脚位置・関節角度に対する目標の計算
	//xref[]：初期姿勢に基づく目標脚位置&関節角度(車輪位置はロボット座標系)
	//CoG_ref[]：重心の目標軌道
	
	p1[0] = xref[0+i*DIMX] * cos(CoG_ref[2]) - xref[1+i*DIMX] * sin(CoG_ref[2]) + CoG_ref[0];
	p1[1] = xref[0+i*DIMX] * sin(CoG_ref[2]) + xref[1+i*DIMX] * cos(CoG_ref[2]) + CoG_ref[1];
	p1[2] = xref[2+i*DIMX];
    p1[3] = xref[3+i*DIMX];
    p1[4] = xref[4+i*DIMX];

}
/*-------------- dPhi/dx -------------- */
void phix(double t, double x[], double phx1[], int i)
{
   	double p[DIMP];

	pfunc(t, p, i);

	phx1[0] = 1.*sf[0] * (-1.*p[0] + x[0]);
	phx1[1] = 1.*sf[1] * (-1.*p[1] + x[1]);
	phx1[2] = 1.*sf[2] * (-1.*p[2] + x[2]);
	phx1[3] = 1.*sf[3] * (-1.*p[3] + x[3]);
	phx1[4] = 1.*sf[4] * (-1.*p[4] + x[4]);

}

/*-------------- State Equation -------------- */
//void xpfunc(double t, double x[], double u[], double xprime[], int dvnum)
void xpfunc(double t, double x[], double u[], double xprime[], int dvnum, int i)
{
    double o[1];

    o[0] = thoff[i] + x[4];
    xprime[0] = cos(o[0])*u[0];
    xprime[1] = sin(o[0])*u[0];
    xprime[2] = u[1];
    xprime[3] = u[2];
    xprime[4] = u[3];

}


/*-------------- Costate Equation -------------- */
void lpfunc(double t, double lmd[], double linp[], double lprime[], int dvnum, int i)
{

	double x[DIMX], u[DIMUC];
    double o[1];
    double p[DIMP];

    vmov(DIMX, linp, x);
	vmov(DIMUC, linp + DIMX, u);

	pfunc(t, p, i);

    o[0] = thoff[i] + x[4];
    lprime[0] = -1.*u[4] - 1.*q[0]*(-1.*p[0] + x[0]);
    lprime[1] = -1.*u[5] - 1.*q[1]*(-1.*p[1] + x[1]);
    lprime[2] = -1.*l[1]*sin(x[2])*u[4] + l[1]*cos(x[2])*u[5] - 1.*q[2]*(-1.*p[2] + x[2]);
    lprime[3] = -1.*l[2]*sin(x[3])*u[4] + l[2]*cos(x[3])*u[5] - 1.*q[3]*(-1.*p[3] + x[3]);
    lprime[4] = -1.*cos(o[0])*lmd[1]*u[0] + lmd[0]*sin(o[0])*u[0] - 1.*l[3]*sin(x[4])*u[4] + l[3]*cos(x[4])*u[5] - 1.*q[4]*(-1.*p[4] + x[4]);
}

/*-------------- Error in Optimality Condition, Hu -------------- */
void hufunc(double t, double x[], double lmd[], double u[], double hui[], int dvnum, int i)
{

	double o[3];

    o[0] = thoff[i] + x[4];
    o[1] = x_b_opt[dvnum][2];
    o[2] = beta[i] + o[1];
    hui[0] = cos(o[0])*lmd[0] + lmd[1]*sin(o[0]) + 1.*r[0]*u[0];
    hui[1] = lmd[2] + 1.*r[1]*u[1];
    hui[2] = lmd[3] + 1.*r[2]*u[2];
    hui[3] = lmd[4] + 1.*r[3]*u[3];
    hui[4] = -1.*l[0]*cos(o[2]) - 1.*l[1]*cos(x[2]) - 1.*l[2]*cos(x[3]) - 1.*l[3]*cos(x[4]) + x[0] - 1.*x_b_opt[dvnum][0];
    hui[5] = -1.*l[0]*sin(o[2]) - 1.*l[1]*sin(x[2]) - 1.*l[2]*sin(x[3]) - 1.*l[3]*sin(x[4]) + x[1] - 1.*x_b_opt[dvnum][1];
}