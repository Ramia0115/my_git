#include <math.h>
#include <stdio.h>

#include "rhfuncu.h"
#include "cgmres.h"
#include "cgmres_body.h"
#include "myfunc.h"
#include "my_math.h"
#include "PRS.h"

/*-------------- Global Variables Defined by User -------------- */
double  q_b[DIMXB] = {1.0, 1.0, 1.0};
double  r_b[DIMUCB] = {0.1, 0.1, 0.1};
double  sf_b[DIMXB] = {3.0, 3.0, 3.0 };
//double w[8] = {1., 1., 1., 1., 1., 1., 1., 1.};
double w[8] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

extern double x_l_opt[WHEEL_NUM][dv+1][DIMX];
extern double x0_b[DIMXB];

/*-------------- Time-Variant Parameters --------------*/
void pfunc_b(double t, double p1[])
{
	static double Vx;
	static double angle_rot = 30*MyPI/180;
	static double angle_saw = 20*MyPI/180;
	
	Vx = 0.25*(1 - exp(-0.5*t));

	p1[0] = x0_b[0];
	p1[1] = x0_b[1];
	p1[2] = x0_b[2];

	///////////////////////////////////ê^Ç¡íºÇÆëñÇÈ/////////////////////////////////////////////
	#ifdef STR
	p1[0] += Vx*t;
	p1[1] += 0.0;
	p1[2] += 0.0;	
	#endif

	///////////////////////////////////ï˚å¸äpïœâª/////////////////////////////////////////////
	#ifdef ROT
	p1[0] += Vx*t;
	p1[1] += 0.0;
	p1[2] += angle_rot*sin(0.5*t) * (1 - exp(-0.25*t));
	#endif

	///////////////////////////////////ÉmÉRÉMÉäãOìπ/////////////////////////////////////////////
	#ifdef SAW
	if( t<=4.0 )
	{
		p1[0] += Vx*cos(angle_saw)*t;
		p1[1] += Vx*sin(angle_saw)*t;
		p1[2] += angle_saw;
	}else if( t<=8.0 )
	{
		p1[0] += Vx*cos(angle_saw)*t;
		p1[1] += Vx*sin(angle_saw)*4.0 - Vx*sin(angle_saw)*(t - 4.0);
		p1[2] += -angle_saw;
	}else if( t<=12.0 )
	{
		p1[0] += Vx*cos(angle_saw)*t;
		p1[1] += Vx*sin(angle_saw)*(t - 8.0);
		p1[2] += angle_saw;
	}else if( t<=16.0 )
	{
		p1[0] += Vx*cos(angle_saw)*t;
		p1[1] += Vx*sin(angle_saw)*4.0 - Vx*sin(angle_saw)*(t - 12.0);
		p1[2] += -angle_saw;
	}else
	{
		p1[0] += Vx*cos(angle_saw)*t;
		p1[1] += Vx*sin(angle_saw)*(t - 16.0);
		p1[2] += angle_saw;
	}
	#endif

}


/*-------------- dPhi/dx -------------- */
void phix_b(double t, double x[], double phx1[])
{
	double p[DIMPB];

	pfunc_b(t, p);

	phx1[0] = 1.*sf_b[0] * (-1.*p[0] + x[0]);
	phx1[1] = 1.*sf_b[1] * (-1.*p[1] + x[1]);
	phx1[2] = 1.*sf_b[2] * (-1.*p[2] + x[2]);
}

/*-------------- State Equation -------------- */
void xpfunc_b(double t, double x[], double u[], double xprime[], int dvnum, int i)
{
	xprime[0] = u[0];
	xprime[1] = u[1];
	xprime[2] = u[2];
}

/*-------------- Costate Equation -------------- */
void lpfunc_b(double t, double lmd[], double linp[], double lprime[], int dvnum, int i)
{
	double x[DIMXB], u[DIMUCB];
    double o[134];
	double p[DIMPB];
	double l[4] = {0.1767766953, 0.13, 0.03677, 0.061};
	double beta[WHEEL_NUM] = {MyPIq, 3*MyPIq, -3*MyPIq, -MyPIq};

	vmov(DIMXB, linp, x);
	vmov(DIMUCB, linp + DIMXB, u);

	pfunc_b(t, p);

	o[0] = beta[0];
    o[1] = o[0] + x[2];
    o[2] = cos(o[1]);
    o[3] = l[0]*o[2];
    o[4] = -1.*o[3];
    o[5] = x_l_opt[0][dvnum][2];
    o[6] = cos(o[5]);
    o[7] = l[1]*o[6];
    o[8] = -1.*o[7];
    o[9] = x_l_opt[0][dvnum][3];
    o[10] = cos(o[9]);
    o[11] = l[2]*o[10];
    o[12] = -1.*o[11];
    o[13] = x_l_opt[0][dvnum][4];
    o[14] = cos(o[13]);
    o[15] = l[3]*o[14];
    o[16] = -1.*o[15];
    o[17] = -1.*x[0];
    o[18] = x_l_opt[0][dvnum][0];
    o[19] = o[4] + o[8] + o[12] + o[16] + o[17] + o[18];
    o[20] = beta[1];
    o[21] = o[20] + x[2];
    o[22] = cos(o[21]);
    o[23] = l[0]*o[22];
    o[24] = -1.*o[23];
    o[25] = x_l_opt[1][dvnum][2];
    o[26] = cos(o[25]);
    o[27] = l[1]*o[26];
    o[28] = -1.*o[27];
    o[29] = x_l_opt[1][dvnum][3];
    o[30] = cos(o[29]);
    o[31] = l[2]*o[30];
    o[32] = -1.*o[31];
    o[33] = x_l_opt[1][dvnum][4];
    o[34] = cos(o[33]);
    o[35] = l[3]*o[34];
    o[36] = -1.*o[35];
    o[37] = x_l_opt[1][dvnum][0];
    o[38] = o[17] + o[24] + o[28] + o[32] + o[36] + o[37];
    o[39] = beta[2];
    o[40] = o[39] + x[2];
    o[41] = cos(o[40]);
    o[42] = l[0]*o[41];
    o[43] = -1.*o[42];
    o[44] = x_l_opt[2][dvnum][2];
    o[45] = cos(o[44]);
    o[46] = l[1]*o[45];
    o[47] = -1.*o[46];
    o[48] = x_l_opt[2][dvnum][3];
    o[49] = cos(o[48]);
    o[50] = l[2]*o[49];
    o[51] = -1.*o[50];
    o[52] = x_l_opt[2][dvnum][4];
    o[53] = cos(o[52]);
    o[54] = l[3]*o[53];
    o[55] = -1.*o[54];
    o[56] = x_l_opt[2][dvnum][0];
    o[57] = o[17] + o[43] + o[47] + o[51] + o[55] + o[56];
    o[58] = beta[3];
    o[59] = o[58] + x[2];
    o[60] = cos(o[59]);
    o[61] = l[0]*o[60];
    o[62] = -1.*o[61];
    o[63] = x_l_opt[3][dvnum][2];
    o[64] = cos(o[63]);
    o[65] = l[1]*o[64];
    o[66] = -1.*o[65];
    o[67] = x_l_opt[3][dvnum][3];
    o[68] = cos(o[67]);
    o[69] = l[2]*o[68];
    o[70] = -1.*o[69];
    o[71] = x_l_opt[3][dvnum][4];
    o[72] = cos(o[71]);
    o[73] = l[3]*o[72];
    o[74] = -1.*o[73];
    o[75] = x_l_opt[3][dvnum][0];
    o[76] = o[17] + o[62] + o[66] + o[70] + o[74] + o[75];
    o[77] = sin(o[1]);
    o[78] = l[0]*o[77];
    o[79] = -1.*o[78];
    o[80] = sin(o[5]);
    o[81] = l[1]*o[80];
    o[82] = -1.*o[81];
    o[83] = sin(o[9]);
    o[84] = l[2]*o[83];
    o[85] = -1.*o[84];
    o[86] = sin(o[13]);
    o[87] = l[3]*o[86];
    o[88] = -1.*o[87];
    o[89] = -1.*x[1];
    o[90] = x_l_opt[0][dvnum][1];
    o[91] = o[79] + o[82] + o[85] + o[88] + o[89] + o[90];
    o[92] = sin(o[21]);
    o[93] = l[0]*o[92];
    o[94] = -1.*o[93];
    o[95] = sin(o[25]);
    o[96] = l[1]*o[95];
    o[97] = -1.*o[96];
    o[98] = sin(o[29]);
    o[99] = l[2]*o[98];
    o[100] = -1.*o[99];
    o[101] = sin(o[33]);
    o[102] = l[3]*o[101];
    o[103] = -1.*o[102];
    o[104] = x_l_opt[1][dvnum][1];
    o[105] = o[89] + o[94] + o[97] + o[100] + o[103] + o[104];
    o[106] = sin(o[40]);
    o[107] = l[0]*o[106];
    o[108] = -1.*o[107];
    o[109] = sin(o[44]);
    o[110] = l[1]*o[109];
    o[111] = -1.*o[110];
    o[112] = sin(o[48]);
    o[113] = l[2]*o[112];
    o[114] = -1.*o[113];
    o[115] = sin(o[52]);
    o[116] = l[3]*o[115];
    o[117] = -1.*o[116];
    o[118] = x_l_opt[2][dvnum][1];
    o[119] = o[89] + o[108] + o[111] + o[114] + o[117] + o[118];
    o[120] = sin(o[59]);
    o[121] = l[0]*o[120];
    o[122] = -1.*o[121];
    o[123] = sin(o[63]);
    o[124] = l[1]*o[123];
    o[125] = -1.*o[124];
    o[126] = sin(o[67]);
    o[127] = l[2]*o[126];
    o[128] = -1.*o[127];
    o[129] = sin(o[71]);
    o[130] = l[3]*o[129];
    o[131] = -1.*o[130];
    o[132] = x_l_opt[3][dvnum][1];
    o[133] = o[89] + o[122] + o[125] + o[128] + o[131] + o[132];
    lprime[0] = -0.5*(-2.*o[19]*w[0] - 2.*o[38]*w[2] - 2.*o[57]*w[4] - 2.*o[76]*w[6]) - 1.*q_b[0]*(-1.*p[0] + x[0]);
    lprime[1] = -0.5*(-2.*o[91]*w[1] - 2.*o[105]*w[3] - 2.*o[119]*w[5] - 2.*o[133]*w[7]) - 1.*q_b[1]*(-1.*p[1] + x[1]);
    lprime[2] = -0.5*(2.*l[0]*o[19]*o[77]*w[0] - 2.*l[0]*o[2]*o[91]*w[1] + 2.*l[0]*o[38]*o[92]*w[2] - 2.*l[0]*o[22]*o[105]*w[3] + 2.*l[0]*o[57]*o[106]*w[4] - 2.*l[0]*o[41]*o[119]*w[5] + 2.*l[0]*o[76]*o[120]*w[6] - 2.*l[0]*o[60]*o[133]*w[7]) - 1.*q_b[2]*(-1.*p[2] + x[2]);
}

/*-------------- Error in Optimality Condition, Hu -------------- */
void hufunc_b(double t, double x[], double lmd[], double u[], double hui[], int dvnum)
{
    hui[0] = lmd[0] + 1.*r_b[0]*u[0];
    hui[1] = lmd[1] + 1.*r_b[1]*u[1];
    hui[2] = lmd[2] + 1.*r_b[2]*u[2];
}
