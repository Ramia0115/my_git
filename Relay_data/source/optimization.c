//memo
//dRealÅ®float
//dQuaternionÅ®float[4]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "agex.h"
#include "agex_body.h"
#include "cgmres.h"
#include "cgmres_body.h"
#include "rhfuncu.h"
#include "myfunc.h"
#include "obstacle_detection.h"
#include "prototype_declaration.h"
#include "setting.h"

#define MAXLOGGINGNUM 1000000

//x for optimization
const double th_offset[4] = {1.5707963267948966, 1.5707963267948966, -1.5707963267948966, -1.5707963267948966};
double x_body[3] = {1.0,1.0,0.0};
double leg_xy[4][2] = 	{{1+0.22086,	1+0.28651},
			 {1-0.22086,	1+0.28651},
			 {1-0.22086,	1-0.28651},
			 {1+0.22806,	1-0.28651}};
double leg_joint_angle[4][3] = 	{{   MyPI/6,	 7*MyPI/12,	 MyPI/2},
				 { 5*MyPI/6,	 5*MyPI/12,	 MyPI/2},
				 {-5*MyPI/6,	-5*MyPI/12,	-MyPI/2},
				 {  -MyPI/6,	-7*MyPI/12,	-MyPI/2}};
double leg_ref[4][2] = 	{{ 0.22086,	 0.28651},
			 {-0.22086,	 0.28651},
			 {-0.22086,	-0.28651},
			 { 0.22806,	-0.28651}};
double u_body[3];
double u_leg[4][4];
double c_error[DIMC*WHEEL_NUM],c_m[DIMC*WHEEL_NUM];
double u_opt[DIMUC];
double du[DIMUC];
double cgmres_error[5] = {0.0};
double c[dv][DIMC];
double p1_1[DIMX],p2[DIMXB],xreference[DIMPB+DIMP*WHEEL_NUM];

int opt_cnt = 0;

double x_b_opt[dv+1][DIMXB];
double x_l_opt[WHEEL_NUM][dv+1][DIMX];
double u_b_opt[dv][DIMUCB];
double u_l_opt[WHEEL_NUM][dv][DIMUC];

extern double x0_b[DIMXB],u0_b[DIMUCB];
extern double x0[WHEEL_NUM][DIMX],u0[WHEEL_NUM][DIMUC],xref[DIMX*WHEEL_NUM];
extern double sim_dt;

int cnt;
int rows,cols;
double **logging;

int init_optimization(void)
{
	int log_cnt;
	int col_cnt;
	int i,j;

	log_cnt = MAXLOGGINGNUM;

	printf("memory allocation start.\n");
	//memory allocation for logging
	cols = log_cnt+1;
	//time,(velocity,omega*3,leg position*2,joint angle*3)*4,c_error*6*4,|F|*4,time*4
	//rows = 1+4*(1+3+2+3)+6*4+4+4;
	rows = 1+3+4*5+3+4*4+5+DIMC*WHEEL_NUM;
	logging = (double **)malloc(cols*sizeof(double));
	for(i=0;i<cols;i++)
	{
		logging[i] = (double *)malloc(rows*sizeof(double));
	}
	//initialize memory for logging
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			logging[j][i] = 0.0;
		}
	}
	col_cnt = 0;
	cnt = 0;
	logging[cnt][col_cnt++] = 0.0;

	logging[cnt][col_cnt++] = x_body[0];
	logging[cnt][col_cnt++] = x_body[1];
	logging[cnt][col_cnt++] = x_body[2];

	for(i=0;i<4;i++)
	{
		logging[cnt][col_cnt++] = leg_xy[i][0];
		logging[cnt][col_cnt++] = leg_xy[i][1];
		logging[cnt][col_cnt++] = leg_joint_angle[i][0];
		logging[cnt][col_cnt++] = leg_joint_angle[i][1];
		logging[cnt][col_cnt++] = leg_joint_angle[i][2];
	}

	logging[cnt][col_cnt++] = u_body[0];
	logging[cnt][col_cnt++] = u_body[1];
	logging[cnt][col_cnt++] = u_body[2];

	/*
	for(i=0;i<4;i++)
	{
		logging[cnt][col_cnt++] = u_leg[i][0];
		logging[cnt][col_cnt++] = u_leg[i][1];
		logging[cnt][col_cnt++] = u_leg[i][2];
		logging[cnt][col_cnt++] = u_leg[i][3];
	}
	*/

	for(i=0;i<5;i++)
	{
		logging[cnt][col_cnt++] = cgmres_error[i];
	}

	/*
	for(i=0;i<DIMC*WHEEL_NUM;i++)
	{
		logging[cnt][col_cnt++] = c_m[i];
	}
	*/

	cnt++;
	printf("finish.\n");
	
	return 1;
}

int body_optimization(double opt_time)
{
	static double x_b_buff[DIMXB];
	static double u_b_buff[DIMUCB];
	static double x_l_buff[WHEEL_NUM][DIMX];
	static double u_l_buff[WHEEL_NUM][DIMUC];

	static double prex_b[DIMXB];
	static double prex_l[WHEEL_NUM][DIMX];

	static double xb_c_buff[DIMXB];
	static double xl_c_buff[WHEEL_NUM][DIMX];
	static double xb_dot[DIMXB];
	static double xl_dot[WHEEL_NUM][DIMX];

	static double time_tmp = 0.0;
	double tauf,htau;
	double ptime = 0.0;

	int i,j,k;

	//store the current state to buffer
	x_b_buff[0] = x_body[0];
	x_b_buff[1] = x_body[1];
	x_b_buff[2] = x_body[2];
	for(i=0;i<WHEEL_NUM;i++)
	{
		x_l_buff[i][0] = leg_xy[i][0];
		x_l_buff[i][1] = leg_xy[i][1];
		x_l_buff[i][2] = leg_joint_angle[i][0];
		x_l_buff[i][3] = leg_joint_angle[i][1];
		x_l_buff[i][4] = leg_joint_angle[i][2];
	}

	//compute and store the error about constraints using the current state
	//comp_constraint(x_b_buff,x_l_buff,c_m);

	//initialization
	if(!opt_cnt)
	{
		//store the initial state for G/GMRES
		x0_b[0] = x_b_buff[0];
		x0_b[1] = x_b_buff[1];
		x0_b[2] = x_b_buff[2];

		//initialization of variables
		cgmres_init_b(x_b_buff,u_b_buff);

		for(i=0;i<WHEEL_NUM+1;i++)
		{
			cgmres_error[i] = 0.0;
		}

		prex_b[0] = x_b_buff[0];
		prex_b[1] = x_b_buff[1];
		prex_b[2] = x_b_buff[2];

		//initialization of buffer for prediction of sub-system
		for(i=0;i<dv;i++)
		{
			for(j=0;j<DIMUCB;j++)
			{
				u_b_opt[i][j] = u0_b[j];
			}
		}
		for(i=0;i<WHEEL_NUM;i++)
		{
			for(j=0;j<dv;j++)
			{
				for(k=0;k<DIMUC;k++)
				{
					u_l_opt[i][j][k] = u0[i][k];
				}
			}
		}
	}

	//set time
	tauf = tf*(1.0-exp(-alpha*opt_time));
	htau = tauf/(double)dv;

	//compute reference
	//maybe this section is useless
	//in the original, this section is written for logging
	/*
	pfunc_b(opt_time,p2);
	xreference[0] = p2[0];
	xreference[1] = p2[1];
	xreference[2] = p2[2];
	for(i=0;i<WHEEL_NUM;i++)
	{
		pfunc(opt_time,p1_1,i);
		for(j=0;j<DIMX;j++)
		{
			xreference[DIMXB+i*DIMX+j] = p1_1[j];
		}
	}
	*/

	//store the current state as initial state in the prediction
	x_b_opt[0][0] = x_b_buff[0];
	x_b_opt[0][1] = x_b_buff[1];
	x_b_opt[0][2] = x_b_buff[2];
	for(i=0;i<WHEEL_NUM;i++)
	{
		x_l_opt[i][0][0] = x_l_buff[i][0];
		x_l_opt[i][0][1] = x_l_buff[i][1];
		x_l_opt[i][0][2] = x_l_buff[i][2];
		x_l_opt[i][0][3] = x_l_buff[i][3];
		x_l_opt[i][0][4] = x_l_buff[i][4];
	}

	for(i=0;i<dv;i++)
	{
		for(j=0;j<WHEEL_NUM;j++)
		{
			nfeulerinp(xpfunc,opt_time,x_l_opt[j][i],u_l_opt[j][i],htau,DIMX,x_l_opt[j][i+1],i,j);
		}
	}

	cgmres_error[0] = unew_b(opt_time,prex_b,x_b_buff,u_b_buff);

	prex_b[0] = x_b_buff[0];
	prex_b[1] = x_b_buff[1];
	prex_b[2] = x_b_buff[2];

	u_body[0] = u_b_buff[0];
	u_body[1] = u_b_buff[1];
	u_body[2] = u_b_buff[2];

	return 0;
}

void leg_update(int leg_number)
{
	int i;

	leg_xy[leg_number][0] = leg_xy[leg_number][0]+u_leg[leg_number][0]*cos(leg_joint_angle[leg_number][2]-th_offset[leg_number])*sim_dt;
	leg_xy[leg_number][1] += u_leg[leg_number][0]*sin(leg_joint_angle[leg_number][2]-th_offset[leg_number])*sim_dt;
	for(i=0;i<3;i++)
	{
		leg_joint_angle[leg_number][i] += u_leg[leg_number][1+i]*sim_dt;
	}
}

void body_update(void)
{
	x_body[0] = x_body[0]+u_body[0]*sim_dt;
	x_body[1] = x_body[1]+u_body[1]*sim_dt;
	x_body[2] = x_body[2]+u_body[2]*sim_dt;
}

void data_logging(double rt_time)
{
	int i,j;
	int col_cnt;

	//logging
	col_cnt = 0;

	logging[cnt][col_cnt++] = rt_time;

	logging[cnt][col_cnt++] = x_body[0];
	logging[cnt][col_cnt++] = x_body[1];
	logging[cnt][col_cnt++] = x_body[2];

	for(i=0;i<4;i++)
	{
		logging[cnt][col_cnt++] = leg_xy[i][0];
		logging[cnt][col_cnt++] = leg_xy[i][1];
		logging[cnt][col_cnt++] = leg_joint_angle[i][0];
		logging[cnt][col_cnt++] = leg_joint_angle[i][1];
		logging[cnt][col_cnt++] = leg_joint_angle[i][2];
	}

	logging[cnt][col_cnt++] = u_body[0];
	logging[cnt][col_cnt++] = u_body[1];
	logging[cnt][col_cnt++] = u_body[2];

	/*
	for(i=0;i<4;i++)
	{
		logging[cnt][col_cnt++] = u_leg[i][0];
		logging[cnt][col_cnt++] = u_leg[i][1];
		logging[cnt][col_cnt++] = u_leg[i][2];
		logging[cnt][col_cnt++] = u_leg[i][3];
	}
	*/

	for(i=0;i<5;i++)
	{
		logging[cnt][col_cnt++] = cgmres_error[i];
	}

	/*
	for(i=0;i<DIMC*WHEEL_NUM;i++)
	{
		logging[cnt][col_cnt++] = c_m[i];
	}
	*/

	cnt++;
}

int savedata(void)
{
	FILE *fp;
	int i,j;
	
	printf("saving data.\n");
	//save data
	fp = fopen("logging.mat","w");
	for(i=0;i<cnt;i++)
	{
		for(j=0;j<rows;j++)
		{
			fprintf(fp,"%lf\t",logging[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	printf("finished.\n");
	
	return 0;
}
