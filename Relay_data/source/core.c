#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <termios.h>
#include <errno.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>

#include "prototype_declaration.h"
#include "setting.h"

extern int cnt;
extern int FlagUDP;
extern int FlagSCI;
extern int StartFlag;
extern int FlagRcv;
extern int sendcnt;
extern int sendcnt_sci;
extern double *buf_tx;
extern double *buf_rx;
extern float *buf_tx_sci;

int FlagMain = 0;
double sim_dt = 0.018;
double leg_xy[4] = {0.0},angle[4][3] = {0.0};
int cols,rows;
double **logging_snd;
int cnt;

int main (void)
{
	double sim_time = 0.0,simulation = 20.0;
	int sim_cnt;
	int check;
	int error_num;
	double *top_addr;
	float *top_addrf;
	int legnum_tmp;

	struct timespec t;
	double start,opt_start;
	double rt_time;
	int sleep_cnt;

	int i,j;

	printf("Initialization start.\n");
	check = initialization();
	if(check == 1)
	{
		FlagMain = 1;
		printf("Initialization succeeded.\n");
	}

	/*
	while(!StartFlag)
	{
		usleep(1000);
	}
	*/

	printf("Start main-loop.\n");
	clock_gettime(CLOCK_MONOTONIC,&t);
	t.tv_nsec += OPTIMAL_CYCLE;
	start = (double)t.tv_sec*NSEC_PER_SEC+(double)t.tv_nsec;
	while(FlagMain)
	{
		clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&t,NULL);
		clock_gettime(CLOCK_MONOTONIC,&t);
		opt_start = (double)t.tv_sec*NSEC_PER_SEC+(double)t.tv_nsec;
		t.tv_nsec += OPTIMAL_CYCLE;
		while(t.tv_nsec >= NSEC_PER_SEC)
		{
			t.tv_nsec -= NSEC_PER_SEC;
			t.tv_sec++;
		}

		rt_time = (opt_start-start)/NSEC_PER_SEC;

		if(FlagRcv == 1)
		{
			top_addr = buf_rx;
			for(i=0;i<4;i++)
			{
				leg_xy[i] = *buf_rx++;
				angle[i][0] = *buf_rx++;
				angle[i][1] = *buf_rx++;
				angle[i][2] = *buf_rx++;
			}
			buf_rx = top_addr;

			top_addrf = buf_tx_sci;
			for(i=0;i<4;i++)
			{
				*buf_tx_sci++ = (float)leg_xy[i];
				*buf_tx_sci++ = (float)angle[i][0];
				*buf_tx_sci++ = (float)angle[i][1];
				*buf_tx_sci++ = (float)angle[i][2];
			}
			buf_tx_sci = top_addrf;
			sendcnt_sci++;
			data_logging(rt_time);
		}
	}

	KB_close();
	check = close_serial_port();
	savedata();

	return 0;
}

int initialization(void)
{
	int flag_udp,flag_kb,flag_sci,flag_log;
	int error_num;
	int check;

	FlagUDP = 1;
	flag_udp = Init_UDP();
	sleep(1);

	FlagSCI = 1;
	flag_sci = open_serial_port();

	flag_kb =1;

	flag_log = init_logging();

	return flag_udp*flag_sci*flag_kb*flag_log;
}

void print_error(int error_code)
{
	printf("print error string : %s\n",strerror(error_code));
	printf("print error code : %d\n",error_code);
}

int init_logging(void)
{
	int log_cnt;
	int col_cnt;
	int i,j;

	log_cnt = MAXLOGGINGNUM;

	printf("Memory allocation start.\n");
	rows = log_cnt+1;
	cols = 4*4+1;
	logging_snd = (double **)malloc(rows*sizeof(double));
	for(i=0;i<rows;i++)
	{
		logging_snd[i] = (double *)malloc(cols*sizeof(double));
	}
	for(i=0;i<cols;i++)
	{
		for(j=0;j<rows;j++)
		{
			logging_snd[j][i] = 0.0;
		}
	}
	col_cnt = 0;
	cnt = 0;

	logging_snd[cnt][col_cnt++] = 0.0;

	for(i=0;i<4;i++)
	{
		logging_snd[cnt][col_cnt++] = 0.0;
		logging_snd[cnt][col_cnt++] = 0.0;
		logging_snd[cnt][col_cnt++] = 0.0;
		logging_snd[cnt][col_cnt++] = 0.0;
	}

	printf("finish.\n");

	return 1;
}

void data_logging(double rt_time)
{
	int i,j;
	int col_cnt;

	col_cnt = 0;

	logging_snd[cnt][col_cnt++] = rt_time;

	for(i=0;i<4;i++)
	{
		logging_snd[cnt][col_cnt++] = leg_xy[i];
		logging_snd[cnt][col_cnt++] = angle[i][0];
		logging_snd[cnt][col_cnt++] = angle[i][1];
		logging_snd[cnt][col_cnt++] = angle[i][2];
	}

	cnt++;
}

void savedata(void)
{
	FILE *fp;
	int i,j;

	printf("saving data.\n");
	fp = fopen("logging.mat","w");
	for(i=0;i<cnt;i++)
	{
		for(j=0;j<cols;j++)
		{
			fprintf(fp,"%lf\t",logging_snd[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	printf("finish,\n");
}
