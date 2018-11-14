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
#include "agex.h"
#include "agex_body.h"
#include "cgmres.h"
#include "cgmres_body.h"
#include "setting.h"

extern double x_body[3];
extern double leg_xy[4][2];
extern double leg_joint_angle[4][3];
extern double u_b_opt[dv][DIMUCB];
extern double u_l_opt[WHEEL_NUM][dv][DIMUC];
extern int cnt;
extern int opt_cnt;
extern int FlagUDP;
extern int StartFlag;
extern int FlagRcv;
extern int sendcnt;
extern double *buf_tx;
extern double *buf_rx;

extern pthread_t ThreadUDPrcv_sq;

int FlagMain = 0;
double sim_dt = 0.018;

int main (void)
{
	double sim_time = 0.0,simulation = 20.0;
	int sim_cnt;
	int check;
	int error_num;
	double *top_addr;
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

	while(!StartFlag)
	{
		usleep(1000);
	}

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

		if(rt_time > 20.0 || !StartFlag)
		{
			FlagUDP = 0;
			FlagMain = 0;
			break;
		}

		top_addr = buf_tx;
		*buf_tx++ = x_body[0];
		*buf_tx++ = x_body[1];
		*buf_tx++ = x_body[2];
		for(i=0;i<dv;i++)
		{
			*buf_tx++ = u_b_opt[i][0];
			*buf_tx++ = u_b_opt[i][1];
			*buf_tx++ = u_b_opt[i][2];
		}
		buf_tx = top_addr;

		sendcnt++;
		opt_cnt++;
	}

	KB_close();
	check = savedata();

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

	flag_sci = 1;

	flag_kb =1;

	flag_log = init_optimization();

	return flag_udp*flag_sci*flag_kb*flag_log;
}

void print_error(int error_code)
{
	printf("Calling function is failed.\n");
	printf("print error string : %s\n",strerror(error_code));
	printf("print error code : %d\n",error_code);
}
