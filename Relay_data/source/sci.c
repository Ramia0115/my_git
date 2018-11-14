#include <stdio.h>
#include <stdlib.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <sys/resource.h>

#include "prototype_declaration.h"
#include "setting.h"

int fd_sci;
struct termios sci_oldtio;

int FlagSCI = 0;
int sendcnt_sci = 0,sendcnt_old_sci = 0;

unsigned char *str_tx_sci;
float *buf_tx_sci;
unsigned char *str_rx_sci;
double *buf_rx_sci;

pthread_t ThreadSnd;
pthread_t ThreadRcv;

int open_serial_port(void)
{
	struct termios newtio;
	int error_num;
	int check;

	fd_sci = open(DEVICE,O_RDWR | O_NOCTTY);
	/*
	if(fd_sci == -1)
	{
		error_num = errno;
		print_error(error_num);

		FlagSCI = 0;
		return 0;
	}
	*/

	check = tcgetattr(fd_sci,&sci_oldtio);
	/*
	if(check == -1)
	{
		error_num = errno;
		print_error(error_num);

		FlagSCI = 0;
		return 0;
	}
	*/

	bzero(&newtio,sizeof(newtio));

	newtio.c_cflag = (BOUDRATE | CS8 | CLOCAL | CREAD);

	newtio.c_iflag = IGNPAR;

	newtio.c_oflag = 0;
	newtio.c_lflag = 0;

	newtio.c_cc[VTIME] = 0;
	newtio.c_cc[VMIN] = 1;

	check = tcflush(fd_sci,TCIFLUSH);
	/*
	if(check == -1)
	{
		error_num = errno;
		print_error(error_num);

		FlagSCI = 0;
		return 0;
	}
	*/

	check = tcsetattr(fd_sci,TCSANOW,&newtio);
	/*
	if(check = -1)
	{
		error_num = errno;
		print_error(error_num);

		FlagSCI = 0;
		return 0;
	}
	*/

	printf("Opening and setting device success.\n");

	check = pthread_create(&ThreadRcv,NULL,SCIrcv,NULL);
	if(check != 0)
	{
		error_num = errno;
		print_error(error_num);

		FlagSCI = 0;
		return 0;
	}
	sleep(1);
	check = pthread_create(&ThreadSnd,NULL,SCIsnd,NULL);
	if(check != 0)
	{
		error_num = errno;
		print_error(error_num);

		FlagSCI = 0;
		return 0;
	}

	return 1;
}

int close_serial_port(void)
{
	int check;
	int error_num;

	check = tcsetattr(fd_sci,TCSANOW,&sci_oldtio);
	if(check == -1)
	{
		printf("@close_serial_port.\n");
		error_num = errno;
		print_error(error_num);
		return 0;
	}

	close(fd_sci);
	return 1;
}

int get_serial_char(unsigned char c[])
{
	int readbyte;
	int error_num;

	readbyte = read(fd_sci,c,1);
	if(readbyte != 1)
	{
		error_num = errno;
		print_error(error_num);
	}

	return readbyte;
}

int put_serial_char(unsigned char c[], int rqbyte)
{
	int sendbyte;
	int error_num;

	sendbyte = write(fd_sci,c,rqbyte);
	if(sendbyte != rqbyte)
	{
		printf("@put_serial_char.\n");
		error_num = errno;
		print_error(error_num);
	}

	return sendbyte;
}

void *SCIsnd(void *SCIsnd)
{
	struct timespec t_snd;
	struct timespec t_snd_fin;

	unsigned char cmd[4] = {"MATL"};
	unsigned char size[4] = {0,4*4,0,1};
	unsigned char *data;
	unsigned char temp;

	int len;
	int i;

	str_tx_sci = (unsigned char *)malloc(sizeof(double)*256);
	buf_tx_sci = (float *)malloc(sizeof(float)*256);
	//data = (unsigned char *)malloc(sizeof(float)*(4*4));

	clock_gettime(CLOCK_MONOTONIC,&t_snd);
	t_snd.tv_nsec += SCI_CYCLE;
	while(FlagSCI)
	{
		clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&t_snd,NULL);
		clock_gettime(CLOCK_MONOTONIC,&t_snd);
		t_snd.tv_nsec += SCI_CYCLE;
		while(t_snd.tv_nsec >= NSEC_PER_SEC)
		{
			t_snd.tv_nsec -= NSEC_PER_SEC;
			t_snd.tv_sec++;
		}

		//if(sendcnt_sci > sendcnt_old_sci && sendcnt_sci%3 == 1)
		if(sendcnt_sci > sendcnt_old_sci)
		{
			len = put_serial_char(cmd,sizeof(cmd));
			len = put_serial_char(size,sizeof(size));
			//memcpy(data,(unsigned char *)buf_tx_sci,sizeof(float)*(4*4));
			//memcpy(str_tx_sci,(unsigned char *)buf_tx_sci,sizeof(float)*(4*4));
			data = (unsigned char *)buf_tx_sci;
			///*
			for(i=0;i<4*4;i++)
			{
				temp = *(data+i*4);
				*(data+i*4) = *(data+i*4+3);
				*(data+i*4+3) = temp;
				temp = *(data+i*4+1);
				*(data+i*4+1) = *(data+i*4+2);
				*(data+i*4+2) = temp;
			}
			memcpy(str_tx_sci,data,sizeof(float)*(4*4));
			//*/
			len = put_serial_char(str_tx_sci,sizeof(float)*(4*4));
			//len = put_serial_char(data,sizeof(float)*(4*4));
			sendcnt_old_sci = sendcnt_sci;
		}
	}

	pthread_exit((void *)0);
}

void *SCIrcv(void *SCIrcv)
{
	struct timespec t_rcv;

	int len;

	unsigned char *rcv_buf;
	int bytecnt;

	int i;

	str_rx_sci = (unsigned char *)malloc(sizeof(double)*256);
	buf_rx_sci = (double *)malloc(sizeof(double)*256);
	rcv_buf = (unsigned char *)malloc(sizeof(double)*256);

	while(FlagSCI)
	{
		len = get_serial_char(str_rx_sci);
		if(len < 0)
		{
			FlagSCI = 0;
			break;
		}
		else if(len > 0)
		{
			rcv_buf[i] = str_rx_sci[0];
			bytecnt += len;
			i++;
			memset(str_rx_sci,0,sizeof(double)*256);
			if(bytecnt = RECEIVE_NUM)
			{
				memcpy(buf_rx_sci,(double *)rcv_buf,sizeof(double)*RECEIVE_NUM);

				i = 0;
			}
		}
	}

	pthread_exit((void *)0);
}
