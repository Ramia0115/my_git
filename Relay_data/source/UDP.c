#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <sys/resource.h>
#include <unistd.h>

#include "prototype_declaration.h"
#include "setting.h"

int multicastPort[5] = {10000,10001,10002,10003,10004};
int FlagUDP;
int StartFlag;
int FlagRcv;
int sendcnt = 0,sendcnt_old = 0;
unsigned char *str_tx;
double *buf_tx;
unsigned char *str_rx;
double *buf_rx;
unsigned char *str_rx_prime;
double *buf_rx_prime;
double stamp;
double time_stamp;

pthread_t ThreadUDPsnd;
pthread_t ThreadUDPrcv;
pthread_t ThreadUDPrcv_sq;

int Init_UDP(void)
{
	int error_num;
	int build_th = 0;
	int build_rcv = 0;
	
	printf("Create thread for udp communication.\n");
	/*
	build_th = pthread_create(&ThreadUDPsnd,NULL,UDPsnd,NULL);
	if(build_th != 0)
	{
		error_num = errno;
		print_error(error_num);
		FlagUDP = 0;
		
		return 0;
	}
	sleep(1);
	*/
	build_rcv = pthread_create(&ThreadUDPrcv,NULL,UDPrcv,NULL);
	if(build_rcv != 0)
	{
		error_num = errno;
		print_error(error_num);
		FlagUDP = 0;
	}

	if(build_th+build_rcv == 0)
	{
		printf("finish.\n");
		return 1;
	}
	else
	{
		printf("failed.\n");
		return 0;
	}
}

void *UDPsnd(void *UDPsnd)
{
	int len;
	int error_num;
	int sock_snd;
	int check;
	struct sockaddr_in server;
	unsigned char multicastTTL = 1;
	in_addr_t ipaddr;

	struct timespec t_snd;
	struct timespec t_m;

	sock_snd = socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP);
	if(sock_snd < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}
	check = setsockopt(sock_snd,IPPROTO_IP,IP_MULTICAST_TTL,(void *)&multicastTTL,sizeof(multicastTTL));
	if(check < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}
	
	server.sin_family = AF_INET;
	server.sin_addr.s_addr = inet_addr("239.114.214.104");
	ipaddr = inet_addr("192.168.1.30");
	server.sin_port = htons(multicastPort[4]);
	
	str_tx = (unsigned char *)malloc(sizeof(double)*256);
	buf_tx = (double *)malloc(sizeof(double)*256);
	
	memset(str_tx,0,sizeof(str_tx));
	memset(buf_tx,0,sizeof(buf_tx));
	
	clock_gettime(CLOCK_MONOTONIC,&t_snd);
	t_snd.tv_nsec += UDP_CYCLE;
	while(FlagUDP)
	{
		clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&t_snd,NULL);
		clock_gettime(CLOCK_MONOTONIC,&t_snd);
		t_snd.tv_nsec += UDP_CYCLE;
		while(t_snd.tv_nsec >= NSEC_PER_SEC)
		{
			t_snd.tv_nsec -= NSEC_PER_SEC;
			t_snd.tv_sec++;
		}

		if(sendcnt > sendcnt_old)
		{
			clock_gettime(CLOCK_MONOTONIC,&t_m);
			stamp = (double)t_m.tv_sec*NSEC_PER_SEC+(double)t_m.tv_nsec;
			memcpy(str_tx,(unsigned char *)buf_tx,sizeof(double)*(10*3+3));
			len = sendto(sock_snd,str_tx,sizeof(double)*(10*3+3),0,(struct sockaddr *)&server,sizeof(server));
			if(len != sizeof(double)*(10*3+3))
			{
				error_num = errno;
				print_error(error_num);

				break;
			}

			sendcnt_old = sendcnt;
		}
	}

	free(str_tx);
	free(buf_tx);

	close(sock_snd);
	pthread_exit((void *)0);
}

void *UDPrcv(void *UDPrcv)
{
	int sock_rcv;
	struct sockaddr_in client;
	struct ip_mreq multicastRequest;
	struct timespec t_m;

	int len;
	int check;
	int error_num;

	sock_rcv = socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP);
	if(sock_rcv < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}

	client.sin_family = AF_INET;
	client.sin_addr.s_addr = INADDR_ANY;
	client.sin_port = htons(multicastPort[0]);
	check = bind(sock_rcv,(struct sockaddr *)&client,sizeof(client));
	if(check < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}

	multicastRequest.imr_multiaddr.s_addr = inet_addr("239.114.214.100");
	multicastRequest.imr_interface.s_addr = inet_addr("192.168.1.30");

	check = setsockopt(sock_rcv,IPPROTO_IP,IP_ADD_MEMBERSHIP,(void *)&multicastRequest,sizeof(multicastRequest));
	if(check < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}

	str_rx = (unsigned char *)malloc(sizeof(double)*256);
	buf_rx = (double *)malloc(sizeof(double)*256);

	while(FlagUDP)
	{
		memset(str_rx,0,sizeof(str_rx));
		len = recvfrom(sock_rcv,str_rx,sizeof(double)*(4*4),0,NULL,0);
		if(len == -1)
		{
			error_num = errno;
			print_error(error_num);

			break;
		}
		else if(len > 0)
		{
			clock_gettime(CLOCK_MONOTONIC,&t_m);
			time_stamp = (double)t_m.tv_sec*NSEC_PER_SEC+(double)t_m.tv_nsec;
			memcpy(buf_rx,(double *)str_rx,sizeof(double)*(4*4));
			FlagRcv = 1;
		}
	}

	free(str_rx);
	free(buf_rx);

	close(sock_rcv);

	pthread_exit((void *)0);
}

void *UDPrcv_sq(void *UDPrcv)
{
	int len;
	int error_num;
	int check;
	int sock_rcv;
	struct sockaddr_in client;
	struct ip_mreq multicastRequest;
	unsigned char *str_rx_sq;
	unsigned char *buf_rx_sq;

	sock_rcv = socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP);
	if(sock_rcv < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}

	client.sin_family = AF_INET;
	client.sin_addr.s_addr = INADDR_ANY;
	client.sin_port = htons(20000);
	check = bind(sock_rcv,(struct sockaddr *)&client,sizeof(client));
	if(check < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}

	multicastRequest.imr_multiaddr.s_addr = inet_addr("239.114.214.86");
	multicastRequest.imr_interface.s_addr = inet_addr("192.168.120.200");

	check = setsockopt(sock_rcv,IPPROTO_IP,IP_ADD_MEMBERSHIP,(void *)&multicastRequest,sizeof(multicastRequest));
	if(check < 0)
	{
		error_num = errno;
		print_error(error_num);

		pthread_exit((void *)0);
	}

	str_rx_sq = (unsigned char *)malloc(3);
	buf_rx_sq = (unsigned char *)malloc(3);
	memset(str_rx_sq,0,sizeof(str_rx_sq));

	while(1)
	{
		len = recvfrom(sock_rcv,str_rx_sq,sizeof(unsigned char),0,NULL,0);
		if(len == -1)
		{
			error_num = errno;
			print_error(error_num);

			break;
		}
		else if(len != 0)
		{
			memcpy(buf_rx_sq,(double *)str_rx_sq,sizeof(unsigned char)*2);
			if(buf_rx_sq[0] == 0x01)
			{
				StartFlag =1;
			}
			else if(buf_rx_sq[0] == 0x00)
			{
				StartFlag = 0;
				break;
			}
		}
	}

	pthread_exit((void *)0);
}
