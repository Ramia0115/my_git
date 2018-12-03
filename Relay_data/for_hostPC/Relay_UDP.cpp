#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <windows.h>
#include <conio.h>
#include <process.h>
#include <tchar.h>
//#include <winsock2.h>
//#include <ws2tcpip.h>

#pragma warning(disable:4996)

int Init_UDP(void);
unsigned __stdcall UDPsnd(void *UDPsnd);
void print_error(int error_code);
int init_logging(void);
void savedata(void);

DWORD dUDPsnd;
HANDLE hUDPsnd;

unsigned char *str_tx;
double *buf_tx;
double **logging;
int cols,rows;
int cnt;

int FlagUDP;
int sendcnt = 0,sendcnt_old = 0;

int Init_UDP(void)
{
	printf("Create thread for UDP communication.\n");
	if((hUDPsnd = (HANDLE)_beginthreadex(NULL,0,UDPsnd,NULL,0,(unsigned int *)&dUDPsnd)) == 0)
	{
		printf("failed.\n");
		return 0;
	}
	else
	{
		printf("finished.\n");
		return 1;
	}
}

unsigned __stdcall UDPsnd(void *UDPsnd)
{
	int len;
	int error_num;
	int check;

	WSAData wsaData;
	SOCKET sock_snd;
	struct sockaddr_in server;
	DWORD ipaddr;
	int multicast_TTL = 1;

	__int64 watch_a = 0;
	__int64 watch_b = 0;
	__int64 watch_c = 0,	freq_snd = 0;
	static __int64 time_buf64 = 0,	time_buf64_log = 0,	time_snd = 0;
	static int init_mstime_snd = 0;

	WSAStartup(MAKEWORD(2,0),&wsaData);

	sock_snd = socket(AF_INET,SOCK_DGRAM,0);
	if(sock_snd < 0)
	{
		printf("Calling function 'socket' is failed @*UDPsnd.\n");
		error_num = errno;
		print_error(error_num);

		_endthreadex(0);
	}

	server.sin_family = AF_INET;
	server.sin_addr.S_un.S_addr = inet_addr("239.114.214.100");
	ipaddr = inet_addr("192.168.109.20");
	server.sin_port = htons(10000);

	check = setsockopt(sock_snd,IPPROTO_IP,IP_MULTICAST_IF,(char *)&ipaddr,sizeof(ipaddr));
	if(check < 0)
	{
		printf("Calling function 'setsockopt' is failed @UDPsnd.\n");
		error_num = errno;
		print_error(error_num);

		_endthreadex(0);
	}

	check = setsockopt(sock_snd,IPPROTO_IP,IP_MULTICAST_TTL,(char *)&multicast_TTL,sizeof(multicast_TTL));
	if(check < 0)
	{
		printf("Calling function 'setsockopt' is failed @UDPsnd.\n");
		error_num = errno;
		print_error(error_num);

		_endthreadex(0);
	}

	str_tx = (unsigned char *)malloc(sizeof(double)*256);
	buf_tx = (double *)malloc(sizeof(double)*256);

	memset(str_tx,0,sizeof(str_tx));
	memset(buf_tx,0,sizeof(buf_tx));

	init_mstime_snd = (int)clock();
	QueryPerformanceCounter((LARGE_INTEGER *)&watch_a);
	while(FlagUDP)
	{
		QueryPerformanceCounter((LARGE_INTEGER *)&watch_b);
		QueryPerformanceFrequency((LARGE_INTEGER *)&freq_snd);
		//if((double)(watch_b-time_buf64)/(double)freq_snd >= 0.018
		if((double)(watch_b-time_buf64)/(double)freq_snd >= 0.002
		)
		{
			QueryPerformanceCounter((LARGE_INTEGER *)&time_buf64);
			if(sendcnt > sendcnt_old)
			{
				memcpy(str_tx,(unsigned char *)buf_tx,sizeof(double)*(4*4));
				len = sendto(sock_snd,(const char *)str_tx,sizeof(double)*(4*4),0,(struct sockaddr *)&server,sizeof(server));
				if(len != sizeof(double)*(4*4))
				{
					printf("Calling function 'sendto' is failed @UDPsnd.\n");
					error_num = errno;
					print_error(error_num);

					break;
				}

				sendcnt_old = sendcnt;
			}
		}
	}

	free(str_tx);
	free(buf_tx);

	closesocket(sock_snd);

	WSACleanup();

	_endthreadex(0);

	return 0;
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

	log_cnt = 100000;

	rows = log_cnt+1;
	cols = 4*4+1;
	logging = (double **)malloc(rows*sizeof(double));
	for(i=0;i<rows;i++)
	{
		logging[i] = (double *)malloc(cols*sizeof(double));
	}
	for(i=0;i<cols;i++)
	{
		for(j=0;j<rows;j++)
		{
			logging[j][i] = 0.0;
		}
	}
	col_cnt = 0;
	cnt = 0;

	logging[cnt][col_cnt++] = 0.0;

	for(i=0;i<4;i++)
	{
		logging[cnt][col_cnt++] = 0.0;
		logging[cnt][col_cnt++] = 0.0;
		logging[cnt][col_cnt++] = 0.0;
		logging[cnt][col_cnt++] = 0.0;
	}

	return 1;
}

void savedata(void)
{
	FILE *myfile;
	int i,j;

	myfile = fopen("logging_host.mat","w");
	for(i=0;i<cnt;i++)
	{
		for(j=0;j<cols;j++)
		{
			fprintf(myfile,"%lf\t",logging[i][j]);
		}
		fprintf(myfile,"\n");
	}
	fclose(myfile);
}
