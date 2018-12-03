#include <stdio.h>
#include <inttypes.h>
#include <tchar.h>
#include <conio.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "setting.h"
#include "prototype_declaration.h"

#pragma warning(disable:4996)

#define UDP_Port 65000
#define TRANSMISSION_NUM 3	//x,y,theta

WSAData wsaData;
SOCKET sock;
struct sockaddr_in addr;

extern struct DATA_SAVE MyDATA_SAVE[MAXLOGGINGNUM][100];
extern int rst_data;

int open_UDP()
{	
	WSAStartup(MAKEWORD(2,0),&wsaData);
	
	sock = socket(AF_INET,SOCK_DGRAM,0);
	
	addr.sin_family = AF_INET;
	addr.sin_port = htons(UDP_Port);
	addr.sin_addr.S_un.S_addr = inet_addr("192.168.109.80");
	
	return 0;
}

int send_data(int cnt)
{
	unsigned char *str;
	double  *buf;

	str = (unsigned char *)malloc(256);
	buf = (double *)malloc(256);

	memset(buf,0,sizeof(double)*TRANSMISSION_NUM);

	buf[0] = MyDATA_SAVE[rst_data][cnt].prime_x;
	buf[1] = MyDATA_SAVE[rst_data][cnt].prime_y;
	buf[2] = MyDATA_SAVE[rst_data][cnt].prime_yaw;

	memcpy(str,(unsigned char *)buf,sizeof(double)*TRANSMISSION_NUM);

	sendto(sock, (const char *)str, sizeof(double)*TRANSMISSION_NUM, 0, (struct sockaddr *)&addr, sizeof(addr));

	return 0;
}

int close_UDP()
{
	closesocket(sock);
	WSACleanup();

	return 0;
}
