#include <termios.h>
#include <unistd.h>
#include <pthread.h>
#include "prototype_declaration.h"

static struct termios old_set;
static struct termios new_set;
static int ReadChar = -1;

pthread_t ThreadCheck_kb;

extern int FlagMain;
extern int FlagUDP;

void KB_open(void)
{
	tcgetattr(0,&old_set);
	new_set = old_set;
	new_set.c_lflag &= ~ICANON;
	new_set.c_lflag &= ~ECHO;
	new_set.c_lflag &= ~ISIG;
	new_set.c_cc[VMIN] = 0;
	new_set.c_cc[VTIME] = 0;
	tcsetattr(0,TCSANOW,&old_set);
}

void KB_close(void)
{
	tcsetattr(0,TCSANOW,&old_set);
}

bool kbhit()
{
	char ch;
	int nread;
	
	if(ReadChar != -1)
	{
		return true;
	}
	
	new_set.c_cc[VMIN] = 0;
	tcsetattr(0,TCSANOW,&new_set);
	nread=read(0,&ch,1);
	new_set.c_cc[VMIN] = 1;
	tcsetattr(0,TCSANOW,&new_set);
	
	if(nread == 1)
	{
		ReadChar = ch;
		return true;
	}
	
	return false;
}

char linux_getch()
{
	char ch;
	
	if(ReadChar != 1)
	{
		ch = ReadChar;
		ReadChar = -1;
		return ch;
	}
	
	read(0,&ch,1);
	return ch;
}

void *Check_kb(void *chkb)
{
	KB_open();
	while(FlagMain)
	{
		if(kbhit())
		{
			FlagUDP = 0;
			FlagMain = 0;
			break;
		}
	}
	KB_close();
	
	pthread_exit((void *)0);
}
