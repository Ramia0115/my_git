#include <stdbool.h>

//UDP.c
int Init_UDP(void);
void *UDPsnd(void *UDPsnd);
void *UDPrcv(void *UDPrcv);
void *UDPrcv4Prime(void *UDPrcv);
void *UDPrcv_sq(void *UDPrcv);

//kbhit.c
void KB_open(void);
void KB_close(void);
bool kbhit();
char linux_getch();
void *Check_kb(void *chkb);

//sci.c
int open_serial_port(void);
int close_serial_port(void);
int get_serial_char(unsigned char c[]);
int put_serial_char(unsigned char c[], int rqbyte);
void *SCIsnd(void *SCIsnd);
void *SCIrcv(void *SCIrcv);

//core.c
int initialization(void);
void print_error(int error_code);
int init_logging(void);
void data_logging(double rt_time);
void savedata(void);
