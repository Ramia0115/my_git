#include <stdbool.h>

//optimization.c
int init_optimization(void);
int body_optimization(double opt_time);
void leg_update(int leg_number);
void body_update(void);
void data_logging(double rt_time);
int savedata(void);

//UDP.c
int Init_UDP(void);
void *UDPsnd(void *UDPsnd);
void *UDPrcv4Leg1(void *UDPrcv);
void *UDPrcv4Leg2(void *UDPrcv);
void *UDPrcv4Leg3(void *UDPrcv);
void *UDPrcv4Leg4(void *UDPrcv);
void *UDPrcv4Prime(void *UDPrcv);
void *UDPrcv_sq(void *UDPrcv);

//kbhit.c
void KB_open(void);
void KB_close(void);
bool kbhit();
char linux_getch();
void *Check_kb(void *chkb);

//core.c
int initialization(void);
void print_error(int error_code);
