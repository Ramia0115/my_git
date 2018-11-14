#ifndef CGMRES_H

#define CGMRES_H



#include "agex.h"

#define WHEEL_NUM	4

#define t_int		0.2

#define lb			0.25

#define D_dash	0.0

void cgmres_initialization();
void init(double x[], double u[], int leg_number);
double unew(double t, double x[], double x1[], double u[], int leg_number);

#endif
