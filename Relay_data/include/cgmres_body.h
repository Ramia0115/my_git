#ifndef CGMRES_HB

#define CGMRES_HB

#include "agex_body.h"

void cgmres_initialization_b();
void cgmres_init_b(double x[], double u[]);
double unew_b(double t, double x[], double x1[], double u[]);

#endif
