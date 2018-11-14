#ifndef AGEXB

#define AGEXB

/*-------------- Global Variagles -------------- */

#define kmax_b	15

#define STR
//#define ROT
//#define SAW

/*-------------- Dimensions -------------- */
#define DIMXB   3
#define DIMUCB  3
#define DIMPB   3
#define DIMCB   8

/*-------------- Time-Variant Parameters --------------*/
void pfunc_b(double t, double p1[]);

/*-------------- dPhi/dx -------------- */
void phix_b(double t, double x[], double phx1[]);

/*-------------- State Equation -------------- */
void xpfunc_b(double t, double x[], double u[], double xprime[], int dvnum, int i);

/*-------------- Costate Equation -------------- */
void lpfunc_b(double t, double lmd[], double linp[], double lprime[], int divnum, int i);

/*-------------- Error in Optimality Condition, Hu -------------- */
void hufunc_b(double t, double x[], double lmd[], double u[], double hui[], int divnum);

#endif