//#ifndef AGEX

#define AGEX

/*-------------- Global Variagles -------------- */
#define tsim0	0.0
#define tf		1.0		
#define ht		0.036//0.018	
#define alpha	0.5
#define zeta	27.8//55.6	
#define hdir	0.002
#define rtol	1.e-6
#define kmax	30		
#define dv		10		

/*-------------- Dimensions -------------- */
#define DIMX   5
#define DIMUC  6
#define DIMP   5
#define DIMC   2

//çSë©èåè
void comp_constraint(double x_b[], double x_l[][DIMX], double constraint[]);

void pfunc(double t, double p1[], int i);

/*-------------- dPhi/dx -------------- */
void phix(double t, double x[], double phx1[], int i);

/*-------------- State Equation -------------- */
void xpfunc(double t, double x[], double u[], double xprime[], int dvnum, int i);

/*-------------- Costate Equation -------------- */
//void lpfunc(double t, double lmd[], double linp[], double lprime[]);
void lpfunc(double t, double lmd[], double linp[], double lprime[], int dvnum, int i);

/*-------------- Error in Optimality Condition, Hu -------------- */
void hufunc(double t, double x[], double lmd[], double u[], double hui[], int dvnum, int i);

//#endif