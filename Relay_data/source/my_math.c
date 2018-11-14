#include <stdio.h>
#include <math.h>

#define MyPI		3.141592653589793
#define MyPI_6		5.235987755982989E-1
#define MyPI_6_inv	1.909859317102744
#define SQRT3		1.732050807568877
#define inv_6		1.666666666666667E-1
#define inv_24		4.166666666666667E-2

double my_cos(double theta)
{
	return	cos(theta);
	/*
	int num;
	double abs_theta = fabs(theta);
	double sin_tmp, cos_tmp;
	double delta_theta;
	const double sqrt3_h = SQRT3*0.5;
	const double const_sin[13] = {0.0, 0.5, sqrt3_h, 1.0, sqrt3_h, 0.5, 0.0, -0.5, -sqrt3_h, -1.0, -sqrt3_h, -0.5, 0.0};
	const double const_cos[13] = {1.0, sqrt3_h, 0.5, 0.0, -0.5, -sqrt3_h, -1.0, -sqrt3_h, -0.5, 0.0, 0.5, sqrt3_h, 1.0};

	while(abs_theta > 2.0*MyPI){abs_theta -= 2.0*MyPI;};

	num = (int)(abs_theta*MyPI_6_inv);
	sin_tmp = const_sin[num];
	cos_tmp = const_cos[num];

	delta_theta = abs_theta - MyPI_6*num;

	return	cos_tmp + (-sin_tmp + ( -cos_tmp*0.5 + (sin_tmp*inv_6 + cos_tmp*inv_24*delta_theta)*delta_theta)*delta_theta)*delta_theta;
	*/
}

double my_sin(double theta)
{
	return	sin(theta);
	/*
	int num;
	double abs_theta = fabs(theta);
	double sin_tmp, cos_tmp;
	double delta_theta;
	const double sqrt3_h = SQRT3*0.5;
	const double const_sin[13] = {0.0, 0.5, sqrt3_h, 1.0, sqrt3_h, 0.5, 0.0, -0.5, -sqrt3_h, -1.0, -sqrt3_h, -0.5, 0.0};
	const double const_cos[13] = {1.0, sqrt3_h, 0.5, 0.0, -0.5, -sqrt3_h, -1.0, -sqrt3_h, -0.5, 0.0, 0.5, sqrt3_h, 1.0};

	while(abs_theta >= 2.0*MyPI){abs_theta -= 2.0*MyPI;};

	num = (int)(abs_theta*MyPI_6_inv);

	sin_tmp = const_sin[num];
	cos_tmp = const_cos[num];

	delta_theta = abs_theta - MyPI_6*num;

	if(theta>0.0)	return    sin_tmp + (cos_tmp + ( -sin_tmp*0.5 + (-cos_tmp*inv_6 + sin_tmp*inv_24*delta_theta)*delta_theta)*delta_theta)*delta_theta;
	else			return  -(sin_tmp + (cos_tmp + ( -sin_tmp*0.5 + (-cos_tmp*inv_6 + sin_tmp*inv_24*delta_theta)*delta_theta)*delta_theta)*delta_theta);
	*/
}

double my_square(double A)
{
	return A*A;
}

double my_sqrt(double A)
{
	return sqrt(A);
	/*
	double x;
	double delta;

	x = (A + 1.0)*0.5;

	do
	{
		delta = 0.5*(x - A/x);
		x = x - delta;
	}while(delta > 0.01);

	return x;
	*/
}
