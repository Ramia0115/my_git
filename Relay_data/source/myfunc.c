#include <math.h>
#include<stdio.h>

#include "cgmres.h"
#include "rhfuncu.h"
#include "myfunc.h"
#include "my_math.h"
#include "PRS.h"

/*-------------- My variables -------------- */

//struct Inertial_Coord obs[OBS_NUM][VER_NUM] = {{{40.7, -0.2, 0.0},{40.7, -0.5, 0.0},{50.0, -0.5, 0.0},{50.0, -0.2, 0.0}}};
struct Inertial_Coord obs[OBS_NUM][VER_NUM] = {{{50, 0.0, 0.0}}};

extern double x_l_opt[WHEEL_NUM][dv + 1][DIMX];

/*-------------- My Functions -------------- */
double signed_aria(struct Inertial_Coord p1, struct Inertial_Coord p2, struct Inertial_Coord p3)
{
	return (p2.X - p1.X)*(p3.Y - p1.Y) - (p3.X - p1.X)*(p2.Y - p1.Y);
}

double comp_D
(
	struct Inertial_Coord point,
	struct Inertial_Coord vertex[4],
	int ver_num,
	int left,
	int right
)
{
	double D;
	double p_l,p_r,r_l;
	double p_l_X,p_r_X,r_l_X;
	double p_l_Y,p_r_Y,r_l_Y;
	int j,k;

	p_l_X = point.X - vertex[left].X;
	p_l_Y = point.Y - vertex[left].Y;

	p_r_X = point.X - vertex[right].X;
	p_r_Y = point.Y - vertex[right].Y;

	p_l = sqrt(p_l_X*p_l_X + p_l_Y*p_l_Y);
	p_r = sqrt(p_r_X*p_r_X + p_r_Y*p_r_Y);

	j = left;
	r_l = 0.0;
	do
	{
		k = j+1;
		if(k > ver_num-1)	k = 0;

		r_l_X = vertex[j].X - vertex[k].X;
		r_l_Y = vertex[j].Y - vertex[k].Y;
		r_l = r_l + sqrt(r_l_X*r_l_X + r_l_Y*r_l_Y);

		if(k == right)	break;

		j++;
		if(j > ver_num-1)	j = 0;
	}while(1);

	D = p_r + p_l - r_l;

	return D;
}

int comp_verge_point
(
	struct Inertial_Coord point,
	struct Inertial_Coord vertex[4],
	int ver_num,
	int *left,
	int *right
)
{
	int j,k,l;

	//左右端の頂点を探索
	*left = 0;
	*right = 0;

	for(k=0;k<ver_num;k++)
	{
		j = k-1;
		if(j < 0)	j = ver_num-1;
		l = k+1;
		if(l > ver_num-1)	l = 0;

		if(signed_aria(vertex[j], vertex[k], point) >= 0.0 && signed_aria(vertex[k], vertex[l], point) < 0.0)	*left = k;
		if(signed_aria(vertex[j], vertex[k], point) <  0.0 && signed_aria(vertex[k], vertex[l], point) >= 0.0)	*right = k;
	}

	return 0;
}

int comp_minimum_D
(
	struct Inertial_Coord wheel_point[WHEEL_NUM],
	struct Inertial_Coord ver1,		//vertexの小さいほう
	struct Inertial_Coord ver2,		//vertexの大きいほう
	struct Inertial_Coord *nearest_point,
	double *Dmin,
	int wheel_number,
	int i,
	int j
)
{
	struct Inertial_Coord near_tmp;
	double Dtmp;
	double a,b;
	double denominator,numerator;
	int i_m,i_p;
	int j_m,j_p;

	if(ver2.X - ver1.X != 0.0)
	{
		a = (ver2.Y - ver1.Y)
			/(ver2.X - ver1.X);
		b = (- ver1.X*ver2.Y + ver2.X*ver1.Y)
			/(ver2.X - ver1.X);

		denominator = a*(wheel_point[i].X - wheel_point[j].X) - wheel_point[i].Y + wheel_point[j].Y;
		if(denominator != 0.0)
		{
			numerator = -b*(wheel_point[i].X + wheel_point[j].X) - wheel_point[j].X*wheel_point[i].Y + wheel_point[i].X*wheel_point[j].Y;
			near_tmp.X = numerator/denominator;
			near_tmp.Y = a*near_tmp.X + b;

			if((near_tmp.X < ver1.X && near_tmp.X < ver2.X)
				||(near_tmp.X > ver1.X && near_tmp.X > ver2.X))
			{
				denominator = (1.0 + a*a)*(2.0*b + a*wheel_point[i].X + a*wheel_point[j].X - wheel_point[i].Y - wheel_point[j].Y);
				numerator = -2.0*a*b*b + b*wheel_point[i].X - a*a*b*wheel_point[i].X + b*wheel_point[j].X - a*a*b*wheel_point[j].X + 2.0*a*wheel_point[i].X*wheel_point[j].X + 2.0*a*b*wheel_point[i].Y - wheel_point[j].X*wheel_point[i].Y + a*a*wheel_point[j].X*wheel_point[i].Y + 2.0*a*b*wheel_point[j].Y - wheel_point[i].X*wheel_point[j].Y + a*a*wheel_point[i].X*wheel_point[j].Y - 2.0*a*wheel_point[i].Y*wheel_point[j].Y;
				near_tmp.X = numerator/denominator;
				near_tmp.Y = a*near_tmp.X + b;

				if((near_tmp.X < ver1.X && near_tmp.X < ver2.X)
					||(near_tmp.X > ver1.X && near_tmp.X > ver2.X))
				{
					return 0;
				}
			}
		}
		else
		{
			denominator = (1.0 + a*a)*(2.0*b + a*wheel_point[i].X + a*wheel_point[j].X - wheel_point[i].Y - wheel_point[j].Y);
			numerator = -2.0*a*b*b + b*wheel_point[i].X - a*a*b*wheel_point[i].X + b*wheel_point[j].X - a*a*b*wheel_point[j].X + 2.0*a*wheel_point[i].X*wheel_point[j].X + 2.0*a*b*wheel_point[i].Y - wheel_point[j].X*wheel_point[i].Y + a*a*wheel_point[j].X*wheel_point[i].Y + 2.0*a*b*wheel_point[j].Y - wheel_point[i].X*wheel_point[j].Y + a*a*wheel_point[i].X*wheel_point[j].Y - 2.0*a*wheel_point[i].Y*wheel_point[j].Y;
			near_tmp.X = numerator/denominator;
			near_tmp.Y = a*near_tmp.X + b;

			if((near_tmp.X < ver1.X && near_tmp.X < ver2.X)
				||(near_tmp.X > ver1.X && near_tmp.X > ver2.X))
			{
				return 0;
			}
		}
	}
	else if(ver2.Y - ver1.Y != 0.0)
	{
		b = (- ver2.X*ver1.Y + ver1.X*ver2.Y)
			/(ver2.Y - ver1.Y);

		denominator = 2.0*b - wheel_point[i].X -wheel_point[j].X;
		numerator = b*(wheel_point[i].Y + wheel_point[j].Y) - wheel_point[j].X*wheel_point[i].Y - wheel_point[i].X*wheel_point[j].Y;

		near_tmp.X = b;
		near_tmp.Y = numerator/denominator;

		if((near_tmp.Y < ver1.Y && near_tmp.Y < ver2.Y)
		|| (near_tmp.Y > ver1.Y && near_tmp.Y > ver2.Y))
		{
			denominator = wheel_point[i].X -wheel_point[j].X;
			numerator = b*(wheel_point[i].Y - wheel_point[j].Y) - wheel_point[j].X*wheel_point[i].Y + wheel_point[i].X*wheel_point[j].Y;
			near_tmp.Y = numerator/denominator;

			if((near_tmp.Y < ver1.Y && near_tmp.Y < ver2.Y)
			|| (near_tmp.Y > ver1.Y && near_tmp.Y > ver2.Y))
			{
				return 0;
			}
		}
	}
	else
	{
		return 0;
	}

	//i,jがそれぞれ左端，右端であるかを判別
	i_m = i - 1;
	if(i_m < 0)	i_m = wheel_number - 1;
	i_p = i + 1;
	if(i_p > wheel_number - 1)	i_p = 0;
	j_m = j - 1;
	if(j_m < 0)	j_m = wheel_number - 1;
	j_p = j + 1;
	if(j_p > wheel_number - 1)	j_p = 0;

	if((signed_aria(wheel_point[i_m], wheel_point[i], near_tmp) >= 0.0 && signed_aria(wheel_point[i], wheel_point[i_p], near_tmp) < 0.0)
	&& (signed_aria(wheel_point[j_m], wheel_point[j], near_tmp) <  0.0 && signed_aria(wheel_point[j], wheel_point[j_p], near_tmp) >= 0.0))
	{
		Dtmp = comp_D(near_tmp, wheel_point, wheel_number, i, j);

		if(Dtmp < *Dmin)
		{
			*Dmin = Dtmp;
			nearest_point->X = near_tmp.X;
			nearest_point->Y = near_tmp.Y;
		}
	}

	return 1;
}

//障害物外郭上の最近傍点を計算
int comp_nearest_point
(
	struct Inertial_Coord wheel_point[WHEEL_NUM],
	struct Inertial_Coord *nearest_point,
	struct Inertial_Coord vertex[VER_NUM],
	int vertex_number
)
{
	double D;
	double D_min;
	int wheel_number = WHEEL_NUM;
	int left,right;
	int nearest_ver_num;
	int previous_num;
	int next_num;
	int i,j;

	int left_edge[2], right_edge[2];

	D_min = 100.0;
	nearest_ver_num = 0;

	//障害物の頂点のうちDが最小のものを探索
	for(i=0;i<vertex_number;i++)
	{
		comp_verge_point(vertex[i], wheel_point, wheel_number, &left, &right);

		D = comp_D(vertex[i], wheel_point, wheel_number, left, right);

		if(D < D_min)
		{
			D_min = D;
			nearest_ver_num = i;
		}
	}

	nearest_point->X = vertex[nearest_ver_num].X;
	nearest_point->Y = vertex[nearest_ver_num].Y;
	
	//障害物外角上の最小のDを探索(i:左側、j:右側)
	//最近傍点の反時計回り側の直線
	previous_num = nearest_ver_num - 1;
	if(previous_num < 0)	previous_num = vertex_number - 1;

	comp_verge_point(vertex[nearest_ver_num], wheel_point, wheel_number, &left_edge[0], &right_edge[0]);
	comp_verge_point(vertex[previous_num], wheel_point, wheel_number, &left_edge[1], &right_edge[1]);

	i = left_edge[0];
	do
	{
		j = right_edge[0];
		do
		{
			//最小のDを計算
			comp_minimum_D(wheel_point, vertex[previous_num], vertex[nearest_ver_num],	nearest_point, &D_min, wheel_number, i, j);

			if(j == right_edge[1])	break;
			
			j++;
			if(j >= wheel_number)	j = 0;
		}while(1);
		
		if(i == left_edge[1])	break;
		
		i++;
		if(i >= wheel_number)	i = 0;
	}while(1);
	
	//最近傍点の時計回り側の直線
	next_num = nearest_ver_num + 1;
	if(next_num > vertex_number - 1)	next_num = 0;

	comp_verge_point(vertex[next_num], wheel_point, wheel_number, &left_edge[0], &right_edge[0]);
	comp_verge_point(vertex[nearest_ver_num], wheel_point, wheel_number, &left_edge[1], &right_edge[1]);
	
	i = left_edge[0];
	do
	{
		j = right_edge[0];
		do
		{
			//最小のDを計算
			comp_minimum_D(wheel_point, vertex[nearest_ver_num], vertex[next_num],	nearest_point, &D_min, wheel_number, i, j);
			
			if(j == right_edge[1])	break;
			
			j++;
			if(j >= wheel_number)	j = 0;
		}while(1);
		
		if(i == left_edge[1])	break;
		
		i++;
		if(i >= wheel_number)	i = 0;
	}while(1);

	/*
	//障害物外角上の最小のDを探索(i:左側、j:右側)
	//最近傍点の反時計回り側の直線
	previous_num = nearest_ver_num - 1;
	if(previous_num < 0)	previous_num = vertex_number - 1;

	for(i=0;i<wheel_number;i++)
	{
		for(j=0;j<wheel_number;j++)
		{
			if(i == j)	continue;
			//最小のDを計算
			if(comp_minimum_D(wheel_point, vertex[previous_num], vertex[nearest_ver_num],	nearest_point, &D_min, wheel_number, i, j) == 0)	continue;
		}
	}
	//最近傍点の時計回り側の直線
	next_num = nearest_ver_num + 1;
	if(next_num > vertex_number - 1)	next_num = 0;

	for(i=0;i<wheel_number;i++)
	{
		for(j=0;j<wheel_number;j++)
		{
			if(i == j)	continue;
			//最小のDを計算
			if(comp_minimum_D(wheel_point, vertex[nearest_ver_num], vertex[next_num],	nearest_point, &D_min, wheel_number, i, j) == 0)	continue;
		}
	}
	*/

	return 0;
}

//障害物に対するポテンシャル
void nabla_oval_potential
(
	struct Inertial_Coord wheel_point[WHEEL_NUM],
	struct Inertial_Coord vertex[VER_NUM],
	double dpot_dx[][DIMX],
	int ver_num
)
{
	
	struct Inertial_Coord obstacle_point;
	double D;
	double dDdx[WHEEL_NUM][DIMX];
	double dDRdx[WHEEL_NUM][DIMX];
	double dDLdx[WHEEL_NUM][DIMX];
	double ddfdx[WHEEL_NUM][DIMX];
	double invD;
	double p_l,p_r,r_l;
	double p_l_X,p_r_X,r_l_X;
	double p_l_Y,p_r_Y,r_l_Y;
	double tmp;
	int wheel_number = WHEEL_NUM;
	int right,left;
	int j,k;
	
	comp_nearest_point(wheel_point, &obstacle_point, vertex, ver_num);
	
	for(j=0;j<WHEEL_NUM;j++)
	{
		for(k=0;k<DIMX;k++)
		{
			dDRdx[j][k] = 0.0;
			dDLdx[j][k] = 0.0;
			ddfdx[j][k] = 0.0;
		}
	}
	
	comp_verge_point(obstacle_point, wheel_point, wheel_number, &left, &right);
	
	p_l_X = obstacle_point.X - wheel_point[left].X;
	p_l_Y = obstacle_point.Y - wheel_point[left].Y;

	p_r_X = obstacle_point.X - wheel_point[right].X;
	p_r_Y = obstacle_point.Y - wheel_point[right].Y;

	p_l = sqrt(p_l_X*p_l_X + p_l_Y*p_l_Y);
	p_r = sqrt(p_r_X*p_r_X + p_r_Y*p_r_Y);

	dDLdx[left][0] = -p_l_X/p_l;
	dDLdx[left][1] = -p_l_Y/p_l;
	dDRdx[right][0] = -p_r_X/p_r;
	dDRdx[right][1] = -p_r_Y/p_r;

	j = left;
	r_l = 0.0;
	
	do
	{
		k = j+1;
		if(k > wheel_number-1)	k = 0;

		r_l_X = wheel_point[j].X - wheel_point[k].X;
		r_l_Y = wheel_point[j].Y - wheel_point[k].Y;
		tmp = sqrt(r_l_X*r_l_X + r_l_Y*r_l_Y);
		r_l = r_l + tmp;

		ddfdx[j][0] +=  r_l_X/tmp;
		ddfdx[j][1] +=  r_l_Y/tmp;
		ddfdx[k][0] += -r_l_X/tmp;
		ddfdx[k][1] += -r_l_Y/tmp;

		if(k == right)	break;

		j++;
		if(j > wheel_number-1)	j = 0;
	}while(1);
	
	D = p_r + p_l - r_l - D_dash;

	invD = 1.0/D;
	tmp = -2.0*invD*invD*invD;
	
	for(j=0;j<WHEEL_NUM;j++)
	{
		for(k=0;k<DIMX;k++)
		{
			dDdx[j][k] = dDLdx[j][k] + dDRdx[j][k] - ddfdx[j][k];
			dpot_dx[j][k] = tmp*dDdx[j][k];
		}
	}
	
}

void comp_potential_for_autogenu
(
	double px[DIMX],
	double pdpdx[DIMX],
	int legnum,
	int dvnum
)
{
	struct Inertial_Coord wheel_point[WHEEL_NUM];
	double dpdx_tmp[WHEEL_NUM][DIMX];
	int i,j;

	for(i=0;i<WHEEL_NUM;i++)
	{
		if(i == legnum)
		{
			wheel_point[i].X = px[0];
			wheel_point[i].Y = px[1];
			wheel_point[i].theta = 0.0;
		}
		else
		{
			wheel_point[i].X = x_l_opt[i][dvnum][0];
			wheel_point[i].Y = x_l_opt[i][dvnum][1];
			wheel_point[i].theta =  0.0;
		}
	}

	//各車輪の勾配を初期化
	for(i=0;i<WHEEL_NUM;i++)
	{
		dpdx_tmp[i][0] = 0.0;
		dpdx_tmp[i][1] = 0.0;
	}

	//最適化対象の車輪の勾配を初期化
	for (j = 0; j<DIMX; j++)	pdpdx[j] = 0.0;

	//各障害物に対して勾配を計算
	for(i=0;i<OBS_NUM;i++)
	{
		nabla_oval_potential(wheel_point, obs[i], dpdx_tmp, VER_NUM);
		
		//for(j=0;j<DIMX;j++)	pdpdx[j] += dpdx_tmp[leg_opt][j];

		pdpdx[0] += dpdx_tmp[legnum][0];
		pdpdx[1] += dpdx_tmp[legnum][1];
	}
}