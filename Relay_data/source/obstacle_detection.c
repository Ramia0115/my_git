#include<stdio.h>
#include <math.h>
#include"obstacle_detection.h"

//�����t���ʐς̔{���v�Z
double compute_signed_aria(struct Inertial_Coord p1, struct Inertial_Coord p2, struct Inertial_Coord p3)
{
	return (p2.X - p1.X)*(p3.Y - p1.Y) - (p3.X - p1.X)*(p2.Y - p1.Y);
}

//�����Ɠ_�̋����̓����v�Z(�����̍��[�C�E�[�C�_�̏�)
double comp_distance_sq(struct Inertial_Coord p1, struct Inertial_Coord p2, struct Inertial_Coord p3)
{
	double a,b,c;

	a = p3.Y - p1.Y;
	b = p1.X - p3.X;
	c = p3.X*p1.Y - p1.X*p3.Y;

	return (a*p2.X + b*p2.Y + c)*(a*p2.X + b*p2.Y + c)/(a*a + b*b);
}

//�ʑ��p�`�ƂȂ钸�_��T���E�Ȃ����0��Ԃ�
int detect_vertex_2
(
	struct Inertial_Coord contour_data[],
	int left,
	int center,
	int right
)
{
	double dis_sq;
	double max_dis_sq = 0.0;
	int middle = 0;
	int i;

	//left�Ecenter�Eright�����ꂼ�ꌋ��2��������ł������ʑ��p�`���`�����钸�_��T��
	for(i=left+1;i<center-1;i++)
	{
		if(compute_signed_aria(contour_data[left],contour_data[i],contour_data[center]) > 0.0)
		{
			dis_sq = comp_distance_sq(contour_data[left],contour_data[i],contour_data[center]);
			if(dis_sq > max_dis_sq)
			{
				middle = i;
				max_dis_sq = dis_sq;
			}
		}
	}
	for(i=center+1;i<right-1;i++)
	{
		if(compute_signed_aria(contour_data[center],contour_data[i],contour_data[right]) > 0.0)
		{
			dis_sq = comp_distance_sq(contour_data[center],contour_data[i],contour_data[right]);
			if(dis_sq > max_dis_sq)
			{
				middle = i;
				max_dis_sq = dis_sq;
			}
		}
	}
	
	//�_�E�����Ԃ̋������\����������Β����Ƃ݂Ȃ�
	if(max_dis_sq < 0.0001)	middle = 0;

	return middle;
}

//�ʑ��p�`�ƂȂ钸�_��T���E�Ȃ����0��Ԃ�
int detect_vertex
(
	struct Inertial_Coord contour_data[],
	int left,
	int right
)
{
	double dis_sq;
	double max_dis_sq = 0.0;
	int middle = 0;
	int i;

	//left�Eright�����Ԑ�������ł������ʑ��p�`���`�����钸�_��T��
	for(i=left+1;i<right-1;i++)
	{
		if(compute_signed_aria(contour_data[left],contour_data[i],contour_data[right]) > 0.0)
		{
			dis_sq = comp_distance_sq(contour_data[left],contour_data[i],contour_data[right]);
			if(dis_sq > max_dis_sq)
			{
				middle = i;
				max_dis_sq = dis_sq;
			}
		}
	}

	//�_�E�����Ԃ̋������\����������Β����Ƃ݂Ȃ�
	if(max_dis_sq < 0.0001)	middle = 0;
	
	return middle;
}

//�֊s��5�p�`�ɋߎ�
void approximate_polygon
(
	struct Contour contour_data[MAX_OBS_NUM+1],	
	struct Inertial_Coord vertex[MAX_OBS_NUM][MAX_VER_NUM],
	int ver_num[],
	int *obs_num
)
{
	int p0,p1,p2,p3,p4;
	int i,j;

	for(i=0;i<*obs_num;i++)
	{
		//���_�̐���������
		ver_num[i] = 0;

		//0�_�ڂ��i�[
		p0 = 0;
		vertex[i][0] = contour_data[i].contour_data[p0];
		ver_num[i]++;

		//4�_�ڂ�T���E�Ȃ���ΏI��
		p4 = contour_data[i].contour_data_num;
		
		if(p4 == 0)
		{
			for(j=ver_num[i];j<MAX_VER_NUM;j++)
			{
				vertex[i][j].X = vertex[i][j-1].X;
				vertex[i][j].Y = vertex[i][j-1].Y;
			}
			continue;
		}
		vertex[i][4] = contour_data[i].contour_data[p4];
		ver_num[i]++;

		//2�_�ڂ�T���E�Ȃ���Β��_�������l�߂ďI��
		p2 = detect_vertex(contour_data[i].contour_data, p0, p4);
		if(p2 == 0)
		{
			vertex[i][1] = vertex[i][4];
			for(j=ver_num[i];j<MAX_VER_NUM;j++)
			{
				vertex[i][j].X = vertex[i][j-1].X;
				vertex[i][j].Y = vertex[i][j-1].Y;
			}
			continue;
		}
		vertex[i][2] = contour_data[i].contour_data[p2];
		ver_num[i]++;

		//1�_�ځC3�_�ڂ�T��
		p1 = detect_vertex(contour_data[i].contour_data, p0, p2);
		p3 = detect_vertex(contour_data[i].contour_data, p2, p4);

		if(p1 == 0)
		{
			//1�_�ځC3�_�ڋ��ɂȂ���Β��_�����l�߂ďI��
			if(p3 == 0)
			{
				vertex[i][1] = vertex[i][2];
				vertex[i][2] = vertex[i][4];
			}
			//3�_�ڂ݂̂���Ύc��̈�_��T��
			else
			{
				p1 = p2;
				vertex[i][1] = vertex[i][2];
				ver_num[i]++;
				//p2��T������
				p2 = detect_vertex_2(contour_data[i].contour_data, p1, p3, p4);
				if(p2 == 0)
				{
					vertex[i][2] = contour_data[i].contour_data[p3];
					vertex[i][3] = contour_data[i].contour_data[p4];
				}
				else if(p2 < p3)
				{
					vertex[i][2] = contour_data[i].contour_data[p2];
					vertex[i][3] = contour_data[i].contour_data[p3];
					ver_num[i]++;
				}
				else
				{
					vertex[i][2] = contour_data[i].contour_data[p3];
					vertex[i][3] = contour_data[i].contour_data[p2];
					ver_num[i]++;
				}
			}
		}
		else
		{
			//1�_�ڂ݂̂���Ύc��̈�_��T��
			if(p3 == 0)
			{
				p3 = p2;
				ver_num[i]++;
				//p2��T������
				p2 = detect_vertex_2(contour_data[i].contour_data, p0, p1, p2);
				if(p2 == 0)
				{
					vertex[i][1] = contour_data[i].contour_data[p1];
					vertex[i][3] = contour_data[i].contour_data[p4];
				}
				else if(p1 < p2)
				{
					vertex[i][3] = vertex[i][2];
					vertex[i][1] = contour_data[i].contour_data[p1];
					vertex[i][2] = contour_data[i].contour_data[p2];
					ver_num[i]++;
				}
				else
				{
					vertex[i][3] = vertex[i][2];
					vertex[i][1] = contour_data[i].contour_data[p2];
					vertex[i][2] = contour_data[i].contour_data[p1];
					ver_num[i]++;
				}
			}
			//1�_�ځC3�_�ڋ��ɂ���ΏI��
			else
			{
				vertex[i][1] = contour_data[i].contour_data[p1];
				ver_num[i]++;
				vertex[i][3] = contour_data[i].contour_data[p3];
				ver_num[i]++;
			}
		}

		for(j=ver_num[i];j<MAX_VER_NUM;j++)
		{
				vertex[i][j].X = vertex[i][j-1].X;
				vertex[i][j].Y = vertex[i][j-1].Y;
		}
	}
}

//�e��Q�����Ƃɗ֊s�f�[�^���i�[
int contour
(
	struct Inertial_Coord measure_obs[DATA_NUM],
	struct Contour contour_data[MAX_OBS_NUM+1],
	int obs_data_num
)
{
	double dx,dy,dp;
	double dis_sq;
	double min_dis_sq = 100.0;
	int obs_num = 0;
	int num = 1;
	int first=0;
	int i,j,k,l;

	//�֊s�f�[�^����������
	for(i=0;i<=MAX_OBS_NUM;i++)	contour_data[i].contour_data_num = 0;

	if(obs_data_num == 0)
	{
		//printf("check\n");
		return 0;
	}
	else
	{
		//��ڂ̏�Q����T��
		for(i=0;i<obs_data_num;i++)
		{
			j = i-1;
			if(j < 0)	j = obs_data_num - 1;
			dx = measure_obs[i].X - measure_obs[j].X;
			dy = measure_obs[i].Y - measure_obs[j].Y;
			dp = measure_obs[i].theta - measure_obs[j].theta;
			if(dp < 0.0)	dp = dp + 2.0*PI;

			//if(dx*dx + dy*dy > THRESHOLD*THRESHOLD)
			if(dx*dx + dy*dy > THRESHOLD*THRESHOLD && dp > 00.2)
			{
				first = i;
				
				break;
			}
		}
	}

	contour_data[0].contour_data[0] = measure_obs[first];
	contour_data[0].contour_data_num++;

	//printf("obs_data_num:%d\tfirst:%d\n",obs_data_num,first);
	i++;
	if(i >= obs_data_num)	i = i - obs_data_num;

	do
	{
		do
		{
			if(i == first)	break;
			j = i+1;
			if(j >= obs_data_num)	j = 0;
			//printf("%d\t%d\n",contour_data[*obs_num].contour_data_num,obs_num);

			dx = measure_obs[i].X - measure_obs[j].X;
			dy = measure_obs[i].Y - measure_obs[j].Y;

			//�f�[�^���i�[���C�֊s�f�[�^���𑝂₷
			contour_data[obs_num].contour_data[contour_data[obs_num].contour_data_num] = measure_obs[i];
			contour_data[obs_num].contour_data_num++;

			//printf("%d\t%d\t%f\t%f\n",obs_num-1,contour_data[*obs_num-1].contour_data_num,measure_obs[i].X,contour_data[obs_num-1].contour_data[contour_data[obs_num-1].contour_data_num-1].X);

			dp = measure_obs[i].theta - measure_obs[j].theta;
			if(dp < 0.0)	dp = dp + 2.0*PI;

			i++;
			if(i >= obs_data_num)	i = 0;

			//if(dx*dx + dy*dy > THRESHOLD*THRESHOLD)	break;
			if(dx*dx + dy*dy > THRESHOLD*THRESHOLD && dp > 0.02)	break;
		}
		while(1);

		obs_num++;

		//��Q���̐����K�萔�ȏ�̏ꍇ�͋K�萔�ɂȂ�܂ŏ�Q���𓝍�
		if(obs_num > MAX_OBS_NUM)
		{
			num = 1;
			min_dis_sq = 100.0;

			//�Ԋu���ŏ��̏�Q����T��
			for(k=0;k<obs_num;k++)
			{
				l = k -1;
				if(l < 0)	l = obs_num - 1;

				dx = contour_data[k].contour_data[0].X - contour_data[l].contour_data[contour_data[l].contour_data_num-1].X;
				dy = contour_data[k].contour_data[0].Y - contour_data[l].contour_data[contour_data[l].contour_data_num-1].Y;
				dis_sq = dx*dx + dy*dy;
				if(dis_sq < min_dis_sq)
				{
					min_dis_sq = dis_sq;
					num = k;
				}
				//printf("%d\t%d\t%d\t%f\t%f\n",*obs_num,k,num,dis_sq,min_dis_sq);
			}

			l = num-1;
			if(l < 0)	l = obs_num - 1;

			//��Q���𓝍��E�ۑ��̈���l�߂�
			for(k=0;k<contour_data[num].contour_data_num;k++)
			{
				contour_data[l].contour_data[contour_data[l].contour_data_num] = contour_data[num].contour_data[k];
				contour_data[l].contour_data_num++;
			}
			//for(k=num;k<*obs_num;k++)
			for(k=num;k<obs_num-1;k++)
			{
				contour_data[k] = contour_data[k+1];
			}
			contour_data[obs_num-1].contour_data_num = 0;
			obs_num--;
			//printf("%d\n",obs_num);
		}

		//printf("%d\t%d\t%d\n",obs_num,contour_data[*obs_num-1].contour_data_num,i,num);
		//printf("%d\n",obs_num);

		if(i == first)	break;
	}
	while(1);

	//printf("%d\n",obs_num);

	//�֊s�f�[�^�������炷�i0�`contour_data_num-1�ɂ���j
	for(i=0;i<obs_num;i++)	contour_data[i].contour_data_num--;

	return obs_num;
}

//�o�u���\�[�g
void bubble_sort(struct Inertial_Coord data[DATA_NUM*2], int num)
{
	int i,j;
	struct Inertial_Coord temp;

	for(i=0;i<num;i++)
	{
		for(j=num-1;j>i;j--)
		{
			if(data[j].theta > data[j-1].theta)
			{
				temp = data[j];
				data[j] = data[j-1];
				data[j-1] = temp;
			}
		}
	}
}


//�N�C�b�N�\�[�g�E�Ίp���傫�����ɕ��בւ���(�o�b�t�@�I�[�o�[�t���[����������\������)
void quick_sort(struct Inertial_Coord data[], int begin, int end)
{
	int i,j;
	double pivot;
	struct Inertial_Coord tmp;

	pivot = data[(begin + end)/2].theta;
	i = begin;
	j = end;

	while(1)
	{
		while(data[i].theta > pivot){++i;}
		while(data[j].theta < pivot){--j;}
		if(i >= j)	break;

		tmp = data[i];
		data[i] = data[j];
		data[j] = tmp;

		i++;
		j--;
	}

	if(begin < i-1)	quick_sort(data, begin, i-1);
	if(j+1 < end)	quick_sort(data, j+1, end);
}

//�Ίp�̌v�Z(���Ԃ͍������C���m�Ȓl�ł͂Ȃ�)
double alpha
(
	struct Inertial_Coord robot_pos,
	struct Inertial_Coord p
)
{
	double dx,dy,ax,ay;
	double t;

	dx = p.X - robot_pos.X;
	ax = fabs(dx);
	dy = p.Y - robot_pos.Y;
	ay = fabs(dy);

	if(ax+ay == 0.0)	t = 0.0;
	else	t = dy/(ax+ay);

	if(dx < 0.0)	t = 2.0 - t;
	else if(dy < 0.0)	t = 4.0 + t;

	return t*PIh;
}

int convert_data
(
	struct Inertial_Coord robot_pos,
	struct Inertial_Coord measure_obs[DATA_NUM],
	double measure_data[2][DATA_NUM]
)
{
	struct Inertial_Coord LRF_pos;
	double Xw,Yw;
	double length_sq;
	double phi;
	double psi[2] = {PI,0.0};
	double delta_phi = 2.0*PI/1024.0;
	int obs_data_num;
	int device;
	int i;

	obs_data_num = 0;

	for(device=0;device<2;device++)
	{
		//�Z���T�ʒu���v�Z
		LRF_pos.X = robot_pos.X + SENSOR_OFFSET*cos(robot_pos.theta + psi[device]);
		LRF_pos.Y = robot_pos.Y + SENSOR_OFFSET*sin(robot_pos.theta + psi[device]);

		for(i=0;i<DATA_NUM;i=i+SKIP)
		{
			phi = (i - 85)*delta_phi;

			Xw = LRF_pos.X + measure_data[device][i]*0.001*cos(robot_pos.theta + psi[device] - PIh + phi);
			Yw = LRF_pos.Y + measure_data[device][i]*0.001*sin(robot_pos.theta + psi[device] - PIh + phi);

			length_sq = (Xw - robot_pos.X)*(Xw - robot_pos.X) + (Yw - robot_pos.Y)*(Yw - robot_pos.Y);

			if(Xw < 0.2 || MAP_WIDTH - 0.2f < Xw || Yw < 0.2 || MAP_HEIGHT - 0.2f < Yw) continue;
			if(length_sq > MIN_MEASURE_OBS*MIN_MEASURE_OBS && length_sq < MAX_MEASURE_OBS*MAX_MEASURE_OBS)
			{
				measure_obs[obs_data_num].X = Xw;
				measure_obs[obs_data_num].Y = Yw;
				obs_data_num++;
			}
		}
	}

	//��ڂ̃f�[�^����̕Ίp���v�Z
	for(i=0;i<obs_data_num;i++)
	{
		measure_obs[i].theta = alpha(robot_pos, measure_obs[i]);
	}

	//�N�C�b�N�\�[�g
	//quick_sort(measure_obs, 0, obs_data_num-1);
	bubble_sort(measure_obs,obs_data_num);

	return obs_data_num;
}

//��_�܂ł̋����̑���
double measure_distance
(
	struct Inertial_Coord p1,
	struct Inertial_Coord p2,
	struct Inertial_Coord p3,
	struct Inertial_Coord p4
)
{
	double distance;
	double x_tmp,y_tmp;
	double ksi,eta,delta;
	double lambda,mu;

	ksi		=  (p4.Y - p3.Y)*(p4.X - p1.X) - (p4.X - p3.X)*(p4.Y - p1.Y);
	eta		= -(p2.Y - p1.Y)*(p4.X - p1.X) + (p2.X - p1.X)*(p4.Y - p1.Y);
	delta	=  (p4.Y - p3.Y)*(p2.X - p1.X) - (p4.X - p3.X)*(p2.Y - p1.Y);

	lambda = ksi/delta;
	mu = eta/delta;

	if((lambda >= 0.0 && lambda <= 1.0) && (mu >= 0.0 && mu <= 1.0))
	{
		x_tmp = p1.X + lambda*(p2.X - p1.X);
		y_tmp = p1.Y + lambda*(p2.Y - p1.Y);
		distance = sqrt((x_tmp - p1.X)*(x_tmp - p1.X) + (y_tmp - p1.Y)*(y_tmp - p1.Y));
	}
	else
	{
		distance = MEASURE_MAX;
	}

	return distance;
}

void measure_obstacle
(
	struct Inertial_Coord robot_pos,
	struct Inertial_Coord vertex[][VER_NUM],
	double measure_data[2][DATA_NUM]
)
{
	struct Inertial_Coord LRF_pos;
	struct Inertial_Coord end_pos;
	struct Inertial_Coord wall_vertex[4];
	double phi;
	double psi[2] = {PI,0.0};
	double delta_phi = 2.0*PI/1024.0;
	double distance;
	double minimum_dis;
	int device;
	int i,j,k,l;
	
	wall_vertex[0].X = 0.0;
	wall_vertex[0].Y = 0.0;
	wall_vertex[1].X = MAP_WIDTH;
	wall_vertex[1].Y = 0.0;
	wall_vertex[2].X = MAP_WIDTH;
	wall_vertex[2].Y = MAP_HEIGHT;
	wall_vertex[3].X = 0.0;
	wall_vertex[3].Y = MAP_HEIGHT;

	//LRF2��̑���
	for(device=0;device<2;device++)
	{
		//�Z���T�ʒu���v�Z
		LRF_pos.X = robot_pos.X + SENSOR_OFFSET*cos(robot_pos.theta + psi[device]);
		LRF_pos.Y = robot_pos.Y + SENSOR_OFFSET*sin(robot_pos.theta + psi[device]);

		for(i=0;i<DATA_NUM;i++)
		{
			phi = (i - 85)*delta_phi;
			end_pos.X = LRF_pos.X + MEASURE_MAX*cos(robot_pos.theta + psi[device] - PIh + phi);
			end_pos.Y = LRF_pos.Y + MEASURE_MAX*sin(robot_pos.theta + psi[device] - PIh + phi);
			//printf("%d\t%f\t%f\t%f\n",i,phi,LRF_pos.X,end_pos.X);
			
			minimum_dis = MEASURE_MAX;

			//�ǖ�
			for(j=0;j<4;j++)
			{
				k = j+1;
				if(k>3)	k = 0;

				distance = measure_distance(LRF_pos, end_pos, wall_vertex[j], wall_vertex[k]);

				if(distance < minimum_dis)	minimum_dis = distance;
			}

			//��Q��
			for(j=0;j<OBS_NUM;j++)
			{
				for(k=0;k<VER_NUM;k++)
				{
					l = k+1;
					if(l>VER_NUM-1)	l = 0;

					distance = measure_distance(LRF_pos, end_pos, vertex[j][k], vertex[j][l]);

					if(distance < minimum_dis)	minimum_dis = distance;
				}
			}

			//�Z���T�f�[�^��mm�P�ʂŎ擾
			measure_data[device][i] = minimum_dis*1000.0;
		}
	}
}

void obstacle_detection
(
	struct Inertial_Coord robot_pos,
	struct obstacle_data *obstacle,
	double measure_data[2][DATA_NUM]
)
{
	struct Inertial_Coord measure_obs[DATA_NUM];
	struct Contour contour_data[MAX_OBS_NUM+1];
	int obs_data_num = 0;

	obs_data_num = convert_data(robot_pos, measure_obs, measure_data);

	//�e��Q�����Ƃ̗֊s�f�[�^��ۑ�
	obstacle->obs_num = contour(measure_obs, contour_data, obs_data_num);

	//5�p�`�ߎ�
	approximate_polygon(contour_data, obstacle->vertex, obstacle->ver_num, &obstacle->obs_num);
}
