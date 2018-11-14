#ifndef OBSTACLE_DETECTION
	#define OBSTACLE_DETECTION
	#include"PRS.h"

	#define SENSOR_OFFSET		0.1					//�Z���T�̎��t���ʒu�̃I�t�Z�b�g[m]
	#define MEASURE_MIN			0.02				//LRF�̍ŏ��v������[m]�i�X�y�b�N�l�j
	#define MEASURE_MAX			5.6					//LRF�̍ő�v������[m]�i�X�y�b�N�l�j
	#define DATA_NUM			683					//LRF�̑���f�[�^��
	#define SKIP				2					//����f�[�^�̃X�L�b�v��
	#define MIN_MEASURE_OBS		0.1					//��Q�����o�����̍ŏ��l[m]
	#define MAX_MEASURE_OBS		5.6//1.5					//��Q�����o�����̍ő�l[m]
	#define THRESHOLD			0.2					//��Q���ԋ����̃X���b�V�����h[m]
	#define PI					3.14159265359
	#define PIh					1.57079632679
	#define PIq					0.78539816340
	#define MAP_HEIGHT 5.0
	#define MAP_WIDTH 7.0
	#define MAP_DELTA 0.05

	double comp_distance_sq
	(
		struct Inertial_Coord p1,
		struct Inertial_Coord p2,
		struct Inertial_Coord p3
	);

	struct Contour
	{
		int contour_data_num;
		struct Inertial_Coord contour_data[DATA_NUM];
	};

	void measure_obstacle
	(
		struct Inertial_Coord robot_pos,
		struct Inertial_Coord vertex[][VER_NUM],
		double measure_data[2][DATA_NUM]
	);

	void obstacle_detection
	(
		struct Inertial_Coord robot_pos,
		struct obstacle_data *obstacle,
		double measure_data[2][DATA_NUM]
	);
#endif
