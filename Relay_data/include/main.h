#ifndef MAIN

#define MAIN

#include <stdio.h>

#define KALMAN				//�J���}���t�B���^�𗘗p
#define MEASUREMENT_DELAY
#ifdef MEASUREMENT_DELAY
#define LRF_DELAY 200		//LRF�f�[�^�擾�̒x��
#endif
#define LRF_ST 100			//LRF�̃T���v�����O�^�C��
#define CS_ST 24			//�L���X�^�Z���T�̃T���v�����O�^�C��

#define SIM_CNT 20000//10000		//�V�~�����[�V��������[ms]
#define SIM_DT 0.001		//�V�~�����[�V�����̎��ԍ��ݕ�[s]
#define ROW 20000//10000			//�f�[�^�擾�p�z��̍s��(=�V�~�����[�V��������)
#define COL 90				//�f�[�^�擾�p�z��̗�

#endif