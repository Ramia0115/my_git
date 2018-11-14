#define MyPI	3.14159265
#define MyPIh	1.57079633
#define MyPIq	0.78539816
#define MyRad70	1.22173048
#define MyRad88	1.53589
#define Tan88	28.6363
#define MyEPS	1.0E-16
#define MyPrc	1.0E-4
#define MyBnd	1.0E-2
#define DegPRad	57.29577951

#define Max_kcnt	30
#define Min_dalpha	1.0E-5

#define GravityAcceleration	9.81

#define WHEEL_NUM 4

//#define CASE1
#define CASE2

#define X0 1.0			//����X���W
#define Y0 1.0			//����Y���W
#define theta0 0.0		//�����p���p
#define SPEED 0.5		//�ڕW���x

#define OBS_NUM 1	//��Q���̐�
#define VER_NUM 1	//��Q���̒��_�̐�

#define MAX_OBS_NUM 4	//��Q���̐��̍ő吔
#define MAX_VER_NUM 5	//��Q���̒��_�̐��̍ő吔

#include "main.h"

//���W�n�Ɋւ���\����
//�������W�n
struct Inertial_Coord
{
	double X;
	double Y;
	double theta;
};

struct PRS_variables
{
	double SelfPosition[3];
	double velocity[3];
};
struct PRSLeg_MPC2 {
	int	fr,LR;
	double th0, th1, th2;			//�֐ߊp�x(�ԗ��Œ���W�n)
	double q0, q1, q2;				//�֐ߊp�x(�������W�n)
	double q0r, q1r, q2r;
	double omega0, omega1, omega2;	//�֐ߊp���x
	double xf,yf;					//�r�ʒu
	//const dReal *xfd,*yfd;					//�r�ʒu
	double L0,L1,L2;
	double HJx,HJy;
	double r0,r1,r2,r3;				//q0,q1,q2,delta�Ɍ����֐ߊp�x�I�t�Z�b�g�l
	double Vij;
	double dcmW;
};

#define Eps_J_sq_p2 1.0E-2

#define Phi_bdy	0.99005
#define Gam_bdy	0.00995017
#define T_bdy	0.1
//iT_leg = 1 / T_leg
#define iT_bdy	10.0

//�֐ߊp�x�̈ꎟ�x��t�B���^
#define T_th	0.2

static const double Phi_th[2][2] = {{9.99987542E-001,9.95012479E-004},{-2.48753120E-002,9.90037417E-001}};
static const double Gam_th[2][1] = {{1.24584114E-005},{2.48753120E-002}};

#define T_RC	0.05
#define Phi_RC	9.80198673E-001
#define Gam_RC	1.98013267E-002

struct PRSBody {
	double fqZ;
	double fx,fy;//�ԑ̌Œ���W�n�ړ��w�����x
	double ux,uy,uphi;//�ԑ̌Œ���W�n�ړ����͑��x
	double Vxh,Vyh;//20070818:�ԑ̌Œ���W�n�ړ��w�����x�t�B���^��ԁF�����x�̐���Ɏg�p����
	double stime;//�������
	double gcz;//���{�b�g�̏d�Sz���W
	double wheelbase,tread;//����,�֋�
	double dx,dy,dqZ;//�ԑ̌Œ���W�n�̑��x�C��]���x
	double fx_last, fy_last, fqZ_last;//20080309
};

struct PRSLeg {
	int	fr,LR;
	double th0,th1,th2;
	double dth0,dth1;//20080619kn:�Q���t�B���^�ɂ��p���x
	double th0r,th1r;//�œK���v�Z�ŋ��߂��œK�p�x
	double lambda;//�֐ߊp�x�̂P���x��t�B���^�̎��萔�̋t��
	double xf,yf;//�v�Z�����r�ʒu
	double xref,yref;//�r�ʒu�ڕW�l
	double L0,L1,L2;
	double HJx,HJy;
	double r0,r1,r2,r3;//q0,q1,q2,delta�Ɍ����֐ߊp�x�I�t�Z�b�g�l
	double sr[3],cr[3],tr3;//20071021
	double WRad;//�ԗ֔��a
	double Vij;
	double dcmW;//Omega of DC motor
	double Vij0;//�v���l
	double dcmW0;//�v���l
	double xa,ya;//ankle��x,y���W
	double pk2d_last;//�O���Dp2�̒l
};


struct obstacle_data
{
	int obs_num;
	int ver_num[MAX_OBS_NUM];
	struct Inertial_Coord vertex[MAX_OBS_NUM][MAX_VER_NUM];
};

int PRS_simulation
(
	struct Inertial_Coord vertex_tmp[][VER_NUM],
	double **logging
);

struct PRSLeg_MPC {
	int	fr,LR;
	double th0,th1,th2;//�֐ߊp�x
	double xf,yf;//�r�ʒu
	double dth0,dth1,dth2;//�֐ߊp���x
	double z1,z2,z3;//�X���b�N�ϐ�
	double xref,yref;//�r�ʒu�ڕW�l
	double L0,L1,L2;
	double HJx,HJy;
	double r0,r1,r2,r3;//q0,q1,q2,delta�Ɍ����֐ߊp�x�I�t�Z�b�g�l
	double sr[3],cr[3],tr3;//20071021
	double WRad;//�ԗ֔��a
	double Vij;
	double dcmW;//Omega of DC motor
	double Vij0;//�v���l
	double dcmW0;//�v���l
	double xa,ya;//ankle��x,y���W
	double E1,E2,E3,E4,E5,E6,E7,E8,E9;//mpc�̃G���[
};

//�֐��̃v���g�^�C�v�錾
//20080728�֐����ύX
//int UpdateDiffEqs(struct PRSLeg*);
int UpdateDiffEqsLeg(struct PRSLeg*);
//int UpdateDiffEqsVxVy(struct PRSBody *);
int UpdateDiffEqsBody(struct PRSBody *);
//int InitDiffEqs(struct PRSLeg *);//20080621
int InitDiffEqsLeg(struct PRSLeg *);//20080621
//int InitDiffEqsVxVy(struct PRSBody *);
int InitDiffEqsBody(struct PRSBody *);
//
int OutputDiffEqsVxVy(struct PRSBody *, double*, double*);
int OutputDiffEqs_duxduyduphi(struct PRSBody *, double *, double *, double *);//20080619kn
int InvKinPRSQN(double *, double *, double *, struct PRSLeg, struct PRSBody);
int InvKin(struct PRSLeg *, struct PRSBody *);
int PRS_Init();
int Comp_th2_Vij_dcmW(double *, double *, double *, double *, double *, const double [], struct PRSLeg *, struct PRSBody *);
int OutputDiffEqs_ddth01(struct PRSLeg *, double *, double *);//20080619kn

int InvKinPRSMPC(double *,double *, double *, struct PRSLeg, struct PRSBody);//20141105
int PRS_Init_MPC();
int InvKin_MPC(struct PRSLeg_MPC *,struct PRSBody *,double time);
//int InvKin_MPC1(struct PRSLeg_MPC *,struct PRSBody *,double time);

int UpdateDiffEqsLeg_MPC(struct PRSLeg_MPC *w);

int CGMRES_InvKin();
int Comp_Vij_dcmW_from_th(struct PRSLeg_MPC *,struct PRSBody *,double time);