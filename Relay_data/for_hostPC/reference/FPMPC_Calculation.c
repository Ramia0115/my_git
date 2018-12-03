#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<windows.h>
#include<process.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<time.h>

#include "FPMPC_Calculation.h"
#include "solver_FPMPC.h"
#include "FPM.h"
#include "switch.h"
#include "prototype_declaration.h"
#include "structure.h"

//center position of obstacles	{X, Y}
extern double constraint[2];
extern double goal_t;
extern double constraint_fin[2];
extern double goal_t_fin;
extern int num_iters_t, num_iters_r;
extern double ref_x;
extern double rot_target;
extern double x_final;
extern double xref[H_FPMPC+1];
extern double yref[H_FPMPC+1];
extern double tref[H_FPMPC+1];
extern double input_opt[H_FPMPC];
extern double input_opt_t[H_FPMPC+1];
extern double rot_target_fin;
extern double rot_target_mid;
extern struct Obstacle_limitation obstacle[OBS_NUM];
extern double const_d;
extern double bound[4];
double Xpredictive[ROW][12];
double Ypredictive[ROW][12];
double Tpredictive[ROW][12];
double Optimal_input[ROW][H_FPMPC];
double Optimal_input_t[ROW][H_FPMPC+1];
double delta_fpmpc[4],omega_fpmpc[4];
double state[DIMX_FPMPC]={X0_FPMPC, Y0_FPMPC, q0_FPMPC};				//Position	{X, Y, Th}
double v[DIMU]={0.0, 0.0, 0.0};					//Velocity	{vx, vy, vq}
double goal, vt;
int cnt_send_fpmpc = 0;
extern int cnt_localization;
extern int cnt_sci;
FILE *fp;

//以下実機実装時に追加
extern DWORD dFPMPC;
extern HANDLE hFPMPC;
extern struct Inertial_Coord vehicle;
//fpmpcの稼働状況を判別するフラグ
//1：動作開始	2：動作終了	3：動作準備	4：開始待ち
int flag_fpmpc;
extern int flag_operation;

//FPMPCを実行するスレッドを立ち上げる関数
int init_FPMPC_Calculation(void)
{
	if((hFPMPC=(HANDLE)_beginthreadex(NULL,0,FPMPC_Calculation,NULL,0,(unsigned int*)&dFPMPC))==0)
	//if((hFPMPC=(HANDLE)_beginthreadex(NULL,0,FPMPC_Calculation,NULL,CREATE_SUSPENDED,(unsigned int*)&dFPMPC))==0)
	{
		return 0;
	}
	else
	{
		return 1;
	}

	//CREATE_SUSPENDEDはスレッドを立ち上げたときに一時停止状態から始めるためのコマンド
	//beginthreadexは立ち上げ失敗すると0をリターンする
	//成功すればこのinit関数は1をリターンする
}

//FPMPCを実行する関数
unsigned __stdcall FPMPC_Calculation(void *lpx){
	//各変数の定義と初期化
	int i,j,t;
	double logging[ROW][COL];
	double time = 0.0;
	//double state[DIMX_FPMPC]={X0, Y0, q0};				//Position	{X, Y, Th}
    //double v[DIMU]={0.0, 0.0, 0.0};					//Velocity	{vx, vy, vq}
    double v_temp[DIMU]={0.0, 0.0, 0.0};
	double vw[DIMVw]={0.0, 0.0, 0.0, 0.0};			//Velocities of each tire	{vw1, vw2, vw3, vw4} 
//	double omega[DIMVw]={0.0, 0.0, 0.0, 0.0};		//Angular velocities of each tire
//	double obs[OBS_NUM][XY]={0.0, 0.0};
	double distance;
	static double movx,movy,mov;
	double Ts, t_old;
	static double u_trans, Th, dTh;
	static double u_rot, rot, init_MPC;
	static int flag;
	//double goal, vt;
	static double theta_old;
	int unconvergedtime_r = 0;
	int unconvergedtime_t = 0;
	struct Cartesian3 w[] =
	{	
		{	0.5*LENGTH_FPMPC,	0.5*WIDTH_FPMPC,	-0.5*HEIGHT_FPMPC	},//n=0:左前輪
		{	0.5*LENGTH_FPMPC,	-0.5*WIDTH_FPMPC,	-0.5*HEIGHT_FPMPC	},//n=1:右前輪
		{  -0.5*LENGTH_FPMPC,	0.5*WIDTH_FPMPC,	-0.5*HEIGHT_FPMPC	},//n=2:左後輪
		{  -0.5*LENGTH_FPMPC,	-0.5*WIDTH_FPMPC,	-0.5*HEIGHT_FPMPC	} //n=3:右後輪
	};

	__int64 watch_a = 0;	//timeinit
	__int64 watch_b = 0,	watch_c = 0,	freq_fpmpc = 0;
	static __int64 time_buf64 = 0,	time_buf64_log = 0,	time_LRF = 0;
	static int init_mstime_fpmpc;

	int cnt_last_localization = 0;

	//2017.09.23追加
	#ifndef FPMPC_SIM
		//ロギングの粗さを決定
		//logging_cnt:ロギングを行うたびに1，増加する．
		//logging_cycle:ロギングを何サイクル目にやるかを決定する．基本的には1サイクル18ms．
		int logging_cnt,logging_cycle;
		logging_cnt = 0;
		logging_cycle = 2;
	#endif
	u_trans = 0.0;
	Th = 0.0;

	u_rot = 0.0;
	rot = q0_FPMPC;
	theta_old = q0_FPMPC;
	
	t_old = 0.0;
	Ts = 0.018;//1.0E-2;//0.01
	flag = 0;

	vt = 0.5;	//速度は計算上では速度0.5m/sで一定

	#ifndef FPMPC_SIM
		//初期化(もしかしたらいらないかも)
		state[0] = vehicle.X;
		state[1] = vehicle.Y;
		state[2] = vehicle.Phi;
	#endif
	//printf("********Initializing...********\n");
	
	//ログ用の変数の初期化
	for(i=0;i<ROW;i++){
		for(j=0;j<COL;j++)		logging[i][j] = 0.0;
	}
	//printf("********During execution of simulation...********\n");

	Sleep(15000);
	flag_fpmpc = 3;
	Sleep(1000);

	SetupController(t_old, state[2]);	//制御器の初期化
	obstacle_detection();	//障害物情報の格納
	init_MPC = ComputeController(time, state, flag, Th);	//FPMPCの初期最適化

	//printf("flag:%d\n",flag);
	t = 0;
	flag_fpmpc = 1;
	init_mstime_fpmpc = (int)clock();
	QueryPerformanceCounter((LARGE_INTEGER *)&watch_a);	//timeinit

	//for(t=0 ; t<ROW ; t++){
	while(1){

		//time = t * dt;
		//printf("T:%lf X:%lf Y:%lf Th:%lf trans:%lf\n", time, state[0], state[1], state[2]*180/MyPI, Th*180/MyPI);
		
		//theta_old = state[2];

		QueryPerformanceCounter((LARGE_INTEGER *)&watch_b);
		QueryPerformanceFrequency((LARGE_INTEGER *)&freq_fpmpc);

		//if(time > t_old + Ts){
		if((double)(watch_b-time_buf64)/(double)freq_fpmpc >= 0.018 
			//&& cnt_localization > cnt_last_localization
		//if((double)(watch_b-time_buf64)/(double)freq_fpmpc >= 0.036
		)
		{
			QueryPerformanceCounter((LARGE_INTEGER *)&time_buf64);	//このif文に入った瞬間を記録

			//printf("T:%lf X:%lf Y:%lf Th:%lf trans:%lf\n", time, state[0], state[1], state[2]*180/MyPI, Th*180/MyPI);
			time = t * dt;

			#ifndef FPMPC_SIM
				state[0] = vehicle.X;
				state[1] = vehicle.Y;
				state[2] = vehicle.Phi;
			#endif

			theta_old = state[2];

			t_old = t_old + Ts;
			flag = 1;

			u_rot = ComputeController(time, state, flag, Th);
			if(work_fpmpc.converged == 0){
				unconvergedtime_r ++;//printf("*********** not converged ***********\n");
			}
			flag = 2;


			state[2] =  state[2] + u_rot * dt;
			
			u_trans = ComputeController(time, state, flag, Th);	
			if(work_fpmpc.converged == 0){
				unconvergedtime_t ++;
			}

			rot = state[2];

			//printf("T:%lf X:%lf Y:%lf Th:%lf trans:%lf\n", time, state[0], state[1], state[2]*180/MyPI, Th*180/MyPI);

			 dTh = u_trans;
			 Th = Th + dTh * dt;
		 
			 //得られる姿勢角 or 姿勢角速度はvx,vy,vqに変換する．
			vt = Speed(time, state[0], state[1]);	

	//		 if(distance < dc){
	////			for(i=0; i<DIMU; i++) v[i] = v[i] * distance / dc;
	//			 vt = vt * distance / dc;
	//		}
		 
			 v[0] = vt * cos(Th);
			 v[1] = vt * sin(Th);
			 v[2] = u_rot;//dTh;
			 /*
	#ifndef FPMPC_SIM
			 v[0] = v[0]*cos(state[2])-v[1]*sin(state[2]);
			 v[1] = v[0]*sin(state[2])+v[1]*cos(state[2]);
	#endif
	//*/

			 distance = sqrt((Xd - state[0])*(Xd - state[0])+(Yd - state[1])*(Yd - state[1]));	
			 ControllerKinematics(delta_fpmpc,omega_fpmpc,v[0],v[1],v[2],w);
		//	 ControllerKinematics(delta_fpmpc,omega_fpmpc,v[0]*cos(state[2])-v[1]*sin(state[2]),v[0]*sin(state[2])+v[1]*cos(state[2]),v[2],w);
		//	 ControllerKinematics(delta_fpmpc,omega_fpmpc,v[0]*cos(theta_old)+v[1]*sin(theta_old),v[0]*sin(theta_old)-v[1]*cos(theta_old),v[2],w);
		//	 ControllerKinematics(delta_fpmpc,omega_fpmpc,v[0]*cos(Th)-v[1]*sin(Th),v[0]*sin(Th)+v[1]*cos(Th),v[2],w);
			 movx = v[0]*dt;
			 movy = v[1]*dt;
			 mov += sqrt(movx*movx+movy*movy);
		 
			//Updating states of the robot
			#ifdef FPMPC_SIM
				for(i = 0; i < DIMX_FPMPC-1; i++)	state[i] = state[i] + v[i]*dt;
				goal = atan2(Yd - state[1], Xd - state[0]);
			#endif

			#ifndef FPMPC_SIM
			//ロギング
			logging_cnt++;

			if(logging_cnt == logging_cycle)	//ロギングの判別
			{
				logging[t][0] = time;
				logging[t][1] = state[0];
				logging[t][2] = state[1];
				logging[t][3] = state[2];
				logging[t][4] = v[0];
				logging[t][5] = v[1];
				logging[t][6] = v[2];
				logging[t][7] = Th;
				logging[t][8] = dTh;
				logging[t][9] = constraint[0];//th_r[1]
				logging[t][10] = constraint[1];//th_l[1]
				logging[t][11] = num_iters_t;//th_l[1]
				//logging[t][12] = goal;
				logging[t][12] = 0.0;
		//		logging[t][13] = obstacle[1].R;
		//		logging[t][14] = obstacle[1].L;
		//		logging[t][15] = obstacle[2].R;
		//		logging[t][16] = obstacle[2].L;
				logging[t][13] = bound[0];
				logging[t][14] = bound[1];
				logging[t][15] = bound[2];
				logging[t][16] = bound[3];
				logging[t][17] = goal_t;
				logging[t][18] = ref_x;
				logging[t][19] = rot_target;
				logging[t][20] = rot;
				logging[t][21] = xref[0];
				logging[t][22] = xref[1];
				logging[t][23] = num_iters_r;
				logging[t][24] = u_trans;
				logging[t][25] = x_final;
				logging[t][26] = rot_target_fin;
				logging[t][27] = rot_target_mid;
				logging[t][28] = const_d;
				logging[t][29] = constraint_fin[0];//th_r[1]
				logging[t][30] = constraint_fin[1];//th_l[1]
				logging[t][31] = goal_t_fin;//th_l[1]
				logging[t][32] = delta_fpmpc[0];
				logging[t][33] = delta_fpmpc[1];
				logging[t][34] = delta_fpmpc[2];
				logging[t][35] = delta_fpmpc[3];
				logging[t][36] = omega_fpmpc[0];
				logging[t][37] = omega_fpmpc[1];
				logging[t][38] = omega_fpmpc[2];
				logging[t][39] = omega_fpmpc[3];
				logging[t][40] = movx;
				logging[t][41] = movy;
				logging[t][42] = mov;
		//		for(i=0; i <= H_FPMPC ;i++) logging[t][i+29] = xref[i];
		//		for(i=0; i <= H_FPMPC ;i++) logging[t][i + H_FPMPC + 29] = yref[i];
				for(i=0; i <= H_FPMPC ;i++) Xpredictive[t][i] = xref[i];
				for(i=0; i <= H_FPMPC ;i++) Ypredictive[t][i] = yref[i];
				for(i=0; i <= H_FPMPC ;i++) Tpredictive[t][i] = tref[i];
				for(i=0; i < H_FPMPC ;i++) Optimal_input[t][i] = input_opt[i];
				for(i=0; i <= H_FPMPC ;i++) Optimal_input_t[t][i] = input_opt_t[i];
				t++;
				logging_cnt = 0;
			}
			#else
				logging[t][0] = time;
				logging[t][1] = state[0];
				logging[t][2] = state[1];
				logging[t][3] = state[2];
				logging[t][4] = v[0];
				logging[t][5] = v[1];
				logging[t][6] = v[2];
				logging[t][7] = Th;
				logging[t][8] = dTh;
				logging[t][9] = constraint[0];//th_r[1]
				logging[t][10] = constraint[1];//th_l[1]
				logging[t][11] = num_iters_t;//th_l[1]
				//logging[t][12] = goal;
				logging[t][12] = 0.0;
		//		logging[t][13] = obstacle[1].R;
		//		logging[t][14] = obstacle[1].L;
		//		logging[t][15] = obstacle[2].R;
		//		logging[t][16] = obstacle[2].L;
				logging[t][13] = bound[0];
				logging[t][14] = bound[1];
				logging[t][15] = bound[2];
				logging[t][16] = bound[3];
				logging[t][17] = goal_t;
				logging[t][18] = ref_x;
				logging[t][19] = rot_target;
				logging[t][20] = rot;
				logging[t][21] = xref[0];
				logging[t][22] = xref[1];
				logging[t][23] = num_iters_r;
				logging[t][24] = u_trans;
				logging[t][25] = x_final;
				logging[t][26] = rot_target_fin;
				logging[t][27] = rot_target_mid;
				logging[t][28] = const_d;
				logging[t][29] = constraint_fin[0];//th_r[1]
				logging[t][30] = constraint_fin[1];//th_l[1]
				logging[t][31] = goal_t_fin;//th_l[1]
				logging[t][32] = delta_fpmpc[0];
				logging[t][33] = delta_fpmpc[1];
				logging[t][34] = delta_fpmpc[2];
				logging[t][35] = delta_fpmpc[3];
				logging[t][36] = omega_fpmpc[0];
				logging[t][37] = omega_fpmpc[1];
				logging[t][38] = omega_fpmpc[2];
				logging[t][39] = omega_fpmpc[3];
				logging[t][40] = movx;
				logging[t][41] = movy;
				logging[t][42] = mov;
		//		for(i=0; i <= H_FPMPC ;i++) logging[t][i+29] = xref[i];
		//		for(i=0; i <= H_FPMPC ;i++) logging[t][i + H_FPMPC + 29] = yref[i];
				for(i=0; i <= H_FPMPC ;i++) Xpredictive[t][i] = xref[i];
				for(i=0; i <= H_FPMPC ;i++) Ypredictive[t][i] = yref[i];
				for(i=0; i <= H_FPMPC ;i++) Tpredictive[t][i] = tref[i];
				for(i=0; i < H_FPMPC ;i++) Optimal_input[t][i] = input_opt[i];
				for(i=0; i <= H_FPMPC ;i++) Optimal_input_t[t][i] = input_opt_t[i];
				t++;
			#endif
			//cnt_last_localization = cnt_localization;
			cnt_send_fpmpc++;
			//#ifndef CONNECT2SH
			//	cnt_sci++;
			//#endif
		}
		/*
		logging[t][0] = time;
		logging[t][1] = state[0];
		logging[t][2] = state[1];
		logging[t][3] = state[2];
		logging[t][4] = v[0];
		logging[t][5] = v[1];
		logging[t][6] = v[2];
		logging[t][7] = Th;
		logging[t][8] = dTh;
		logging[t][9] = constraint[0];//th_r[1]
		logging[t][10] = constraint[1];//th_l[1]
		logging[t][11] = num_iters_t;//th_l[1]
		logging[t][12] = goal;
//		logging[t][13] = obstacle[1].R;
//		logging[t][14] = obstacle[1].L;
//		logging[t][15] = obstacle[2].R;
//		logging[t][16] = obstacle[2].L;
		logging[t][13] = bound[0];
		logging[t][14] = bound[1];
		logging[t][15] = bound[2];
		logging[t][16] = bound[3];
		logging[t][17] = goal_t;
		logging[t][18] = ref_x;
		logging[t][19] = rot_target;
		logging[t][20] = rot;
		logging[t][21] = xref[0];
		logging[t][22] = xref[1];
		logging[t][23] = num_iters_r;
		logging[t][24] = u_trans;
		logging[t][25] = x_final;
		logging[t][26] = rot_target_fin;
		logging[t][27] = rot_target_mid;
		logging[t][28] = const_d;
		logging[t][29] = constraint_fin[0];//th_r[1]
		logging[t][30] = constraint_fin[1];//th_l[1]
		logging[t][31] = goal_t_fin;//th_l[1]
		logging[t][32] = delta[0];
		logging[t][33] = delta[1];
		logging[t][34] = delta[2];
		logging[t][35] = delta[3];
		logging[t][36] = omega[0];
		logging[t][37] = omega[1];
		logging[t][38] = omega[2];
		logging[t][39] = omega[3];
		logging[t][40] = movx;
		logging[t][41] = movy;
		logging[t][42] = mov;
//		for(i=0; i <= H ;i++) logging[t][i+29] = xref[i];
//		for(i=0; i <= H ;i++) logging[t][i + H + 29] = yref[i];
		for(i=0; i <= H ;i++) Xpredictive[t][i] = xref[i];
		for(i=0; i <= H ;i++) Ypredictive[t][i] = yref[i];
		for(i=0; i <= H ;i++) Tpredictive[t][i] = tref[i];
		for(i=0; i < H ;i++) Optimal_input[t][i] = input_opt[i];
		for(i=0; i <= H ;i++) Optimal_input_t[t][i] = input_opt_t[i];
		*/
		//停止条件
		//計算が発散したら(正確には姿勢角が非数だったら)強制終了
		if(_isnan(state[2])!=0)
		{
			break;
		}
		//何かキーボードが押されたら終了
		if(_kbhit())
		{
			break;
		}
		if(t>=ROW)
		{
			break;
		}
		//t++;
	
	}

	flag_fpmpc = 2;

	printf("********Completed********\n");
	printf("unconvergedtime. in rot = %d. in tra = %d.\n",unconvergedtime_r,unconvergedtime_t);
	fopen_s(&fp,"sim.mat","w");


	fprintf(fp, "# %d %d \n", COL, t);
	for(i=0;i<t;i++)
	{
		for(j=0;j<COL;j++)
		{
			fprintf(fp, "%lf\t", logging[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

		fopen_s(&fp,"xpre.mat","w");
	fprintf(fp, "# %d %d \n", 12, t);
	for(i=0;i<t;i++)
	{
		for(j=0;j<12;j++)
		{
			fprintf(fp, "%lf\t", Xpredictive[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fopen_s(&fp,"ypre.mat","w");
	fprintf(fp, "# %d %d \n", 12, t);
	for(i=0;i<t;i++)
	{
		for(j=0;j<12;j++)
		{
			fprintf(fp, "%lf\t", Ypredictive[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fopen_s(&fp,"tpre.mat","w");
	fprintf(fp, "# %d %d \n", 12, t);
	for(i=0;i<t;i++)
	{
		for(j=0;j<12;j++)
		{
			fprintf(fp, "%lf\t", Tpredictive[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fopen_s(&fp,"ur.mat","w");
	fprintf(fp, "# %d %d \n", 11, t);
	for(i=0;i<t;i++)
	{
		for(j=0;j<12;j++)
		{
			fprintf(fp, "%lf\t", Optimal_input[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fopen_s(&fp,"ut.mat","w");
	fprintf(fp, "# %d %d \n", 12, t);
	for(i=0;i<t;i++)
	{
		for(j=0;j<12;j++)
		{
			fprintf(fp, "%lf\t", Optimal_input_t[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	flag_fpmpc = 0;
	flag_operation = 0;

	_endthreadex(0); //スレッドを終了

	return 0;
}





