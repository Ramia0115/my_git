#ifndef MAIN

#define MAIN

#include <stdio.h>

#define KALMAN				//カルマンフィルタを利用
#define MEASUREMENT_DELAY
#ifdef MEASUREMENT_DELAY
#define LRF_DELAY 200		//LRFデータ取得の遅れ
#endif
#define LRF_ST 100			//LRFのサンプリングタイム
#define CS_ST 24			//キャスタセンサのサンプリングタイム

#define SIM_CNT 20000//10000		//シミュレーション時間[ms]
#define SIM_DT 0.001		//シミュレーションの時間刻み幅[s]
#define ROW 20000//10000			//データ取得用配列の行数(=シミュレーション時間)
#define COL 90				//データ取得用配列の列数

#endif