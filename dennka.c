/*
 *   diff_method --
 *     差分法による無限長同軸円筒電極の電界計算プログラム
 *
 *                                         2020/7
 *                                           小野　亮
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// パラメータの定義
#define N	100	// r方向の格子数（格子点はN + 1個）。
#define M	2	// z方向の格子数（格子点はM + 1個）（NとMは2以上にすること)
#define RMAX	0.1	// r座標の最大値[m]。ZMAX = RMAX/N * M
#define PHI0	Vmax	// 全格子点の電位の初期値 (0, Vmax, 0.3*Vmaxなど)
#define ERR_LIM	1e-8	// 収束判定誤差

/* ポインタの二次元配列を見やすくするための置き換え。
   例えば、i = 0~N、j = 0~M */
#define FLAG(i, j)  (*(flag + (i)*(M+1) + (j)))  // 格子点(i, j)の境界条件の有無
#define PHI(i, j)   (*(phi + (i)*(M+1) + (j)))   // (i, j)の電位
#define ES(i, j)    (*(Es + (i)*(M+1) + (j)))    // (i, j)の電界

double *phi, *Es;
double a, b, c, lin, lout, t, r, z, h, Er, Ez, Emax, Vmax;
double phi_max, error_max, error, phi_prev, error_limit;
int    *flag;
int    i, j, loop;

double electric_field(int, int, int);
void infinite_length_cylinder();
void cylinder_edge();

int
main()
{
  double low = 0;
  double low0 = 3e-4;
  double eps0 = 8.854e-12;  
  // 変数の設定。
  a = 0.05; // 外側電極の半径[m]。a < RMAXとすること
  b = 2e-3; // 内側電極の半径[m]。b < aとすること
  t = 2e-3; // 外側電極の厚さ[m]
  Vmax = 10e3; // 内側電極(高電圧電極)の電位[V]

  // 発展課題の円筒電極端の計算でのみ必要となる変数。
  c = 0.01; // 中心電極端の金属球の半径。
  lin = 0.06;
  lout = 0.03;

  // 格子点の刻み幅h[m]
  h = RMAX/(double)N;

  /*
   *  配列変数phi, Es, flagの領域確保。
   *  プログラム中では見やすくするために
   *     PHI(i, j) = *(phi + i*(M+1) + j)
   *  のように表示を変更する。
   *
   *  PHI(i, j): 格子点(i, j)の電位
   *  Es(i, j):  格子点(i, j)の電界強度
   *  FLAG(i, j): 格子点(i, j)の境界条件の種類を表すフラグ
   */
  phi = (double *)calloc((N+1)*(M+1), sizeof(double));
  Es  = (double *)calloc((N+1)*(M+1), sizeof(double));
  flag = (int *)calloc((N+1)*(M+1), sizeof(int));

  // 各格子点の初期設定。各格子点の初期値とフラグを設定する。
  // 計算する対象に合わせて、どちらか1つコメントアウトする。
  infinite_length_cylinder(); // 無限長同軸円筒電極の計算。N=500, M=2など。
  //cylinder_edge();  // 円筒電極端の計算。N=300, M=300など。

  phi_max = fabs(Vmax); // 計算領域内の最大電位 (= 高圧電極の電位)
  error_limit = ERR_LIM; // 収束判定誤差を代入
  error_max = 1.;
  
  // 反復法によるPHI(i, j)の計算。
  // 全ての格子点電位の前回ループとの誤差がerror_limitを下回るまでループする。
  for (loop = 0; error_max > error_limit; loop++) {

    // 10000ループ毎に経過表示
    if(!(loop%1000))
      fprintf(stderr, "%5d  %.3e\n", loop, error_max);

    // 全ての格子点電位の前回ループとの誤差の最大値 error_max を初期化
    error_max = 0;
    
    // 格子点(i, j)をループ
    for (i = 0; i <= N; i++) {
      for (j = 0; j <= M; j++) {
	// 格子点(i, j)の座標(r, z)
	r = i*h;
	z = j*h;

	// 前回ループの PHI(i, j) を phi_prev にいれておく
	phi_prev = PHI(i, j);
	
    if(b < r && r < a){
        low = low0;
    } else {
        low = 0;
    }
	// 境界でない格子について、ラプラス方程式でPHI(i, j)を計算する
	if (FLAG(i, j) == 0){
	  PHI(i, j) = 0.25*(PHI(i+1, j) + PHI(i-1, j) + PHI(i, j+1)
			    + PHI(i, j-1)
			    + 0.5*h*(PHI(i+1, j) - PHI(i-1, j))/r
                + h*h*low/eps0);
    }
	// ディリクレ境界については PHI(i, j) の計算は不要
	else if (FLAG(i, j) == 1) {
	}

	// z方向ノイマン境界(上側境界)について PHI(i, j) を計算する
	else if (FLAG(i, j) == 2)
	  PHI(i, j) = 0.25*(PHI(i+1, j) + PHI(i-1, j) + 2*PHI(i, j-1)
			    + 0.5*h*(PHI(i+1, j) - PHI(i-1, j))/r
                + h*h*low/eps0);

	// z方向ノイマン境界(下側境界)について PHI(i, j) を計算する
	else if (FLAG(i, j) == 3)
	  PHI(i, j) = 0.25*(PHI(i+1, j) + PHI(i-1, j) + 2*PHI(i, j+1)
			    + 0.5*h*(PHI(i+1, j) - PHI(i-1, j))/r
                + h*h*low/eps0);

	// r中心軸上のノイマン境界について PHI(i, j) を計算する
	else if (FLAG(i, j) == 4)
	  PHI(i, j) = 1/6.*(PHI(i, j+1) + PHI(i, j-1) + 4*PHI(i+1, j)+ h*h*low/eps0);

	// r方向とz方向のノイマン条件の格子点の例外処理。
	else if (FLAG(i, j) == 5)
	  PHI(i, j) = 1/6.*(4*PHI(i+1, j) + 2*PHI(i, j-1)+ h*h*low/eps0);


    
	// 前回ループと新しい答えの差を phi_max で規格化
	error = (fabs(PHI(i, j) - phi_prev))/phi_max;
	
	// ループ内の最大誤差を常にerror_maxに持つようにする
	if (error_max < error)
	  error_max = error;
      }
    }
  }
  
  // 電位PHI(i, j)から電界強度ES(i, j)を差分法で計算。
  for (i = 0; i <= N; i++)
    for (j = 0; j <= M; j++)
      ES(i, j) = electric_field(i, j, 1); // 径方向電界を出力する場合
      //ES(i, j) = electric_field(i, j, 2); // 電界の大きさを出力する場合

  // ここから下は計算結果の出力。出力したい結果に応じて、適宜コメントアウトする。
  
  // z = const としたときの電位や電界強度の径方向分布出力。
  // 長さcm、電位kV、電界kV/cmとしている。
  j = 1; // 0-Mの任意の値を設定する。
  float phi_th = 0;
  float err= 0;
  for (i = 0; i <= N; i++) {
 

//課題その2
  //  phi_th = Vmax/(1e3)*log(a*100/i*h*1e2)/log(a/b);
    //if(b < i*h && i*h < a){
    //  err += fabs(phi_th - PHI(i,j))/PHI(i,j);
    //}
  
   printf("%.5lg\t%.5lg\t%.5lg\r\n", i*h*1e2, PHI(i, j)*1e-3,ES(i, j)*1e-5);
    //printf("%.5lg\t%.5lg\r\n", i*h*1e2, ES(i, j)*1e-5);

  }

  //printf("%.5lg\t%.5lg\r\n",ERR_LIM,err);

  // GraphR のフォーマットでポテンシャル出力。長さcm、電位kV、電界kV/cmとしている。
  /*printf("データ形式,2,\r\n,,\r\n,,\r\n");
  for (i = 0; i <= N; i++)
    for (j = 0; j <= M; j++)
      //printf("%.5lg,%.5lg,%.5lg\r\n", i*h*1e2, j*h*1e2, PHI(i, j)*1e-3);
      printf("%.5lg,%.5lg,%.5lg\r\n", i*h*1e2, j*h*1e2, ES(i, j)*1e-5);*/
}

/*
 *  electric_field(i, j, sw) --
 *    格子点(i, j)の電界を返す。電界ベクトルを(Er, Ez)とおくと、
 *    sw = 1のとき：　径方向の電界 Er を返す。r方向を正とする。
 *    sw = 2のとき：　電界の大きさ sqrt(Er*Er + Ez*Ez) を返す。
 */
double
electric_field(int i, int j, int sw)
{
  // r方向電界 Er
  if (i == 0)
    Er = -(PHI(i+1, j) - PHI(i, j))/h;
  else if (i == N)
    Er = -(PHI(i, j) - PHI(i-1, j))/h;
  else
    Er = -(PHI(i+1, j) - PHI(i-1, j))/(2*h);
  
  // z方向電界 Ez
  if (j == 0)
    Ez = -(PHI(i, j+1) - PHI(i, j))/h;
  else if (j == M)
    Ez = -(PHI(i, j) - PHI(i, j-1))/h;
  else
    Ez = -(PHI(i, j+1) - PHI(i, j-1))/(2*h);

  if (sw == 1)
    return Er;
  else
    return sqrt(Er*Er + Ez*Ez);
}

/*
 *  infinite_length_cylinder() --
 *    無限長同軸円筒電極の計算の初期設定。
 *    各格子点の電位PHI(i, j)の初期値とフラグFLAG(i, j)を設定する。
 */
void
infinite_length_cylinder()
{
  for (i = 0; i <= N; i++)
    for (j = 0; j <= M; j++) {
      // 格子点(i, j)の座標(r, z)
      r = i*h;
      z = j*h;
      
      // 内側電極(高圧電極)の境界条件設定(ディリクレ条件: V = Vmax)
      if (r <= b) {
	FLAG(i, j) = 1; // ディリクレ条件の格子点は FLAG = 1
	PHI(i, j) = Vmax;
      }
      // 外側電極(接地電極)のおよびr方向外側境界の境界条件設定(ディリクレ条件: V = 0)
      else if ((r >= a && r <= a + t) || i == N) {
	FLAG(i, j) = 1;
	PHI(i, j) = 0.0;
      }
      // z方向ノイマン条件(上側境界)
      else if (j == M) {
	FLAG(i, j) = 2; // z方向ノイマン条件(上側境界)の格子点は FLAG = 2
	PHI(i, j) = PHI0;
      }
      // z方向ノイマン条件(下側境界)
      else if (j == 0) {
	FLAG(i, j) = 3; // z方向ノイマン条件(下側境界)の格子点は FLAG = 3
	PHI(i, j) = PHI0;
      }
      // 境界条件のない格子点
      else {
	FLAG(i, j) = 0; // 境界条件のない格子点は FLAG = 0
	PHI(i, j) = PHI0;
      }
    }
}

/*
 *  cylinder_edge() --
 *    円筒電極端の計算の初期設定。
 *    各格子点の電位PHI(i, j)の初期値とフラグFLAG(i, j)を設定する。
 */
void
cylinder_edge()
{
  for (i = 0; i <= N; i++)
    for (j = 0; j <= M; j++) {
      // 格子点(i, j)の座標(r, z)
      r = i*h;
      z = j*h;
      
      // 内側電極(高圧電極)の境界条件設定(ディリクレ条件: V = Vmax)
      if ((r <= b && z <= lin) || r*r + (z-lin)*(z-lin) <= c*c) {
	FLAG(i, j) = 1; // ディリクレ条件の格子点は FLAG = 1
	PHI(i, j) = Vmax;
      }
      // 外側電極(接地電極)のおよびr方向外側境界の境界条件設定(ディリクレ条件: V = 0)
      else if ((r >= a && r <= a + t && z <= lout) || i == N) {
	FLAG(i, j) = 1;
	PHI(i, j) = 0.0;
      }
      // z方向ノイマン条件(上側境界)
      else if (j == M) {
	FLAG(i, j) = 2; // z方向ノイマン条件(上側境界)の格子点は FLAG = 2
	PHI(i, j) = PHI0;
      }
      // z方向ノイマン条件(下側境界)
      else if (j == 0) {
	FLAG(i, j) = 3; // z方向ノイマン条件(下側境界)の格子点は FLAG = 3
	PHI(i, j) = PHI0;
      }
      // r中心軸上のノイマン条件
      else if (i == 0) {
	FLAG(i, j) = 4; // r中心軸上のノイマン条件の格子点は FLAG = 4
	PHI(i, j) = PHI0;
      }
      // 境界条件のない格子点
      else {
	FLAG(i, j) = 0; // 境界条件のない格子点は FLAG = 0
	PHI(i, j) = PHI0;
      }
      // r方向とz方向のノイマン条件の格子点の例外処理。
      FLAG(0, M) = 5;
    }
}
