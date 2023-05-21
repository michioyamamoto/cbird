//☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
//     Sparse Reduced Latent Class AnalysisのALSによる最適化
//
//  ファイル名：OptimSRLCA.c
//  ファイル内容：
//  作成者：YAMAMOTO, Michio
//  作成日：2013年4月12日
//  最終更新日：2017年02月06日
//  コメント：ペナルティ項をsubgradientを利用して最適化
//           個人の得点に直交制約の有無を選択できるようにした（140303）
//           GPアルゴリズムにおけるn_ite_GPの初期化忘れを修正（150114）
//           パッケージ公開に関連して関数名等を変更した（170206）
//           - SRLCA->CBIRD, EstimateScores->EstScore
//☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Parse.h>
#include <R_ext/Lapack.h>

void LogDensity(double *Q, double *F, double *A, double *MU, int N_sub, int N_var, int N_comp, int N_clust, double *LogDens);
void F_Identity(int N, double *I);
double F_LossCBIRD(double *F, double *A, double *U, double *mu, double *X, double lambda, int N_sub, int N_var, int N_comp, int N_clust);
double PenLogLikelihood(double *Q, double *F, double *A, double *MU, double *g, int N_sub, int N_var, int N_comp, int N_clust, double lambda);
void F_Grad(double *F, double *A, double *mu, double *X_s, double *G, int N_sub, int N_var, int N_comp, int N_clust, double *N_k);
void F_Grad_ind(double *F, int row1, int col1, double *A, int row2, int col2, double *mu, int length1, double *X_s, int row3, int col3, double *G);

void F_Identity2(int N, double *I);
double F_LossSLPCA(double *F, int row1, int col1, double *A, int row2, int col2, double *mu, int length1, double *X, int row3, int col3);

/*
  INFO: c(N.sub, N.var, N.comp, N.ite, N.ite.GP, N.ite.alpha, N.clust)
 */

SEXP EstScore_Oblique_C(SEXP R_X, SEXP R_F, SEXP R_A, SEXP R_MU, SEXP R_Q, SEXP INFO)
{
  int    N_sub = (int) REAL(INFO)[0];
  int    N_var = (int) REAL(INFO)[1];
  int    N_comp = (int) REAL(INFO)[2];
  int    N_ite = (int) REAL(INFO)[3];
  int    N_ite_GP = (int) REAL(INFO)[4];
  int    N_ite_alpha = (int) REAL(INFO)[5];

  double alpha, alpha_ini = 100;
  double blas_one = 1.0;
  double blas_zero = 0.0;
  int    blas_one_int = 1;
  double diff_loss;
  double eps = 0.001;
  double *F_current;
  double *F_old;
  double *F_target;
  double *G;
  int    i, j;
  double I_comp[N_comp * N_comp];
  int    info;
  int    ipiv[N_comp];
  double loss_old, loss = 1.0E+100, loss_GP;
  double Mat_comp_comp[N_comp * N_comp];
  double *Mat_sv;
  int    n_ite = 0;
  double N_sub_double = (double) N_sub;
  double scalar;
  double sum = 0.0;
  double temp;
  double vec_1_sub[N_sub];
  double *Z;

  /* double *rans; */
  SEXP ans;
  SEXP Theta;
  SEXP Z_s;
  SEXP ans_ind;


  F_current = (double *) malloc(sizeof(double) * N_sub * N_comp);
  F_old     = (double *) malloc(sizeof(double) * N_sub * N_comp);
  F_target  = (double *) malloc(sizeof(double) * N_sub * N_comp);
  G         = (double *) malloc(sizeof(double) * N_sub * N_comp);
  Mat_sv    = (double *) malloc(sizeof(double) * N_sub * N_var);
  Z         = (double *) malloc(sizeof(double) * N_sub * N_var);

  PROTECT(Theta = allocVector(REALSXP, N_sub * N_var));
  PROTECT(Z_s = allocVector(REALSXP, N_sub * N_var));
  PROTECT(ans = allocVector(VECSXP, 2));
  PROTECT(ans_ind = allocVector(REALSXP, 2));

  // vec.1.sub を求めておく
  for (i = 0; i < N_sub; i++)
	vec_1_sub[i] = 1.0;

  F_Identity2(N_comp, I_comp); // 単位行列を求めておく

  /*******************************/
  /*     ALS algorithm start     */
  /*******************************/
  do {
  	n_ite++;

  	loss_old = loss;

  	/***** Theta の計算 *****/
	for (i = 0; i < (N_sub * N_var); i++)
	  REAL(Theta)[i] = 0.0;
  	dger_(&N_sub, &N_var, &blas_one, vec_1_sub, &blas_one_int, REAL(R_MU), &blas_one_int, REAL(Theta), &N_sub); // vec.1.sub %*% t(mu)
	dgemm_("N", "T", &N_sub, &N_var, &N_comp, &blas_one, REAL(R_F), &N_sub, REAL(R_A), &N_var, &blas_one, REAL(Theta), &N_sub); // Theta <- F %*% t(A) + vec.1.sub %*% t(mu)


  	/***** Z の計算 *****/
  	for (i = 0; i < (N_sub * N_var); i++)
  	  Z[i] = REAL(Theta)[i] + 4 * REAL(R_Q)[i] * (1 - 1 / (1 + exp(-1.0 * REAL(R_Q)[i] * REAL(Theta)[i])));


  	/***** F の更新  *****/
	for (i = 0; i < (N_sub * N_var); i++)
	  Mat_sv[i] = 0.0;
  	dger_(&N_sub, &N_var, &blas_one, vec_1_sub, &blas_one_int, REAL(R_MU), &blas_one_int, Mat_sv, &N_sub); // vec.1.sub %*% t(mu)
  	for (i = 0; i < (N_sub * N_var); i++)
  	  REAL(Z_s)[i] = Z[i] - Mat_sv[i];

	dgemm_("T", "N", &N_comp, &N_comp, &N_var, &blas_one, REAL(R_A), &N_var, REAL(R_A), &N_var, &blas_zero, Mat_comp_comp, &N_comp); // t(A) %*% A
	dgesv_(&N_comp, &N_comp, Mat_comp_comp, &N_comp, ipiv, I_comp, &N_comp, &info); // solve(t(A) %*% A)
	dgemm_("N", "N", &N_sub, &N_comp, &N_var, &blas_one, REAL(Z_s), &N_sub, REAL(R_A), &N_var, &blas_zero, Mat_sv, &N_sub); // Z_s %*% A
	dgemm_("N", "N", &N_sub, &N_comp, &N_comp, &blas_one, Mat_sv, &N_sub, I_comp, &N_comp, &blas_zero, REAL(R_F), &N_sub); // Z_s %*% A %*% solve(t(A) %*% A)

  	/* 収束チェック用の損失関数の計算  */
  	loss = F_LossSLPCA(REAL(R_F), N_sub, N_comp, REAL(R_A), N_var, N_comp, REAL(R_MU), N_var, Z, N_sub, N_var);
  	diff_loss = fabs(loss_old - loss);

  } while (diff_loss > eps && n_ite < N_ite); // End of ALS algorithm


  /* 解の受け渡し */
  REAL(ans_ind)[0] = (double) n_ite;
  REAL(ans_ind)[1] = loss;

  SET_VECTOR_ELT(ans, 0,     R_F);
  SET_VECTOR_ELT(ans, 1, ans_ind);

  UNPROTECT(4);
  free(F_current);
  free(F_old);
  free(F_target);
  free(G);
  free(Mat_sv);
  free(Z);

  return(ans);
}


SEXP EstScore_Orthogonal_C(SEXP R_X, SEXP R_F, SEXP R_A, SEXP R_MU, SEXP R_Q, SEXP INFO)
{
  int     N_sub = (int) REAL(INFO)[0];
  int     N_var = (int) REAL(INFO)[1];
  int     N_comp = (int) REAL(INFO)[2];
  int     N_ite = (int) REAL(INFO)[3];
  int     N_ite_GP = (int) REAL(INFO)[4];
  int     N_ite_alpha = (int) REAL(INFO)[5];

  double  alpha, alpha_ini = 100;
  double  blas_one = 1.0;
  double  blas_zero = 0.0;
  int     blas_one_int = 1;
  double  diff_loss;
  double  diff_loss_GP;
  double  eps = 0.001;
  double *F_current;
  double *F_old;
  double *F_target;
  double *G;
  int     i, j;
  double  I_comp[N_comp * N_comp];
  int     info;
  int     ipiv[N_comp];
  char    jobu = 'S', jobvt = 'S';
  double  loss_old, loss = 1.0E+100, loss_GP, loss_current_GP;
  int     lwork;
  double  Mat_comp_comp[N_comp * N_comp];
  double *Mat_sv;
  int     n_ite = 0;
  int     n_ite_alpha;
  int     n_ite_GP = 0;
  double  N_sub_double = (double) N_sub;
  double  scalar;
  double  singular_value[N_comp];
  double  sum = 0.0;
  double  temp;
  double *U_svd;
  double *vec_1_sub;
  double *Vt_svd;
  double *work;
  double *Z;
  double *Z_s;

  /* double *rans; */
  SEXP ans;
  SEXP Theta;
  SEXP ans_ind;

  // malloc によるメモリ割り当て
  F_current = (double *) malloc(sizeof(double) * N_sub * N_comp);
  F_old     = (double *) malloc(sizeof(double) * N_sub * N_comp);
  F_target  = (double *) malloc(sizeof(double) * N_sub * N_comp);
  G         = (double *) malloc(sizeof(double) * N_sub * N_comp);
  Mat_sv    = (double *) malloc(sizeof(double) * N_sub * N_var);
  U_svd     = (double *) malloc(sizeof(double) * N_sub * N_comp);
  vec_1_sub = (double *) malloc(sizeof(double) * N_sub);
  Vt_svd    = (double *) malloc(sizeof(double) * N_comp * N_comp);
  work      = (double *) malloc(sizeof(double) * N_sub * N_sub);
  Z         = (double *) malloc(sizeof(double) * N_sub * N_var);
  Z_s       = (double *) malloc(sizeof(double) * N_sub * N_var);

  // PROTECT によるメモリ割り当て
  PROTECT(Theta = allocVector(REALSXP, N_sub * N_var));
  PROTECT(ans = allocVector(VECSXP, 2));
  PROTECT(ans_ind = allocVector(REALSXP, 2));

  // vec.1.sub を求めておく
  for (i = 0; i < N_sub; i++)
	vec_1_sub[i] = 1.0;

  // 単位行列を求めておく
  F_Identity2(N_comp, I_comp);

  // dgesvd関数用の変数
  lwork = N_sub * N_sub;

  /*******************************/
  /*     ALS algorithm start     */
  /*******************************/
  do {
  	n_ite++;

  	loss_old = loss;

  	/***** Theta の計算 *****/
	for (i = 0; i < (N_sub * N_var); i++)
	  REAL(Theta)[i] = 0.0;
  	dger_(&N_sub, &N_var, &blas_one, vec_1_sub, &blas_one_int, REAL(R_MU), &blas_one_int, REAL(Theta), &N_sub); // vec.1.sub %*% t(mu)
	dgemm_("N", "T", &N_sub, &N_var, &N_comp, &blas_one, REAL(R_F), &N_sub, REAL(R_A), &N_var, &blas_one, REAL(Theta), &N_sub); // Theta <- F %*% t(A) + vec.1.sub %*% t(mu)


  	/***** Z の計算 *****/
  	for (i = 0; i < (N_sub * N_var); i++)
  	  Z[i] = REAL(Theta)[i] + 4 * REAL(R_Q)[i] * (1 - 1 / (1 + exp(-1.0 * REAL(R_Q)[i] * REAL(Theta)[i])));


  	/***** F の更新  *****/
	for (i = 0; i < (N_sub * N_var); i++)
	  Mat_sv[i] = 0.0;
  	dger_(&N_sub, &N_var, &blas_one, vec_1_sub, &blas_one_int, REAL(R_MU), &blas_one_int, Mat_sv, &N_sub); // vec.1.sub %*% t(mu)
  	for (i = 0; i < (N_sub * N_var); i++)
  	  Z_s[i] = Z[i] - Mat_sv[i];
  	loss_GP = F_LossSLPCA(REAL(R_F), N_sub, N_comp, REAL(R_A), N_var, N_comp, REAL(R_MU), N_var, Z, N_sub, N_var);

	n_ite_GP = 0;
  	do {
  	  n_ite_GP++;

  	  for (i = 0; i < (N_sub * N_comp); i++)
  		F_current[i] = REAL(R_F)[i];

  	  //Calculation of Gradient of Loss function at F_ind
  	  F_Grad_ind(REAL(R_F), N_sub, N_comp, REAL(R_A), N_var, N_comp, REAL(R_MU), N_var, Z_s, N_sub, N_var, G);

  	  /* 最適なalphaの探索 */
  	  // 現在の目的関数の値を計算しておく
  	  n_ite_alpha = 0;
  	  alpha = 2 * alpha_ini;
	  loss_current_GP = F_LossSLPCA(F_current, N_sub, N_comp, REAL(R_A), N_var, N_comp, REAL(R_MU), N_var, Z, N_sub, N_var);

  	  do {
  	  	n_ite_alpha++;
  	  	alpha = alpha / 2;
  	  	scalar = -1 * alpha;

  	  	/* Fの更新 */
  	  	for (i = 0; i < (N_sub * N_comp); i++)
  	  	  F_target[i] = F_current[i];

  	  	dgemm_("N", "N", &N_sub, &N_comp, &N_comp, &scalar, G, &N_sub, I_comp, &N_comp, &blas_one, F_target, &N_sub); // F - alpha * G
  	  	dgesvd_(&jobu, &jobvt, &N_sub, &N_comp, F_target, &N_sub, singular_value, U_svd, &N_sub, Vt_svd, &N_comp, work, &lwork, &info);
  	  	if (info != 0)
  	  	  error(("error code %d from Lapack routine '%s'"), info, "dgesvd");
  	  	dgemm_("N", "N", &N_sub, &N_comp, &N_comp, &blas_one, U_svd, &N_sub, Vt_svd, &N_comp, &blas_zero, REAL(R_F), &N_sub); // F <- U %*% t(V)

  	  	/* 更新値の目的関数の値を計算しておく */
		loss_GP = F_LossSLPCA(REAL(R_F), N_sub, N_comp, REAL(R_A), N_var, N_comp, REAL(R_MU), N_var, Z, N_sub, N_var);

  	  } while (loss_GP > loss_current_GP && n_ite_alpha < N_ite_alpha); // alpha の探索終了

  	  //alpha = 0の場合，元の F を代入しておく
  	  if (n_ite_alpha == N_ite_alpha) {
  	  	alpha = 0;
  	  	for (i = 0; i < (N_sub * N_comp); i++)
  	  	  REAL(R_F)[i] = F_current[i];
  	  	loss_GP = loss_current_GP;
  	  }

  	  /* 収束判定 */
  	  diff_loss_GP = loss_current_GP - loss_GP;

  	} while (diff_loss_GP > eps && n_ite_GP < N_ite_GP); // End of GP algorithm

  	/* 収束チェック用の損失関数の計算  */
  	loss = F_LossSLPCA(REAL(R_F), N_sub, N_comp, REAL(R_A), N_var, N_comp, REAL(R_MU), N_var, Z, N_sub, N_var);
  	diff_loss = fabs(loss_old - loss);

  } while (diff_loss > eps && n_ite < N_ite); // End of ALS algorithm


  /* 解の受け渡し */
  REAL(ans_ind)[0] = (double) n_ite;
  REAL(ans_ind)[1] = loss;

  SET_VECTOR_ELT(ans, 0,     R_F);
  SET_VECTOR_ELT(ans, 1, ans_ind);

  UNPROTECT(3);
  free(F_current);
  free(F_old);
  free(F_target);
  free(G);
  free(Mat_sv);
  free(U_svd);
  free(vec_1_sub);
  free(Vt_svd);
  free(work);
  free(Z);
  free(Z_s);

  return(ans);
}


SEXP OptimCBIRD_C(SEXP Y, SEXP F, SEXP A, SEXP g, SEXP U, SEXP MU, SEXP Q, SEXP INFO, SEXP LAMBDA, SEXP PLL_RECORD)
{
  // Rからの引渡し
  int     N_sub = (int) REAL(INFO)[0];
  int     N_var = (int) REAL(INFO)[1];
  int     N_comp = (int) REAL(INFO)[2];
  int     N_ite = (int) REAL(INFO)[3];
  int     N_ite_GP = (int) REAL(INFO)[4];
  int     N_ite_alpha = (int) REAL(INFO)[5];
  int     N_clust = (int) REAL(INFO)[6];
  double  eps = (double) REAL(INFO)[7];
  double  lambda = (double) REAL(LAMBDA)[0];

  //変数定義
  double *A_target;
  double *A_target_inv;
  double *A_target_var;
  double  alpha, alpha_ini = 100;
  double  blas_one = 1.0, blas_zero = 0.0, blas_minus_one = 1.0;
  int     blas_one_int = 1;
  double  C;
  //  double *D_mat;
  double  diff_loss;
  double  diff_loss_GP;
  double *eigen_value;
  double *eigen_value2; //昇順から降順への並び替え用
  double *F_current;
  double *F_s;
  double *F_target;
  int     flag;
  double *G;
  int     i, j, k, l, l2;
  double *I_comp;
  int     info;
  char    jobu = 'S', jobvt = 'S'; // for dgesvd_
  char    jobz = 'V', uplo = 'U'; // for dsyev_
  double  loss = 1.0E+100, loss_old;
  double  loss_GP, loss_current_GP;
  double *LogDens;
  int     lwork;
  double *Mat_sv;
  double *Mat_cv;
  double  max_C;
  double  max_density;
  int     n_ite = 0;
  int     n_ite_alpha;
  int     n_ite_GP = 0;
  int     N_clust_clust = N_clust * N_clust;
  int     N_clust_comp = N_clust * N_comp;
  int     N_clust_var = N_clust * N_var;
  int     N_comp_comp = N_comp * N_comp;
  double  N_k[N_clust];
  double  N_sub_double = (double) N_sub;
  int     N_sub_var = N_sub * N_var;
  int     N_sub_clust = N_sub * N_clust;
  int     N_sub_comp = N_sub * N_comp;
  int     N_sub_sub = N_sub * N_sub;
  int     N_var_comp = N_var * N_comp;
  double  scalar;
  double  sign_C;
  double  singular_value[N_comp];
  double  sum;
  double  temp;
  double  temp_double;
  double *Theta;
  double *U_old;
  double *U_svd;
  double *vec_1_clust;
  double  vec_comp[N_comp];
  double *Vt_svd;
  double *W;
  double *work;
  double *X;
  double *X_par;
  double *Z;
  double *Z_s;

  SEXP    ans;
  SEXP    ans_F;
  SEXP    ans_A;
  SEXP    ans_g;
  SEXP    ans_U;
  SEXP    ans_MU;
  SEXP    ans_ind;
  SEXP    ans_record;


  // malloc によるメモリ割り当て
  A_target    = (double *) malloc(sizeof(double) * N_comp_comp);
  A_target_inv= (double *) malloc(sizeof(double) * N_comp_comp);
  A_target_var= (double *) malloc(sizeof(double) * N_comp_comp);
  //  D_mat       = (double *) malloc(sizeof(double) * N_var_comp);
  eigen_value = (double *) malloc(sizeof(double) * N_comp);
  eigen_value2= (double *) malloc(sizeof(double) * N_comp);
  F_current   = (double *) malloc(sizeof(double) * N_sub_var);
  F_target    = (double *) malloc(sizeof(double) * N_sub_var);
  F_s         = (double *) malloc(sizeof(double) * N_var_comp);
  G           = (double *) malloc(sizeof(double) * N_clust_comp);
  I_comp      = (double *) malloc(sizeof(double) * N_comp_comp);
  LogDens     = (double *) malloc(sizeof(double) * N_sub_clust);
  Mat_sv      = (double *) malloc(sizeof(double) * N_sub_var);
  Mat_cv      = (double *) malloc(sizeof(double) * N_clust_var);
  Theta       = (double *) malloc(sizeof(double) * N_clust_var);
  U_old       = (double *) malloc(sizeof(double) * N_sub_clust);
  U_svd       = (double *) malloc(sizeof(double) * N_sub_comp);
  vec_1_clust = (double *) malloc(sizeof(double) * N_clust);
  Vt_svd      = (double *) malloc(sizeof(double) * N_comp_comp);
  W           = (double *) malloc(sizeof(double) * N_comp_comp);
  work        = (double *) malloc(sizeof(double) * N_sub_sub);
  X           = (double *) malloc(sizeof(double) * N_sub*N_var*N_clust);
  X_par       = (double *) malloc(sizeof(double) * N_clust_var);
  Z           = (double *) malloc(sizeof(double) * N_sub_var);
  Z_s         = (double *) malloc(sizeof(double) * N_sub_var);


  // PROTECT によるメモリ割り当て
  //  PROTECT(test_ret = allocVector(REALSXP, N_var_comp));
  PROTECT(ans = allocVector(VECSXP, 7));
  PROTECT(ans_F = allocMatrix(REALSXP, N_clust, N_comp));
  PROTECT(ans_A = allocMatrix(REALSXP, N_var, N_comp));
  PROTECT(ans_g = allocVector(REALSXP, N_clust));
  PROTECT(ans_U = allocMatrix(REALSXP, N_sub, N_clust));
  PROTECT(ans_MU = allocVector(REALSXP, N_var));
  PROTECT(ans_ind = allocVector(REALSXP, 2));
  PROTECT(ans_record = allocVector(REALSXP, N_ite));


  // testでの最終的な解を返す変数
  //  rans = REAL(test_ret);

  // vec.1.clust を求めておく
  for (i = 0; i < N_clust; i++)
	vec_1_clust[i] = 1.0;

  // 単位行列を求めておく
  F_Identity(N_comp, I_comp);

  // dgesvd関数用の変数
  lwork = N_sub_sub;

  /*******************************/
  /*     ALS algorithm start     */
  /*******************************/
  do {
  	n_ite++;
  	loss_old = loss;

	for (i = 0; i < N_sub_clust; i++)
	  U_old[i] = REAL(U)[i];


  	/***** D_mat の計算 *****/
  	/* for (i = 0; i < (N_var_comp); i++) { */
  	/*   temp = lambda / (2 * fabs(REAL(A)[i])); */
  	/*   if (temp > DBL_MAX) { */
  	/* 	D_mat[i] = 1.0E+200; */
  	/*   } else { */
  	/* 	D_mat[i] = temp; */
  	/*   } */
  	/* } */


	//---------------------------------
	//
	//          E step
	//
	//  U について条件付き期待値をとる
	//---------------------------------
	// Y.n の k を与えたもとでの条件付き密度を求める
	LogDensity(REAL(Q), REAL(F), REAL(A), REAL(MU), N_sub, N_var, N_comp, N_clust, LogDens);

	// U についての条件付き期待値を計算
	max_density = LogDens[0];
	for (i = 1; i < N_sub_clust; i++)
	  if (max_density < LogDens[i])
		max_density = LogDens[i];

	for (i = 0; i < N_sub; i++) {
	  sum = 0.0;
	  for (k = 0; k < N_clust; k++)
		sum = sum + exp(log(REAL(g)[k]) + LogDens[i + N_sub * k] - max_density);
	  for (k = 0; k < N_clust; k++)
		REAL(U)[i + N_sub * k] = exp(log(REAL(g)[k]) + LogDens[i + N_sub * k] - max_density) / sum;
	}

	for (k = 0; k < N_clust; k++) {
	  sum = 0.0;
	  for (i = 0; i < N_sub; i++)
		sum = sum + REAL(U)[i + N_sub * k];
	  N_k[k] = sum;
	}

	// 空のクラスターが生じた場合は終了する
	flag = 0;
	for (k = 0; k < N_clust; k++)
	  if (N_k[k] < 1.0E-200)
		flag = 1;
	if (flag == 1)
	  break;

	//---------------------------------
	//
	//          M step
	//
	//---------------------------------
	/***** Update g *****/
	sum = 0.0;
	for (k = 0; k < (N_clust - 1); k++) {
	  REAL(g)[k] = N_k[k] / N_sub_double;
	  sum = sum + REAL(g)[k];
	}
	REAL(g)[N_clust - 1] = 1 - sum;


	/***** Theta の計算 *****/
	for (i = 0; i < N_clust_var; i++)
	  Theta[i] = 0.0;
  	dger_(&N_clust, &N_var, &blas_one, vec_1_clust, &blas_one_int, REAL(MU), &blas_one_int, Theta, &N_clust); // vec.1.clust %*% t(mu)
	dgemm_("N", "T", &N_clust, &N_var, &N_comp, &blas_one, REAL(F), &N_clust, REAL(A), &N_var, &blas_one, Theta, &N_clust); // Theta <- F %*% t(A) + vec.1.clust %*% t(mu)


  	/***** X の計算 *****/
	for (i = 0; i < N_sub; i++)
	  for (k = 0; k < N_clust; k++)
		  for (j = 0; j < N_var; j++)
		  X[i + N_sub * k + N_sub * N_clust * j] = Theta[k + N_clust * j] + 4 * REAL(Q)[i + N_sub * j] * (1 - 1 / (1 + exp(-REAL(Q)[i + N_sub * j] * Theta[k + N_clust * j])));


  	/***** mu の更新  *****/
	for (k = 0; k < N_clust; k++) {
	  for (j = 0; j < N_var; j++) {
		sum = 0.0;
		for (i = 0; i < N_sub; i++)
		  sum = sum + X[i + N_sub * k + N_sub * N_clust * j] * REAL(U)[i + N_sub * k];
		X_par[k + N_clust * j] = sum / N_k[k];
	  }
	}
	dgemm_("N", "T", &N_clust, &N_var, &N_comp, &blas_one, REAL(F), &N_clust, REAL(A), &N_var, &blas_zero, Mat_cv, &N_clust); // F %*% t(A)
  	for (j = 0; j < N_var; j++) {
  	  sum = 0.0;
  	  for (k = 0; k < N_clust; k++)
  		sum = sum + N_k[k] * (X_par[k + N_clust * j] - Mat_cv[k + N_clust * j]);
  	  REAL(MU)[j] = sum / N_sub_double;
  	}

  	/***** F の更新  *****/
	for (i = 0; i < N_clust_var; i++)
	  Mat_cv[i] = 0.0;
  	dger_(&N_clust, &N_var, &blas_one, vec_1_clust, &blas_one_int, REAL(MU), &blas_one_int, Mat_cv, &N_clust); // vec.1.clust %*% t(mu)
	for (k = 0; k < N_clust; k++) {
	  for (j = 0; j < N_var; j++) {
		sum = 0.0;
		for (i = 0; i < N_sub; i++)
		  sum = sum + (X[i + N_sub * k + N_sub * N_clust * j] - Mat_cv[k + N_clust * j]) * REAL(U)[i + N_sub * k];
		X_par[k + N_clust * j] = sum / N_k[k];
	  }
	}

	//	loss_GP = F_LossCBIRD(REAL(F), REAL(A), REAL(U), REAL(MU), X, D_mat, N_sub, N_var, N_comp, N_clust);
	loss_GP = F_LossCBIRD(REAL(F), REAL(A), REAL(U), REAL(MU), X, lambda, N_sub, N_var, N_comp, N_clust);

	n_ite_GP = 0;
  	do {
  	  n_ite_GP++;

  	  for (i = 0; i < N_clust_comp; i++)
  		F_current[i] = REAL(F)[i];

  	  //Calculation of Gradient of Loss function at F
  	  F_Grad(REAL(F), REAL(A), REAL(MU), X_par, G, N_sub, N_var, N_comp, N_clust, N_k);

  	  // 最適なalphaの探索
  	  // 現在の目的関数の値を計算しておく
  	  n_ite_alpha = 0;
  	  alpha = 2 * alpha_ini;
	  // loss_current_GP = F_LossCBIRD(F_current, REAL(A), REAL(U), REAL(MU), X, D_mat, N_sub, N_var, N_comp, N_clust);
	  loss_current_GP = F_LossCBIRD(F_current, REAL(A), REAL(U), REAL(MU), X, lambda, N_sub, N_var, N_comp, N_clust);

  	  do {
  	  	n_ite_alpha++;
  	  	alpha = alpha / 2;
  	  	scalar = -1 * alpha;

  	  	for (i = 0; i < N_clust_comp; i++)
  	  	  F_target[i] = F_current[i];

		dgemm_("N", "N", &N_clust, &N_comp, &N_comp, &scalar, G, &N_clust, I_comp, &N_comp, &blas_one, F_target, &N_clust); // F - alpha * G
		dgesvd_(&jobu, &jobvt, &N_clust, &N_comp, F_target, &N_clust, singular_value, U_svd, &N_clust, Vt_svd, &N_comp, work, &lwork, &info);
  	  	if (info != 0)
  	  	  error(("error code %d from Lapack routine '%s'"), info, "dgesvd");
  	  	dgemm_("N", "N", &N_clust, &N_comp, &N_comp, &blas_one, U_svd, &N_clust, Vt_svd, &N_comp, &blas_zero, REAL(F), &N_clust); // F <- U %*% t(V)

  	  	/* 更新値の目的関数の値を計算しておく */
		//		loss_GP = F_LossCBIRD(REAL(F), REAL(A), REAL(U), REAL(MU), X, D_mat, N_sub, N_var, N_comp, N_clust);
		loss_GP = F_LossCBIRD(REAL(F), REAL(A), REAL(U), REAL(MU), X, lambda, N_sub, N_var, N_comp, N_clust);

  	  } while (loss_GP > loss_current_GP && n_ite_alpha < N_ite_alpha); // alpha の探索終了

  	  //alpha = 0の場合，元の F を代入しておく
  	  if (n_ite_alpha == N_ite_alpha) {
  	  	alpha = 0;
  	  	for (i = 0; i < (N_clust_comp); i++)
  	  	  REAL(F)[i] = F_current[i];
  	  	loss_GP = loss_current_GP;
  	  }

  	  /* 収束判定 */
  	  diff_loss_GP = loss_current_GP - loss_GP;

  	} while (diff_loss_GP > eps && n_ite_GP < N_ite_GP); // End of GP algorithm



  	/***** A の更新  *****/
	// W = sum_k N_k f_k f_k' を計算しておく
	for (i = 0; i < N_comp; i++) {
	  for (j = 0; j < N_comp; j++) {
		sum = 0.0;
		for (k = 0; k < N_clust; k++)
		  sum = sum + REAL(F)[k + N_clust * i] * REAL(F)[k + N_clust * j] * N_k[k];
		W[i + N_comp * j] = sum;
	  }
	}

	// F.s を計算しておく
	for (j = 0; j < N_var; j++) {
	  for (l = 0; l < N_comp; l++) {
		sum = 0.0;
		for (i = 0; i < N_sub; i++) {
		  for (k = 0; k < N_clust; k++) {
			sum = sum + REAL(U)[i + N_sub * k] * (X[i + N_sub * k + N_sub * N_clust * j] - REAL(MU)[j]) * REAL(F)[k + N_clust * l];
		  }
		}
		F_s[j + N_var * l] = sum;
	  }
	}

	// 変数ごとに推定値を求める
	scalar = 4 * N_sub_double * lambda;

	for (j = 0; j < N_var; j++) {
	  for (l = 0; l < N_comp; l++) {
		C = 0.0;
		for (l2 = 0; l2 < N_comp; l2++)
		  C = C - W[l + N_comp * l2] * REAL(A)[j + N_var * l2];
		C = C + W[l + N_comp * l] * REAL(A)[j + N_var * l] + F_s[j + N_var * l];
		sign_C = C < 0.0 ? -1.0 : 1.0; // sign_C <- sign(C)
		temp_double = fabs(C) - scalar;
		max_C = 0.0 > temp_double ? 0.0 : temp_double; // max(0, abs(C) - 4 * N.sub * lambda)
		REAL(A)[j + N_var * l] = 1 / W[l + N_comp * l] * sign_C * max_C;
	  }
	}


  	/* 収束チェック用の罰則付き対数尤度の計算  */
	loss = PenLogLikelihood(REAL(Q), REAL(F), REAL(A), REAL(MU), REAL(g), N_sub, N_var, N_comp, N_clust, lambda);

	REAL(PLL_RECORD)[n_ite - 1] = loss;

  	diff_loss = fabs(loss_old - loss);

  } while (diff_loss > eps && n_ite < N_ite); // End of ALS algorithm

  /* for (i = 0; i < N_var_comp; i++) */
  /* 	rans[i] = REAL(A)[i]; */

  /* 解の受け渡し */
  /* j = j + N_var; */
  /* REAL(test_ret)[j] = (double) n_ite; */
  for (i = 0; i < N_clust_comp; i++)
	REAL(ans_F)[i] = REAL(F)[i];
  for (i = 0; i < N_var_comp; i++)
	REAL(ans_A)[i] = REAL(A)[i];
  for (i = 0; i < N_clust; i++)
	REAL(ans_g)[i] = REAL(g)[i];
  for (i = 0; i < N_sub_clust; i++)
	REAL(ans_U)[i] = REAL(U)[i];
  for (i = 0; i < N_var; i++)
	REAL(ans_MU)[i] = REAL(MU)[i];
  REAL(ans_ind)[0] = (double) n_ite;
  REAL(ans_ind)[1] = loss;
  for (i = 0; i < N_ite; i++)
	REAL(ans_record)[i] = REAL(PLL_RECORD)[i];

  SET_VECTOR_ELT(ans, 0, ans_F);
  SET_VECTOR_ELT(ans, 1, ans_A);
  SET_VECTOR_ELT(ans, 2, ans_g);
  SET_VECTOR_ELT(ans, 3, ans_U);
  SET_VECTOR_ELT(ans, 4, ans_MU);
  SET_VECTOR_ELT(ans, 5, ans_ind);
  SET_VECTOR_ELT(ans, 6, ans_record);


  // メモリ開放
  UNPROTECT(8);

  free(A_target);
  free(A_target_inv);
  free(A_target_var);
  //  free(D_mat);
  free(LogDens);
  free(F_current);
  free(F_s);
  free(F_target);
  free(G);
  free(I_comp);
  free(Mat_sv);
  free(Theta);
  free(U_svd);
  free(vec_1_clust);
  free(Vt_svd);
  free(W);
  free(work);
  free(Z);
  free(Z_s);

  return(ans);
}


/* 罰則付き対数尤度を求める */
double PenLogLikelihood(double *Q, double *F, double *A, double *mu, double *g, int N_sub, int N_var, int N_comp, int N_clust, double lambda)
{
  // 基本的な変数の定義
  double  blas_one = 1, blas_zero = 0, blas_minus_one = -1;
  int     blas_one_int = 1;
  int     i, j, k;
  double *LogDens;
  double  lossfunc;
  double  mean_density;
  int     N_clust_var = N_clust * N_var;
  double  N_sub_double = (double) N_sub;
  double  N_sub_clust = N_sub * N_clust;
  double  N_sub_clust_double = (double) N_sub_clust;
  double  prod;
  double  penll;
  double  sum;
  double  term1, term2;

  // malloc によるメモリ割り当て
  LogDens     = (double *) malloc(sizeof(double) * N_sub_clust);

  //条件付き密度の計算
  LogDensity(Q, F, A, mu, N_sub, N_var, N_comp, N_clust, LogDens);

  //Thetaの平均値を求める
  mean_density = 0;
  for (i = 1; i < N_sub_clust; i++)
  	mean_density = mean_density + LogDens[i];
  mean_density = mean_density / N_sub_clust_double;


  /***** 第1項の計算 *****/
  term1 = 0.0;
  for (i = 0; i < N_sub; i++) {
	sum = 0.0;
	for (k = 0; k < N_clust; k++) {
	  sum = sum + exp(log(g[k]) + LogDens[i + N_sub * k] - mean_density);
	}
	term1 = term1 + log(sum) + mean_density;
  }

  /***** 第2項の計算 *****/
  sum = 0;
  for (i = 0; i < (N_var * N_comp); i++)
	sum = sum + fabs(A[i]);
  term2 = sum * N_sub_double * lambda;

  /***** 罰則付き対数尤度の計算 *****/
  penll = term1 - term2;

  free(LogDens);

  return(penll);
}


/* k が与えられたもとでの Yn の条件付き密度 */
void LogDensity(double *Q, double *F, double *A, double *MU, int N_sub, int N_var, int N_comp, int N_clust, double *LogDens)
{
  // 基本的な変数の定義
  double  const_temp;
  double  density;
  int     i, j, k, l;
  int     N_clust_var = N_clust * N_var;
  double  temp;
  double *theta_mat;

  theta_mat = (double *) malloc(sizeof(double) * N_clust_var);


  for (i = 0; i < N_clust_var; i++)
	theta_mat[i] = 0.0;

  for (i = 0; i < N_sub; i++) {
	for (k = 0; k < N_clust; k++) {
	  density = 0;

	  for (j = 0; j < N_var; j++) {
		temp = 0;
		for (l = 0; l < N_comp; l++)
		  temp = temp + F[k + l * N_clust] * A[j + l * N_var];

		const_temp = 1 + exp(-Q[i + j * N_sub] * (MU[j] + temp));
		if (const_temp > 1.0E+20) { //ものすごく大きい場合
		  theta_mat[k + N_clust * j] = 1.0E-20;
		  //		  theta_mat[k + N_clust * j] = 1.0 - 1.0E-20;
		} else if (fabs(1.0 - const_temp) < 1.0E-20) { //ものすごく小さい場合
		  theta_mat[k + N_clust * j] = 1.0 - 1.0E-20;
		} else {
		  theta_mat[k + N_clust * j] = 1 / const_temp;
		}
		density = density + log(theta_mat[k + N_clust * j]);
	  }

	  LogDens[i + k * N_sub] = density;
	}
  }
}


/* Loss Function を計算する */
double F_LossCBIRD(double *F, double *A, double *U, double *mu, double *X, double lambda, int N_sub, int N_var, int N_comp, int N_clust)
{
  // 基本的な変数の定義
  double  blas_one = 1, blas_zero = 0, blas_minus_one = -1;
  int     blas_one_int = 1;
  int     i, j, k;
  double  lossfunc;
  int     N_clust_var = N_clust * N_var;
  double  N_sub_double = (double) N_sub;
  double  sum;
  double  term1, term2;
  double *Theta;
  double  vec_1_clust[N_clust];

  // malloc によるメモリ割り当て
  Theta       = (double *) malloc(sizeof(double) * N_clust_var);

  for (i = 0; i < N_clust; i++)
	vec_1_clust[i] = 1.0;

  /***** Theta の計算 *****/
  for (i = 0; i < N_clust_var; i++)
	Theta[i] = 0.0;
  dger_(&N_clust, &N_var, &blas_one, vec_1_clust, &blas_one_int, mu, &blas_one_int, Theta, &N_clust); // vec.1.clust %*% t(mu)
  dgemm_("N", "T", &N_clust, &N_var, &N_comp, &blas_one, F, &N_clust, A, &N_var, &blas_one, Theta, &N_clust); // Theta <- F %*% t(A) + vec.1.clust %*% t(mu)

  /***** 第1項 *****/
  term1 = 0.0;
  for (i = 0; i < N_sub; i++) {
	for (k = 0; k < N_clust; k++) {
	  sum = 0.0;
	  for (j = 0; j < N_var; j++)
		sum = sum + pow(Theta[k + N_clust * j] - X[i + N_sub * k + N_sub * N_clust * j], 2);
	  term1 = term1 + sum * U[i + N_sub * k];
	}
  }
  term1 = term1 / 8;

  /***** 第2項 *****/
  sum = 0.0;
  for (i = 0; i < (N_var * N_comp); i++)
	sum = sum + fabs(A[i]);
  term2 = sum * N_sub_double * lambda;

  /***** Loss function *****/
  lossfunc = term1 + term2;

  // メモリ開放
  free(Theta);

  return(lossfunc);
}


/* Loss function の F_ind における Gradient を計算する */
void F_Grad_ind(double *F, int row1, int col1, double *A, int row2, int col2, double *mu, int length1, double *X_s, int row3, int col3, double *G)
{
  // 基本的な変数の定義
  int N_sub = row1, N_var = row2, N_comp = col1;

  // 変数の定義
  double  blas_one = 1, blas_zero = 0, blas_minus_one = -1, blas_quarter = 0.25;
  int     i, j;
  double *I_comp;
  double *Mat_sc;
  double *Mat_sv;
  int     N_comp_comp = N_comp * N_comp;
  int     N_sub_var = N_sub * N_var;
  int     N_sub_comp = N_sub * N_comp;

  // malloc によるメモリ割り当て
  I_comp = (double *) malloc(sizeof(double) * N_comp_comp);
  Mat_sc = (double *) malloc(sizeof(double) * N_sub_comp);
  Mat_sv = (double *) malloc(sizeof(double) * N_sub_var);

  /* 単位行列I_varを求めておく */
  F_Identity(N_comp, I_comp);

  /* Gradient の計算 */
  dgemm_("N", "T", &N_sub,  &N_var, &N_comp,     &blas_one,      F, &N_sub,      A,  &N_var,      &blas_zero, Mat_sv, &N_sub); // F %*% t(A)
  dgemm_("N", "N", &N_sub, &N_comp,  &N_var, &blas_quarter, Mat_sv, &N_sub,      A,  &N_var,      &blas_zero, Mat_sc, &N_sub); // 1/4 * F %*% t(A) %*% A
  dgemm_("N", "N", &N_sub, &N_comp,  &N_var, &blas_quarter,    X_s, &N_sub,      A,  &N_var,      &blas_zero,      G, &N_sub); // 1/4 * Z.s %*% A
  dgemm_("N", "N", &N_sub, &N_comp, &N_comp,     &blas_one, Mat_sc, &N_sub, I_comp, &N_comp, &blas_minus_one,      G, &N_sub); // G <- 1/4 * F %*% t(A) %*% A - 1/4 * Z.s %*% A

  // メモリの開放
  free(I_comp);
  free(Mat_sc);
  free(Mat_sv);
}


/* Loss function の F における Gradient を計算する */
void F_Grad(double *F, double *A, double *mu, double *X_s, double *G, int N_sub, int N_var, int N_comp, int N_clust, double *N_k)
{
  // 変数の定義
  double  blas_one = 1.0, blas_zero = 0.0, blas_minus_one = 1.0;
  int     i, k, l;
  double *Mat_cv;
  int     N_clust_comp = N_clust * N_comp;
  int     N_clust_var = N_clust * N_var;

  // malloc によるメモリ割り当て
  Mat_cv = (double *) malloc(sizeof(double) * N_clust_var);


  dgemm_("N", "T", &N_clust, &N_var, &N_comp, &blas_one, F, &N_clust, A, &N_var, &blas_zero, Mat_cv, &N_clust); // F %*% t(A)
  for (i = 0; i < N_clust_var; i++)
	Mat_cv[i] = Mat_cv[i] - X_s[i]; // F %*% t(A) - X.F
  dgemm_("N", "N", &N_clust, &N_comp, &N_var, &blas_one, Mat_cv, &N_clust, A, &N_var, &blas_zero, G, &N_clust); // (F %*% t(A) - X.f) %*% A
  for (k = 0; k < N_clust; k++)
	for (l = 0; l < N_comp; l++)
	  G[k + N_clust * l] = N_k[k] / 4.0 * G[k + N_clust * l];

  // メモリの開放
  free(Mat_cv);
}


/* 単位行列 I_N を計算する */
void F_Identity(int N, double *I)
{
  int i, j;
  for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	  if (i == j) {
		I[j + i * N] = 1;
	  } else {
		I[j + i * N] = 0;
	  }
	}
  }
}

//double F_LossSLPCA(double *F, int row1, int col1, double *A, int row2, int col2, double *mu, int length1, double *X, int row3, int col3, double *D_mat, int row4, int col4, double lambda)
double F_LossSLPCA(double *F, int row1, int col1, double *A, int row2, int col2, double *mu, int length1, double *X, int row3, int col3)
{
  double blas_one = 1, blas_zero = 0, blas_minus_one = -1;
  int    blas_one_int = 1;
  int    i, j;
  double *I_var;
  double lossfunc;
  double *Mat_sv;
  double *Mat_vv;
  int    N_sub = row1, N_var = row2, N_comp = col1;
  double N_sub_double = (double) N_sub;
  double sum;
  double *temp_X;
  double term1, term2;
  double *vec_1_sub;

  I_var     = (double *) malloc(sizeof(double) * N_var * N_var);
  Mat_sv    = (double *) malloc(sizeof(double) * N_sub * N_var);
  Mat_vv    = (double *) malloc(sizeof(double) * N_var * N_var);
  temp_X    = (double *) malloc(sizeof(double) * N_sub * N_var);
  vec_1_sub = (double *) malloc(sizeof(double) * N_sub);


  /* 単位行列I_varを求めておく */
  F_Identity2(N_var, I_var);

  // vec.1.sub を求めておく
  for (i = 0; i < N_sub; i++)
	vec_1_sub[i] = 1.0;

  for (i = 0; i < (N_sub * N_var); i++)
	temp_X[i] = X[i];

  for (i = 0; i < (N_sub * N_var); i++)
	Mat_sv[i] = 0.0;
  dger_(&N_sub, &N_var, &blas_one, vec_1_sub, &blas_one_int, mu, &blas_one_int, Mat_sv, &N_sub); // vec.1.sub %*% t(mu)
  dgemm_("N", "T", &N_sub, &N_var, &N_comp, &blas_minus_one, F, &N_sub, A, &N_var, &blas_one, temp_X, &N_sub); // X - F %*% t(A)
  dgemm_("N", "N", &N_sub, &N_var, &N_var, &blas_one, temp_X, &N_sub, I_var, &N_var, &blas_minus_one, Mat_sv, &N_sub); // temp <- X - F %*% t(A) - vec.1.sub %*% t(mu)

  /* 第1項の計算 */
  dgemm_("T", "N", &N_var, &N_var, &N_sub, &blas_one, Mat_sv, &N_sub, Mat_sv, &N_sub, &blas_zero, Mat_vv, &N_var); // t(temp) %*% temp
  sum = 0;
  for (i = 0; i < N_var; i++) {
	for (j = 0; j < N_var; j++) {
	  if (i == j) {
		sum = sum + Mat_vv[j + i * N_var];
	  }
	}
  }
  term1 = sum / 8.0;

  lossfunc = term1;

  free(I_var);
  free(Mat_sv);
  free(Mat_vv);
  free(temp_X);
  free(vec_1_sub);

  return(lossfunc);
}


/* 単位行列 I_N を計算する */
void F_Identity2(int N, double *I)
{
  int i, j;
  for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	  if (i == j) {
		I[j + i * N] = 1;
	  } else {
		I[j + i * N] = 0;
	  }
	}
  }
}
