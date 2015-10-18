#include "math.h"
#include <functional>

// test data
//static const int N = 20;
//static const int N1 = 10;
//static const int q = 3;
//static const int qq = 5;
//static const double hx = 1.0 / N1;
//static const double hy = 1.0 / N;
//static const int m = 2;
//--------------------------------------

// real data
static const int N = 1200;
static const int N1 = 800;
static const int q = 101;
static const int qq = 20;
static const double hx = 1.0 / 100;
static const double hy = 1.0 / 200;
// real data * 2
//static const int N = 2400;
//static const int N1 = 1600;
//static const int q = 202;
//static const int qq = 40;
//static const double hx = 1.0 / 200;
//static const double hy = 1.0 / 400;
//static const int m = 25000;
static const int m = 1;
//--------------------------------------
static const int M1 = N1 + 1;
static const int M = N + 1;
static const int M2 = M1 * M;
static const double tau = 0.0005;
static const double tg = 2;
static const int cntr = N / 2;
static const int w = q;
//M1 - количество узлов по оси х
//M - количество узлов по оси y
//(2*q-1) - количество узлов в основании клина
//(qq,cntr) - номер узла вершины клина
//tg = hx/hy
//m - количество шагов по времени

static const double Gamma=1.4;
static const double Re=10000;
static const double Mah=4;
static const double Pr=0.72;
static const double omega=0.8;
static const double epsilon=0.0000000001;
static double **A;
static double *D;
static double *f;
static double *Sigma_k;
static double *Sigma_k1;
static double *u_k;
static double *u_k1;
static double *v_k;
static double *v_k1;
static double *B;
static double *u2;
static double *v2;
static double *e_k;
static double *e_k1;
static double *e_k_mu;
static double *e2;
static double *T;
static double *Sigma_kk;
static double *u_kk;
static double *v_kk;
static double *e_kk;
static double *SigmaX_k;
static double *uX_k;
static double *vY_k;
static double *eR_k;

// В массивах с _k1 хранятся значения функций на d-ом шаге по времени
// В массивах с _k хранятся значения функций с предыдущей итерации по нелинейности
// В массивах с _kk хранятся значения функций c (d-1) шага по времени
// Массивы с "2" использутся в итерациях метода Зейделя
// В массивах с X_k хранятся значения функций, вычисленных методом траекторий

int s_m, s_e, s_itr, s_end, itr = 5, itr_tr = itr;


inline double Mu(double e_k)
{
	return pow(Gamma * (Gamma - 1) * Mah * Mah * e_k * e_k, omega);
}

inline double P(double sigma_k, double e_k)
{
	return (Gamma - 1) * sigma_k * sigma_k * e_k * e_k;
}

inline void init_arrays(const int array_element_count, const int param_array_element_count)
{
	int double_size_array = 2 * array_element_count;
	A = new double*[double_size_array];
	for (int i = 0; i < double_size_array; ++i)
	{
		A[i] = new double[param_array_element_count];
	}
	for (int i = 0; i < double_size_array; ++i)
	{
		std::fill_n(A[i], param_array_element_count, 0.);
	}
	B = new double[double_size_array];
	D = new double[double_size_array];
	f = new double[double_size_array];
	Sigma_k = new double[array_element_count];
	u_k = new double[array_element_count];
	v_k = new double[array_element_count];
	Sigma_k1 = new double[array_element_count];
	u_k1 = new double[array_element_count];
	v_k1 = new double[array_element_count];
	u2 = new double[array_element_count];
	v2 = new double[array_element_count];
	e_k = new double[array_element_count];
	e_k_mu = new double[array_element_count];
	e_k1 = new double[array_element_count];
	e2 = new double[array_element_count];
	T = new double[array_element_count];
	Sigma_kk = new double[array_element_count];
	u_kk = new double[array_element_count];
	v_kk = new double[array_element_count];
	e_kk = new double[array_element_count];
	SigmaX_k = new double[array_element_count];
	uX_k = new double[array_element_count];
	vY_k = new double[array_element_count];
	eR_k = new double[array_element_count];
	std::fill_n(B, double_size_array, 0.);
	std::fill_n(D, double_size_array, 0.);
	std::fill_n(f, double_size_array, 0.);
	std::fill_n(Sigma_k, array_element_count, 0.);
	std::fill_n(u_k, array_element_count, 0.);
	std::fill_n(v_k, array_element_count, 0.);
	std::fill_n(Sigma_k1, array_element_count, 0.);
	std::fill_n(u_k1, array_element_count, 0.);
	std::fill_n(v_k1, array_element_count, 0.);
	std::fill_n(u2, array_element_count, 0.);
	std::fill_n(v2, array_element_count, 0.);
	std::fill_n(e_k, array_element_count, 0.);
	std::fill_n(e_k1, array_element_count, 0.);
	std::fill_n(e_k_mu, array_element_count, 0.);
	std::fill_n(e2, array_element_count, 0.);
	std::fill_n(T, array_element_count, 0.);
	std::fill_n(Sigma_kk, array_element_count, 0.);
	std::fill_n(u_kk, array_element_count, 0.);
	std::fill_n(v_kk, array_element_count, 0.);
	std::fill_n(e_kk, array_element_count, 0.);
	std::fill_n(SigmaX_k, array_element_count, 0.);
	std::fill_n(uX_k, array_element_count, 0.);
	std::fill_n(vY_k, array_element_count, 0.);
	std::fill_n(eR_k, array_element_count, 0.);
}