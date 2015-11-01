#include <stdio.h>
#include <math.h>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 1
#define omp_get_num_threads() 1
#endif

#ifdef __GNUC__
#define __pure
/*почему то на GCC не работает и ломает результаты расчета*/
//#define __pure __attribute__((const))
// for memory alignment
#define _ALIGN(N)  __attribute__((aligned(N))) 
#elif __INTEL_COMPILER
#define __pure __declspec(const)
// for memory alignment
#define _ALIGN(N)  __declspec(align(N))
#elif __NVCC__
#define __pure
#else
#define __pure 
#endif

//need add restrict compiler independent
//__restrict
typedef double float_type;
typedef const float_type* __restrict const cnst_arr_t;
typedef float_type* __restrict const cnst_ptr_arr_t;
typedef float_type** __restrict const cnst_ptr_2d_arr_t;

//C_M1 - количество узлов по оси х
//C_M - количество узлов по оси y
//(2*C_q-1) - количество узлов в основании клина
//(C_qq,C_cntr) - номер узла вершины клина
//C_tg = C_hx/C_hy

// test case
//static const int C_N = 20;
//static const int C_N1 = 10;
//static const int C_q = 3;
//static const int C_qq = 5;
//static const float_type C_hx = 1.0 / C_N1;
//static const float_type C_hy = 1.0 / C_N;
//static const int time_steps_nbr = 2; // time_steps_nbr - количество шагов по времени
//-----------------------
// real test 
static const int C_N = 1200;
static const int C_N1 = 800;
static const int C_q = 101;
static const int C_qq = 20;
static const float_type C_hx = 1.0 / 100;
static const float_type C_hy = 1.0 / 200;

// real test  * 2
//static const int C_N = 2400;
//static const int C_N1 = 1600;
//static const int C_q = 202;
//static const int C_qq = 40;
//static const float_type C_hx = 1.0 / 200;
//static const float_type C_hy = 1.0 / 400;

//static const int time_steps_nbr = 25000; // time_steps_nbr - количество шагов по времени
static const int time_steps_nbr = 1; // time_steps_nbr - количество шагов по времени
//-----------------------
static const int C_M = C_N + 1;
static const int C_M1 = C_N1 + 1;
static const int C_M2 = C_M1 * C_M;
static const int C_w = C_q;
static const int C_cntr = C_N / 2;
static const int C_br = (C_N1 - 1) * (C_N - 1);
static const float_type C_tau = 0.0005;
static const float_type C_tg = 2;
static const float_type C_Re = 10000;
static const float_type C_PrRe = 0.72 * C_Re; // Pr * C_Re
static const float_type C_Mah2 = 16; // Mah * Mah
static const float_type C_epsilon = 0.0000000001;
static const float_type C_gamma = 1.4;
static const float_type C_gamma_Mah2 = C_gamma * (C_gamma - 1) * C_Mah2;

// В массивах с _k хранятся значения функций с предыдущей итерации по нелинейности
// В массивах с _kk хранятся значения функций c (d-1) шага по времени
// В массивах с _k1 хранятся значения функций на d-ом шаге по времени
// Массивы с "2" использутся в итерациях метода Якоби

static float_type** A;
static float_type* f;

static float_type* sigma_k;
static float_type* sigma_k1;
static float_type* sigma_kk;

static float_type* u_k;
static float_type* u_k1;
static float_type* u_kk;
static float_type* u2;

static float_type* v_k;
static float_type* v_k1;
static float_type* v_kk;
static float_type* v2;

static float_type* e_k;
static float_type* e_k1;
static float_type* e_kk;
static float_type* e2;
static float_type* e_k_mu;

inline void print_new_line(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure)
{
	fprintf(out, "\n\n");
	fprintf(density, "\n\n");
	fprintf(velocity, "\n\n");
	fprintf(temperature, "\n\n");
	fprintf(pressure, "\n\n");
}

inline void flush_file(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	fflush(out);
	fflush(density);
	fflush(velocity);
	fflush(temperature);
	fflush(pressure);
	fflush(out_itr);
}

inline __pure float_type Mu(float_type value)
{
	const float_type omega = 0.8;
	return pow(C_gamma_Mah2 * value * value, omega);
}

#pragma omp declare simd
inline __pure float_type P(float_type gamma_value, float_type sigma_k_value, float_type e_k_value)
{
	return (gamma_value - 1) * sigma_k_value * sigma_k_value * e_k_value * e_k_value;
}

inline void print_to_file(const float_type gamma, int s_m, int s_e,
                          int current_ts, int s_itr, int s_end,
						  const float_type tau, const float_type hx,
						  const float_type hy,
                          const int m, const int m1, const int n,
						  const float_type mah2,
                          FILE* out, FILE* density, FILE* density_new, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{	
	if (current_ts / 1. == 1
		|| current_ts / 100. == 1
		|| current_ts / 500. == 1
		|| current_ts / 1000. == 1
		|| current_ts / 2000. == 1
		|| current_ts / 3000. == 1
		|| current_ts / 3500. == 1
		|| current_ts / 4000. == 1
		|| current_ts / 4500. == 1
		|| current_ts / 5000. == 1
		|| current_ts / 6000. == 1
		|| current_ts / 6500. == 1
		|| current_ts / 7000. == 1
		|| current_ts / 7500. == 1
		|| current_ts / 8000. == 1
		|| current_ts / 8500. == 1
		|| current_ts / 9000. == 1
		|| current_ts / 9500. == 1
		|| current_ts / 10000. == 1
		|| current_ts / 11000. == 1
		|| current_ts / 12000. == 1
		|| current_ts / 13000. == 1
		|| current_ts / 14000. == 1
		|| current_ts / 15000. == 1
		|| current_ts / 16000. == 1
		|| current_ts / 17000. == 1
		|| current_ts / 18000. == 1
		|| current_ts / 19000. == 1
		|| current_ts / 20000. == 1
		|| current_ts / 23000. == 1
		|| current_ts / 25000. == 1
		|| current_ts / 28000. == 1
		|| current_ts / 30000. == 1
		|| current_ts / 33000. == 1
		|| current_ts / 35000. == 1
		|| current_ts / 38000. == 1
		|| current_ts / 40000. == 1
		|| current_ts / 43000. == 1
		|| current_ts / 45000. == 1
		|| current_ts / 48000. == 1
		|| current_ts / 50000. == 1
		|| current_ts / 53000. == 1
		|| current_ts / 55000. == 1
		|| current_ts / 58000. == 1
		|| current_ts / 60000. == 1
		|| current_ts / 63000. == 1
		|| current_ts / 65000. == 1
		|| current_ts / 68000. == 1
		|| current_ts / 70000. == 1
		|| current_ts / 73000. == 1
		|| current_ts / 75000. == 1
		|| current_ts / 78000. == 1
		|| current_ts / 80000. == 1
		|| current_ts / 83000. == 1
		|| current_ts / 85000. == 1
		|| current_ts / 88000. == 1
		|| current_ts / 90000. == 1
		|| current_ts / 93000. == 1
		|| current_ts / 95000. == 1
		|| current_ts / 98000. == 1
		|| current_ts / 100000. == 1
		|| current_ts / 103000. == 1
		|| current_ts / 105000. == 1
		|| current_ts / 108000. == 1
		|| current_ts / 110000. == 1
		|| current_ts / 113000. == 1
		|| current_ts / 115000. == 1
		|| current_ts / 118000. == 1
		|| current_ts / 120000. == 1
		|| current_ts / 123000. == 1
		|| current_ts / 125000. == 1
		|| current_ts / 128000. == 1
		|| current_ts / 130000. == 1
		|| current_ts / 133000. == 1
		|| current_ts / 135000. == 1
		|| current_ts / 138000. == 1
		|| current_ts / 140000. == 1
		|| current_ts / 143000. == 1
		|| current_ts / 145000. == 1
		|| current_ts / 148000. == 1
		|| current_ts / 150000. == 1
		|| current_ts / 153000. == 1
		|| current_ts / 155000. == 1
		|| current_ts / 158000. == 1
		|| current_ts / 160000. == 1
		|| current_ts / 163000. == 1
		|| current_ts / 165000. == 1
		|| current_ts / 168000. == 1
		|| current_ts / 170000. == 1
		|| current_ts / 173000. == 1
		|| current_ts / 175000. == 1
		|| current_ts / 178000. == 1
		|| current_ts / 180000. == 1
		|| current_ts / 183000. == 1
		|| current_ts / 185000. == 1
		|| current_ts / 188000. == 1
		|| current_ts / 190000. == 1
		|| current_ts / 193000. == 1
		|| current_ts / 195000. == 1
		|| current_ts / 198000. == 1
		|| current_ts / 200000. == 1)
	{
		float_type ts_tau = current_ts * tau;

		fprintf(out, "\t\t d = %i\t d*C_tau = %.5f\n\n", current_ts, ts_tau);
		fprintf(out, "C_q = %i\t C_w = %i\n\n", C_q, C_w);
		fprintf(out_itr, "\t\t d = %i\t d*C_tau = %.5f\n\n", current_ts, ts_tau);
		if (s_end == 0)
		{
			fprintf(out_itr, "s_m = %i\t s_e = %i\t s_itr = %i\n\n", s_m, s_e, s_itr - 1);
			fprintf(out, "s_m = %i\t s_e = %i\t s_itr = %i\n\n", s_m, s_e, s_itr - 1);
		}
		else
		{
			fprintf(out_itr, "s_m = %i\t s_e = %i\t s_end = %i\n\n", s_m, s_e, s_end);
			fprintf(out, "s_m = %i\t s_e = %i\t s_end = %i\n\n", s_m, s_e, s_end);
		}

		fprintf(out, "\t\t rho\t\t u\t\t v\t\t e\t\t T\t\t P\t\t Mu\n\n");
		fprintf(density, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", n, ts_tau, m1, m);
		fprintf(velocity, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", n, ts_tau, m1, m);
		fprintf(temperature, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", n, ts_tau, m1, m);
		fprintf(pressure, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", n, ts_tau, m1, m);

		for (int i = 0; i < m1; i++)		
		{
			for (int j = 0; j < m; j++)
			{
				int a = i * m + j;
				float_type ihx = i * hx;
				float_type jhy = j * hy;
				fprintf(density_new, "%.3f\t %.4f\t %.15e\n", ihx, jhy, sigma_k1[a] * sigma_k1[a]);
				fprintf(density, "%.3f\t %.4f\t %.10f\n", ihx, jhy, sigma_k1[a] * sigma_k1[a]);
				fprintf(velocity, "%.3f\t %.4f\t%.10f\t %.10f\n", ihx, jhy, u_k1[a], v_k1[a]);
				fprintf(temperature, "%.3f\t %.4f\t %.10f\n", ihx, jhy, e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * mah2));
				fprintf(pressure, "%.3f\t %.4f\t %.10f\n", ihx, jhy, sigma_k1[a] * sigma_k1[a] * e_k1[a] * e_k1[a] * (gamma - 1));
				if (current_ts == 1)
				{
					fprintf(out, "i=%i j=%i\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\n", i, j,
					        sigma_k1[a] * sigma_k1[a], u_k1[a], v_k1[a], e_k1[a] * e_k1[a], e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * mah2), P(gamma, sigma_k1[a], e_k1[a]), Mu(e_k1[a]));
				}
			}
		}		
		print_new_line(out, density, velocity, temperature, pressure);
		flush_file(out, density, velocity, temperature, pressure, out_itr);
	}
}

inline void print_file_header(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr, const float_type tau, const float_type hx, const float_type hy)
{
	fprintf(out, "Cone_2D\n\n");
	fprintf(out, "C_N = %i\t C_hx = %.5f\t C_hy = %.5f\t C_tau = %.5f\n\n", C_N, hx, hy, tau);
	fprintf(out_itr, "Cone_2D\n\n");
	fprintf(out_itr, "C_N = %i\t C_hx = %.5f\t C_hy = %.5f\t C_tau = %.5f\n\n", C_N, hx, hy, tau);
	fprintf(density, "TITLE=\"density\"\n\nVARIABLES=\"x\",\"y\",\"Ro\"\n\n");
	fprintf(velocity, "TITLE=\"velocity\"\n\nVARIABLES=\"x\",\"y\",\"u\",\"v\"\n\n");
	fprintf(temperature, "TITLE=\"temperature\"\n\nVARIABLES=\"x\",\"y\",\"T\"\n\n");
	fprintf(pressure, "TITLE=\"pressure\"\n\nVARIABLES=\"x\",\"y\",\"P\"\n\n");
}

// Initial boundary conditions with t = 0
inline void set_initial_boundary_conditions_parallel()
{
	const float_type sqrt1 = sqrt(1 / (C_gamma * (C_gamma - 1) * C_Mah2));

	for (int i = 0; i < C_M1; i++)
	{
		for (int j = 0; j < C_M; j++)
		{
			if (i < C_qq || i >= C_qq + C_w - 1)
			{
				sigma_k1[i * C_M + j] = 1;
				e_k1[i * C_M + j] = sqrt1;
				u_k1[i * C_M + j] = i == 0 ? 1 : 0;
				e2[i * C_M + j] = e_k1[i * C_M + j];
				u2[i * C_M + j] = u_k1[i * C_M + j];
			}
		}
	}

	for (int i = C_qq; i < C_qq + C_w - 1; i++)
	{
		for (int j = C_cntr + i - C_qq; j < C_M; j++)
		{
			sigma_k1[i * C_M + j] = 1;
			e_k1[i * C_M + j] = sqrt1;
			u_k1[i * C_M + j] = 0;
			e2[i * C_M + j] = e_k1[i * C_M + j];
			u2[i * C_M + j] = u_k1[i * C_M + j];
		}
		for (int j = C_cntr - i + C_qq; j >= 0; j--)
		{
			sigma_k1[i * C_M + j] = 1;
			e_k1[i * C_M + j] = sqrt1;
			u_k1[i * C_M + j] = 0;
			e2[i * C_M + j] = e_k1[i * C_M + j];
			u2[i * C_M + j] = u_k1[i * C_M + j];
		}
	}
}

#pragma omp declare simd
inline __pure float_type trajectory(const float_type tau_d, const float_type hx_d, const float_type hy_d, int i, int j, cnst_arr_t arr, const float_type u_k_value, const float_type v_k_value, const int m_i)
{
	int idx = i * m_i + j;
	int idx2 = i * m_i + (j - 1);
	int idx3 = (i - 1) * m_i + j;
	int idx4 = (i + 1) * m_i + j;
	int idx5 = i * m_i + (j + 1);

	if (u_k_value == 0. && v_k_value == 0.)
	{
		return 0;
	}

	if (u_k_value > 0. && v_k_value > 0.)
	{
		return u_k_value * tau_d * ((arr[idx3] - arr[idx]) / ((i - 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx2] - arr[idx]) / ((j - 1) * hy_d - j * hy_d));
	}

	if (u_k_value < 0. && v_k_value > 0.)
	{
		return u_k_value * tau_d * ((arr[idx4] - arr[idx]) / ((i + 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx2] - arr[idx]) / ((j - 1) * hy_d - j * hy_d));
	}

	if (u_k_value < 0. && v_k_value < 0.)
	{
		return u_k_value * tau_d * ((arr[idx4] - arr[idx]) / ((i + 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx5] - arr[idx]) / ((j + 1) * hy_d - j * hy_d));
	}
	if (u_k_value > 0. && v_k_value < 0.)
	{
		return u_k_value * tau_d * ((arr[idx3] - arr[idx]) / ((i - 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx5] - arr[idx]) / ((j + 1) * hy_d - j * hy_d));
	}
	if (u_k_value > 0. && v_k_value == 0.)
	{
		return u_k_value * tau_d * ((arr[idx3] - arr[idx]) / ((i - 1) * hx_d - i * hx_d));
	}
	if (u_k_value == 0. && v_k_value > 0.)
	{
		return v_k_value * tau_d * ((arr[idx2] - arr[idx]) / ((j - 1) * hy_d - j * hy_d));
	}
	if (u_k_value < 0. && v_k_value == 0.)
	{
		return u_k_value * tau_d * ((arr[idx4] - arr[idx]) / ((i + 1) * hx_d - i * hx_d));
	}
	if (u_k_value == 0. && v_k_value < 0.)
	{
		return v_k_value * tau_d * ((arr[idx5] - arr[idx]) / ((j + 1) * hy_d - j * hy_d));
	}
	return 0;
}

inline void continuity(cnst_arr_t sigma_kk_arr, cnst_ptr_arr_t sigma_k1_arr, cnst_arr_t u_k_arr, cnst_arr_t v_k_arr)
{
	//Для внутренних узлов
//#pragma omp parallel
	{
//#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) / C_tau / (1 / C_tau + (u_k_arr[(i + 1) * C_M + j] - u_k_arr[(i - 1) * C_M + j]) / (4 * C_hx)
					+ (v_k_arr[i * C_M + j + 1] - v_k_arr[i * C_M + j - 1]) / (4 * C_hy));
			}
		}
//#pragma omp for nowait
		for (int j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
		{
			//Для Г5. l = C_w-1; m = 1,...,C_q-2;
			int i = C_qq + C_w - 1;
			sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) / (2 * C_tau) / (1 / (2 * C_tau) + (u_k_arr[(i + 1) * C_M + j] - u_k_arr[i * C_M + j]) / (4 * C_hx)
				+ (v_k_arr[i * C_M + j + 1] - v_k_arr[i * C_M + j - 1]) / (8 * C_hy));
		}
//#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//Для Г6.
			int j = C_cntr + i - C_qq;
			sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) * (1 / (4 * C_tau) + 1 / (4 * C_tau)) / (1 / (4 * C_tau) + 1 / (4 * C_tau) + (u_k_arr[i * C_M + j] - u_k_arr[(i - 1) * C_M + j]) / (8 * C_hx)
				- u_k_arr[(i - 1) * C_M + j] / (16 * C_hx) + (v_k_arr[i * C_M + j + 1] - v_k_arr[i * C_M + j]) / (8 * C_hy) + v_k_arr[i * C_M + j + 1] / (16 * C_hy));
			//Для Г7.
			j = C_cntr - i + C_qq;
			sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) * (1 / (4 * C_tau) + 1 / (4 * C_tau)) / (1 / (4 * C_tau) + 1 / (4 * C_tau) + (u_k_arr[i * C_M + j] - u_k_arr[(i - 1) * C_M + j]) / (8 * C_hx)
				- u_k_arr[(i - 1) * C_M + j] / (16 * C_hx) + (v_k_arr[i * C_M + j] - v_k_arr[i * C_M + j - 1]) / (8 * C_hy) - v_k_arr[i * C_M + j - 1] / (16 * C_hy));
		}
//#pragma omp single
		{
			//Для S_w-1q-1.
			int i = C_qq + C_w - 1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			sigma_k1_arr[a] = (sigma_kk_arr[a] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[a], v_k_arr[a], C_M)) * (3 / (4 * C_tau) + 1 / (8 * C_tau)) / (3 / (4 * C_tau) + 1 / (8 * C_tau) + (2 * u_k_arr[(i + 1) * C_M + j] - u_k_arr[(i - 1) * C_M + j] - u_k_arr[a]) / (8 * C_hx)
				+ (u_k_arr[a] - u_k_arr[(i - 1) * C_M + j]) / (16 * C_hx) + (2 * v_k_arr[a + 1] - v_k_arr[a - 1] - v_k_arr[a]) / (8 * C_hy) + v_k_arr[a] / (16 * C_hy));
			//Для S.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			sigma_k1_arr[a] = (sigma_kk_arr[a] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[a], v_k_arr[a], C_M)) * (3 / (4 * C_tau) + 1 / (8 * C_tau)) / (3 / (4 * C_tau) + 1 / (8 * C_tau) + (2 * u_k_arr[(i + 1) * C_M + j] - u_k_arr[(i - 1) * C_M + j] - u_k_arr[a]) / (8 * C_hx)
				+ (u_k_arr[a] - u_k_arr[(i - 1) * C_M + j]) / (16 * C_hx) + (v_k_arr[a + 1] - 2 * v_k_arr[a - 1] + v_k_arr[a]) / (8 * C_hy) - v_k_arr[a] / (16 * C_hy));
			//Для S_qq_0
			i = C_qq;
			j = C_cntr;
			a = i * C_M + j;
			sigma_k1_arr[a] = (sigma_kk_arr[a] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[a], v_k_arr[a], C_M)) * (1 / (4 * C_tau) + 1 / (2 * C_tau)) / (1 / (4 * C_tau) + 1 / (2 * C_tau) + (u_k_arr[a] - u_k_arr[(i - 1) * C_M + j]) / (4 * C_hx)
				+ (u_k_arr[(i + 1) * C_M + j] - u_k_arr[a]) / (8 * C_hx) + (v_k_arr[a + 1] - v_k_arr[a - 1]) / (8 * C_hy) + (v_k_arr[a + 1] - v_k_arr[a - 1]) / (16 * C_hy));
		}
	} // //#pragma omp parallel
}

/*Energy*/
/*----- Функция заполняет элементы матрицы для уравнения энергии----*/
inline void nrg_calculate_common(cnst_ptr_2d_arr_t a_arr, cnst_arr_t sigma_k_arr, cnst_arr_t sigma_k1_arr, cnst_arr_t e_k_arr, cnst_arr_t e_k_mu_arr, cnst_arr_t u_k_arr, cnst_arr_t v_k_arr, 
	cnst_ptr_arr_t f_arr)
{
	const float_type c_coef1 = 2 * C_hx * C_hx * C_PrRe;
	const float_type c_coef2 = 2 * C_hy * C_hy * C_PrRe;
	const float_type c_coef3 = 4 * C_hy * C_hy * C_PrRe;
	const float_type c_coef4 = 8 * C_hy * C_hy * C_PrRe;
	const float_type c_coef5 = 2 * C_tg * C_hx * C_hy * C_PrRe;
	const float_type c_coef6 = 4 * C_hx * C_hx * C_PrRe;
	const float_type c_coef7 = 8 * C_hx * C_hx * C_PrRe;

//#pragma omp parallel 
	{
//#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;
					a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
					a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
					a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
						C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
						C_gamma / c_coef1 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[(i - 1) * C_M + j]) +
						C_gamma / c_coef2 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a + 1] + e_k_mu_arr[a - 1]);
					a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
					a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));
					f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (4 * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +
						e_k_mu_arr[a] / (6 * C_hx * C_hx * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +
						e_k_mu_arr[a] / (6 * C_hy * C_hy * C_Re * e_k_arr[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +
						e_k_mu_arr[a] / (8 * C_Re * e_k_arr[a]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
							(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
							(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
							(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
						e_k_mu_arr[a] / (12 * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
							(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
							(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
							(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
				}
			}
		}
//#pragma omp for nowait 
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
			a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1])) - C_gamma / c_coef3 * (e_k_mu_arr[a - 1] + e_k_mu_arr[a]) - C_gamma / c_coef4 * (e_k_mu_arr[a - 1] + 2 * e_k_mu_arr[a]);
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
				C_gamma / c_coef7 * (8 * e_k_mu_arr[a] + 3 * e_k_mu_arr[(i + 1) * C_M + j] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (8 * e_k_mu_arr[a] + 4 * e_k_mu_arr[a + 1] + 3 * e_k_mu_arr[a - 1]);
			a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a])) - C_gamma / c_coef6 * (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]) - C_gamma / c_coef7 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]);

			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
			a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
				C_gamma / c_coef7 * (8 * e_k_mu_arr[a] + 3 * e_k_mu_arr[(i + 1) * C_M + j] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (8 * e_k_mu_arr[a] + 3 * e_k_mu_arr[a + 1] + 4 * e_k_mu_arr[a - 1]);
			a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a])) - C_gamma / c_coef3 * (e_k_mu_arr[a] + e_k_mu_arr[a + 1]) - C_gamma / c_coef4 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a + 1]);
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a])) - C_gamma / c_coef6 * (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]) - C_gamma / c_coef7 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]);
		}
//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				int a = i * C_M + j;
				a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
				a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
				a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
					C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
					C_gamma / c_coef1 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[(i - 1) * C_M + j]) +
					C_gamma / c_coef2 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a + 1] + e_k_mu_arr[a - 1]);
				a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
				a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				int a = i * C_M + j;
				a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
				a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
				a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
					C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
					C_gamma / c_coef1 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[(i - 1) * C_M + j]) +
					C_gamma / c_coef2 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a + 1] + e_k_mu_arr[a - 1]);
				a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
				a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));
			}
		}
//#pragma omp for nowait
		for (int j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
		{
			//Для Г5. l = C_q-1; C_M = 1,...,C_q-1;	
			int i = C_qq + C_w - 1;
			int a = i * C_M + j;
			a_arr[a][1] = C_gamma / c_coef3 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / (2 * C_tau) - C_gamma / c_coef6 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - 2 * e_k_arr[(i + 1) * C_M + j]) -
				C_gamma / c_coef3 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
				C_gamma / c_coef6 * (2 * e_k_mu_arr[a] + 2 * e_k_mu_arr[(i + 1) * C_M + j]) +
				C_gamma / c_coef3 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a - 1] + e_k_mu_arr[a + 1]);
			a_arr[a][3] = -C_gamma / c_coef3 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[a]));
			f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / (2 * C_tau) - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (8 * e_k_arr[a]) * ((2 * u_k1[(i + 1) * C_M + j] - 2 * u_k1[a]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

				e_k_mu_arr[a] / (6 * C_hx * C_hx * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a])) +

				e_k_mu_arr[a] / (12 * C_hy * C_hy * C_Re * e_k_arr[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

				e_k_mu_arr[a] / (8 * C_Re * e_k_arr[a]) * (((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) * ((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) +
					((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy) * ((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy)) +

				e_k_mu_arr[a] / (12 * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		}
//#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//Для Г6. l = 1,...,C_q-1; C_M = C_q-1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef6 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu_arr[a] + e_k_mu_arr[(i - 1) * C_M + j]) - C_gamma / c_coef7 * (e_k_mu_arr[a] + 2 * e_k_mu_arr[(i - 1) * C_M + j]);
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (4 * e_k_arr[a] - 3 * e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (4 * e_k_arr[a] - 3 * e_k_arr[a + 1]) +
				C_gamma / c_coef7 * (4 * e_k_mu_arr[a] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (4 * e_k_mu_arr[a] + 4 * e_k_mu_arr[a + 1]);
			a_arr[a][3] = -C_gamma / c_coef3 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) - C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a])
				- C_gamma / c_coef3 * (e_k_mu_arr[a + 1] + e_k_mu_arr[a]) - C_gamma / c_coef4 * (2 * e_k_mu_arr[a + 1] + e_k_mu_arr[a]);
			a_arr[a][5] = -C_gamma / c_coef5 * e_k_mu_arr[a];
			a_arr[a][8] = C_gamma / c_coef5 * e_k_mu_arr[a];
			f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (8 * e_k_arr[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a] / C_hy) -
				-P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (16 * e_k_arr[a]) * (-u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a + 1] / C_hy) +

				e_k_mu_arr[a] / (24 * C_hx * C_hx * C_Re * e_k_arr[a]) * (-u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

				e_k_mu_arr[a] / (24 * C_hy * C_hy * C_Re * e_k_arr[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + v_k1[a] * v_k1[a]) +

				e_k_mu_arr[a] / (16 * C_Re * e_k_arr[a]) * ((-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy)) +

				e_k_mu_arr[a] / (24 * C_Re * e_k_arr[a]) * ((-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy));

			// Для Г7.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef6 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu_arr[a] + e_k_mu_arr[(i - 1) * C_M + j]) - C_gamma / c_coef7 * (e_k_mu_arr[a] + 2 * e_k_mu_arr[(i - 1) * C_M + j]);
			a_arr[a][1] = C_gamma / c_coef3 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) + C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1])
				- C_gamma / c_coef3 * (e_k_mu_arr[a] + e_k_mu_arr[a - 1]) - C_gamma / c_coef4 * (e_k_mu_arr[a] + 2 * e_k_mu_arr[a - 1]);
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (4 * e_k_arr[a] - 3 * e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (4 * e_k_arr[a] - 3 * e_k_arr[a - 1]) +
				C_gamma / c_coef7 * (4 * e_k_mu_arr[a] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (4 * e_k_mu_arr[a] + 4 * e_k_mu_arr[a - 1]);
			a_arr[a][6] = -C_gamma / c_coef5 * e_k_mu_arr[a];
			a_arr[a][7] = C_gamma / c_coef5 * e_k_mu_arr[a];
			f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (8 * e_k_arr[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				-P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (16 * e_k_arr[a]) * (-u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a - 1] / C_hy) +

				e_k_mu_arr[a] / (24 * C_hx * C_hx * C_Re * e_k_arr[a]) * (-u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

				e_k_mu_arr[a] / (24 * C_hy * C_hy * C_Re * e_k_arr[a]) * (3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + -v_k1[a] * -v_k1[a]) +

				e_k_mu_arr[a] / (16 * C_Re * e_k_arr[a]) * ((-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy)) +

				e_k_mu_arr[a] / (24 * C_Re * e_k_arr[a]) * ((-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy));
		}
//#pragma omp single nowait
		{
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
			a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
				C_gamma / c_coef1 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef2 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a + 1] + e_k_mu_arr[a - 1]);
			a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));

			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
			a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - C_gamma / c_coef1 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[(i + 1) * C_M + j] - e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu_arr[a] / e_k_arr[a] * (2 * e_k_arr[a] - e_k_arr[a + 1] - e_k_arr[a - 1]) +
				C_gamma / c_coef1 * (2 * e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef2 * (2 * e_k_mu_arr[a] + e_k_mu_arr[a + 1] + e_k_mu_arr[a - 1]);
			a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));

			//Для S_qq, C_N / 2 + C_q.
			j = C_cntr + i - C_qq;
			a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef6 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]) - C_gamma / c_coef7 * (2 * e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]);
			a_arr[a][1] = C_gamma / c_coef3 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (7 * e_k_arr[a] - 4 * e_k_arr[(i + 1) * C_M + j] - 3 * e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (7 * e_k_arr[a] - 4 * e_k_arr[a + 1] - 2 * e_k_arr[a - 1]) +
				C_gamma / c_coef7 * (7 * e_k_mu_arr[a] + 4 * e_k_mu_arr[(i + 1) * C_M + j] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (7 * e_k_mu_arr[a] + 4 * e_k_mu_arr[a + 1] + 2 * e_k_mu_arr[a - 1]) + C_gamma / c_coef5 * e_k_mu_arr[a];
			a_arr[a][3] = -C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));
			a_arr[a][5] = -C_gamma / c_coef5 * e_k_mu_arr[a];
			f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (8 * e_k_arr[a]) * (2 * u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + 2 * v_k1[a + 1] / C_hy - v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (16 * e_k_arr[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) +
				e_k_mu_arr[a] / (24 * C_hx * C_hx * C_Re * e_k_arr[a]) * (4 * (u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +
				e_k_mu_arr[a] / (24 * C_hy * C_hy * C_Re * e_k_arr[a]) * (4 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 2 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + v_k1[a] * v_k1[a]) +
				e_k_mu_arr[a] / (16 * C_Re * e_k_arr[a]) * (2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy)) +
				e_k_mu_arr[a] / (24 * C_Re * e_k_arr[a]) * (2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy));

			//Для S_qq, C_N/2 - C_q.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef6 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]) - C_gamma / c_coef7 * (2 * e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]);
			a_arr[a][1] = C_gamma / c_coef2 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) - (e_k_mu_arr[a - 1] + e_k_mu_arr[a]));
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (7 * e_k_arr[a] - 4 * e_k_arr[(i + 1) * C_M + j] - 3 * e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (7 * e_k_arr[a] - 2 * e_k_arr[a + 1] - 4 * e_k_arr[a - 1]) +
				C_gamma / c_coef7 * (7 * e_k_mu_arr[a] + 4 * e_k_mu_arr[(i + 1) * C_M + j] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (7 * e_k_mu_arr[a] + 2 * e_k_mu_arr[a + 1] + 4 * e_k_mu_arr[a - 1]) + C_gamma / c_coef5 * e_k_mu_arr[a];
			a_arr[a][3] = -C_gamma / c_coef3 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[a + 1]));
			a_arr[a][4] = -C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[(i + 1) * C_M + j] - e_k_arr[a]) + (e_k_mu_arr[a] + e_k_mu_arr[(i + 1) * C_M + j]));
			a_arr[a][6] = -C_gamma / c_coef5 * e_k_mu_arr[a];
			f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (8 * e_k_arr[a]) * (2 * u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx + v_k1[a + 1] / C_hy + v_k1[a] / C_hy - u_k1[(i - 1) * C_M + j] / C_hx - 2 * v_k1[a - 1] / C_hy) -
				P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (16 * e_k_arr[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy) +
				e_k_mu_arr[a] / (24 * C_hx * C_hx * C_Re * e_k_arr[a]) * (4 * (u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +
				e_k_mu_arr[a] / (24 * C_hy * C_hy * C_Re * e_k_arr[a]) * (4 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + 2 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + -v_k1[a] * -v_k1[a]) +
				e_k_mu_arr[a] / (16 * C_Re * e_k_arr[a]) * (2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy) +
					2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
				e_k_mu_arr[a] / (24 * C_Re * e_k_arr[a]) * (2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) +
					2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));

			//Для S_qq,C_N/2
			i = C_qq;
			j = C_cntr;
			a = i * C_M + j;
			a_arr[a][0] = C_gamma / c_coef1 * (e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[(i - 1) * C_M + j]) - (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[a]));
			a_arr[a][1] = C_gamma / c_coef3 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1]) + C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a] - e_k_arr[a - 1])
				- C_gamma / c_coef3 * (e_k_mu_arr[a - 1] + e_k_mu_arr[a]) - C_gamma / c_coef4 * (2 * e_k_mu_arr[a - 1] + e_k_mu_arr[a]);
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - C_gamma / c_coef7 * e_k_mu_arr[a] / e_k_arr[a] * (6 * e_k_arr[a] - 4 * e_k_arr[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (6 * e_k_arr[a] - 3 * e_k_arr[a + 1] - 3 * e_k_arr[a - 1]) +
				C_gamma / c_coef7 * (6 * e_k_mu_arr[a] + 4 * e_k_mu_arr[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (6 * e_k_mu_arr[a] + 4 * e_k_mu_arr[a + 1] + 4 * e_k_mu_arr[a - 1]) - C_gamma / (C_tg * C_hx * C_hy * C_PrRe) * e_k_mu_arr[a];
			a_arr[a][3] = -C_gamma / c_coef3 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a]) - C_gamma / c_coef4 * e_k_mu_arr[a] / e_k_arr[a] * (e_k_arr[a + 1] - e_k_arr[a])
				- C_gamma / c_coef3 * (e_k_mu_arr[a] + e_k_mu_arr[a + 1]) - C_gamma / c_coef4 * (e_k_mu_arr[a] + 2 * e_k_mu_arr[a + 1]);
			a_arr[a][7] = C_gamma / c_coef5 * e_k_mu_arr[a];
			a_arr[a][8] = C_gamma / c_coef5 * e_k_mu_arr[a];
			f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (8 * e_k_arr[a]) * (2 * u_k1[a] / C_hx - 2 * u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) -
				P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (16 * e_k_arr[a]) * (-2 * u_k1[a] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) +
				e_k_mu_arr[a] / (12 * C_hx * C_hx * C_Re * e_k_arr[a]) * (2 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j]) + -u_k1[a] * -u_k1[a]) +
				e_k_mu_arr[a] / (24 * C_hy * C_hy * C_Re * e_k_arr[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +
				e_k_mu_arr[a] / (16 * C_Re * e_k_arr[a]) * (2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
				e_k_mu_arr[a] / (24 * C_Re * e_k_arr[a]) * (2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		} // //#pragma omp single		
	} // //#pragma omp parallel

//#pragma omp parallel 
	{
//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				int a = i * C_M + j;
				f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (4 * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					e_k_mu_arr[a] / (6 * C_hx * C_hx * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

					e_k_mu_arr[a] / (6 * C_hy * C_hy * C_Re * e_k_arr[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					e_k_mu_arr[a] / (8 * C_Re * e_k_arr[a]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					e_k_mu_arr[a] / (12 * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
			for (int j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				int a = i * C_M + j;
				f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (4 * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					e_k_mu_arr[a] / (6 * C_hx * C_hx * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

					e_k_mu_arr[a] / (6 * C_hy * C_hy * C_Re * e_k_arr[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					e_k_mu_arr[a] / (8 * C_Re * e_k_arr[a]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					e_k_mu_arr[a] / (12 * C_Re * e_k_arr[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}
	}
}

//Вектор B = A*Xk1
inline void nrg_calculate_jakobi(cnst_arr_t e_k_arr, cnst_ptr_arr_t e2_arr)
{
//#pragma omp parallel 
	{
		//Для внутренних узлов
//#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;
					float_type b = A[a][0] * e_k_arr[(i - 1) * C_M + j] + A[a][1] * e_k_arr[i * C_M + (j - 1)] + A[a][2] * e_k_arr[a] +
						A[a][3] * e_k_arr[i * C_M + (j + 1)] + A[a][4] * e_k_arr[(i + 1) * C_M + j];
					e2_arr[a] = e_k_arr[a] - 1 / A[a][2] * (b - f[a]);
				}
			}
		}
//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				int a = i * C_M + j;
				float_type b = A[a][0] * e_k_arr[(i - 1) * C_M + j] + A[a][1] * e_k_arr[i * C_M + (j - 1)] + A[a][2] * e_k_arr[a] +
					A[a][3] * e_k_arr[i * C_M + (j + 1)] + A[a][4] * e_k_arr[(i + 1) * C_M + j];
				e2_arr[a] = e_k_arr[a] - 1 / A[a][2] * (b - f[a]);
			}
			for (int j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				int a = i * C_M + j;
				float_type b = A[a][0] * e_k_arr[(i - 1) * C_M + j] + A[a][1] * e_k_arr[i * C_M + (j - 1)] + A[a][2] * e_k_arr[a] +
					A[a][3] * e_k_arr[i * C_M + (j + 1)] + A[a][4] * e_k_arr[(i + 1) * C_M + j];
				e2_arr[a] = e_k_arr[a] - 1 / A[a][2] * (b - f[a]);
			}
		}
//#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//Для Г6. l = 1,...,C_q-1; C_M = C_q-1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			float_type b = A[a][0] * e_k_arr[(i - 1) * C_M + j] + A[a][2] * e_k_arr[i * C_M + j] + A[a][3] * e_k_arr[a + 1]
				+ A[a][5] * e_k_arr[(i - 1) * C_M + j - 1] + A[a][8] * e_k_arr[(i + 1) * C_M + j + 1];
			e2_arr[a] = e_k_arr[a] - 1 / A[a][2] * (b - f[a]);
			//Для Г7.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			b = A[a][0] * e_k_arr[(i - 1) * C_M + j] + A[a][1] * e_k_arr[i * C_M + j - 1] + A[a][2] * e_k_arr[a]
				+ A[a][6] * e_k_arr[(i - 1) * C_M + j + 1] + A[a][7] * e_k_arr[(i + 1) * C_M + j - 1];
			e2_arr[a] = e_k_arr[a] - 1 / A[a][2] * (b - f[a]);
		}
	} // //#pragma omp parallel 
}

inline int energy(
	cnst_arr_t sigma_k_arr,
	cnst_arr_t sigma_k1_arr,
	cnst_arr_t e_k_arr,
	cnst_ptr_arr_t e_k1_arr,
	cnst_ptr_arr_t e2_arr,	
	cnst_arr_t e_k_mu_arr,
	cnst_arr_t u_k_arr,
	cnst_arr_t v_k_arr, cnst_ptr_arr_t f_arr)
{	
	nrg_calculate_common(A, sigma_k_arr, sigma_k1_arr, e_k_arr, e_k_mu_arr, u_k_arr, v_k_arr, f_arr);
	
	int c; 
	int s_e = 0;
	for (s_e = 0; s_e <= 20; ++s_e)
	{
		nrg_calculate_jakobi(e_k1_arr, e2_arr);
		c = 0;

//#pragma omp parallel for collapse(2) reduction(+:c)
		for (int i = 1; i < C_M1 - 1; i++)
			for (int j = 1; j < C_M - 1; j++)
				if (fabs(e_k1_arr[i * C_M + j] - e2_arr[i * C_M + j]) <= C_epsilon) ++c;

		if (c == C_br) break;

//#pragma omp parallel for collapse(2)
		for (int i = 1; i < C_M1 - 1; i++)
			for (int j = 1; j < C_M - 1; j++)
				e_k1_arr[i * C_M + j] = e2_arr[i * C_M + j];
	}
	return s_e;
}

/*End of energy*/

/* Motion */
/* ----- Функция заполняет элементы матрицы, составленной для двух уравнении движения.---- */
inline void mtn_calculate_common(cnst_ptr_2d_arr_t a_arr, cnst_arr_t sigma_k_arr, cnst_arr_t sigma_k1_arr, cnst_arr_t e_arr, cnst_arr_t e_k_mu_arr, cnst_arr_t u_k_arr, cnst_arr_t v_k_arr, cnst_ptr_arr_t f_arr)
{
	const float_type c_coef1 = C_hx * C_hy * C_Re;
	const float_type c_coef2 = 3 * C_hx * C_hx * C_Re;
	const float_type c_coef3 = 2 * C_hy * C_hy * C_Re;
	const float_type c_coef4 = 2 * C_hx * C_hx * C_Re;
	const float_type c_coef5 = 3 * C_hy * C_hy * C_Re;
	const float_type c_coef6 = 4 * C_hy * C_hy * C_Re;
	const float_type c_coef7 = 6 * C_hx * C_hx * C_Re;
	const float_type c_coef8 = 8 * C_hy * C_hy * C_Re;

	// Уравнение для u	
//#pragma omp parallel 
	{
//#pragma omp for collapse(2) nowait

		for (int i = 1; i < C_M1 - 1; ++i)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; ++j)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * C_M + j;
					a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
					a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3;
					a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
					a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
					a_arr[a][4] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2;
					a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
					a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
					a_arr[a][7] = (e_k_mu_arr[i * C_M + j - 1] / 4 - e_k_mu_arr[(i + 1) * C_M + j] / 6) / c_coef1;
					a_arr[a][8] = (e_k_mu_arr[(i + 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j + 1)] / 4) / c_coef1;
					f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);

					// v
					a = i * C_M + j + C_M2;
					a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
					a_arr[a][1] = -2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5;
					a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
					a_arr[a][3] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
					a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4;
					a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
					a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
					a_arr[a][7] = (e_k_mu_arr[(i + 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j - 1] / 6) / c_coef1;
					a_arr[a][8] = (e_k_mu_arr[i * C_M + j + 1] / 6 - e_k_mu_arr[(i + 1) * C_M + j] / 4) / c_coef1;
				}
			}
		}
//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			// u
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef6 - (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j]) / c_coef8;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3
				+ (e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2 + (e_k_mu_arr[(i + 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j]) / c_coef7
				+ (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j - 1)]) / c_coef6 + (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j - 1)]) / c_coef8;
			a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
			a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 - (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef7;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[(i + 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j + 1)] / 4) / c_coef1;
			a_arr[a][9] = (1. / 4. - 1. / 8.) * e_k_mu_arr[i * C_M + (j - 1)] / c_coef1;
			a_arr[a][11] = (1. / 12. - 1. / 6.) * e_k_mu_arr[(i + 1) * C_M + j] / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
			a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5 - (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j]) / (6 * C_hy * C_hy * C_Re);
			a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5
				+ (e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / (4 * C_hx * C_hx * C_Re) + (e_k_mu_arr[(i + 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j]) / (8 * C_hx * C_hx * C_Re)
				+ (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j - 1)]) / c_coef5 + (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j - 1)]) / (6 * C_hy * C_hy * C_Re);
			a_arr[a][3] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
			a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / (4 * C_hx * C_hx * C_Re) - (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / (8 * C_hx * C_hx * C_Re);
			a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[i * C_M + j + 1] / 6 - e_k_mu_arr[(i + 1) * C_M + j] / 4) / c_coef1;
			a_arr[a][9] = (1. / 12. - 1. / 6.) * e_k_mu_arr[i * C_M + (j - 1)] / c_coef1;
			a_arr[a][11] = (1. / 4. - 1. / 8.) * e_k_mu_arr[(i + 1) * C_M + j] / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3
				+ (e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2 + (e_k_mu_arr[(i + 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j]) / c_coef7
				+ (e_k_mu_arr[i * C_M + (j + 1)] + e_k_mu_arr[i * C_M + j]) / c_coef6 + (e_k_mu_arr[i * C_M + (j + 1)] + 2 * e_k_mu_arr[i * C_M + j]) / c_coef8;
			a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef6 - (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef8;
			a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 - (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef7;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[i * C_M + j - 1] / 4 - e_k_mu_arr[(i + 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][10] = (1. / 8. - 1. / 4.) * e_k_mu_arr[i * C_M + (j + 1)] / c_coef1;
			a_arr[a][11] = (1. / 6. - 1. / 12.) * e_k_mu_arr[(i + 1) * C_M + j] / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
			a_arr[a][1] = -2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5;
			a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5
				+ (e_k_mu_arr[(i + 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / (4 * C_hx * C_hx * C_Re) + (e_k_mu_arr[(i + 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j]) / (8 * C_hx * C_hx * C_Re)
				+ (e_k_mu_arr[i * C_M + (j + 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5 + (e_k_mu_arr[i * C_M + (j + 1)] + 2 * e_k_mu_arr[i * C_M + j]) / (6 * C_hy * C_hy * C_Re);
			a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5 - (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / (6 * C_hy * C_hy * C_Re);
			a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / (4 * C_hx * C_hx * C_Re) - (2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / (8 * C_hx * C_hx * C_Re);
			a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[(i + 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j - 1] / 6) / c_coef1;
			a_arr[a][10] = (1. / 6. - 1. / 12.) * e_k_mu_arr[i * C_M + (j + 1)] / c_coef1;
			a_arr[a][11] = (1. / 8. - 1. / 4.) * e_k_mu_arr[(i + 1) * C_M + j] / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
		}
//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				// u
				int a = i * C_M + j;
				a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
				a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3;
				a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
				a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
				a_arr[a][4] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2;
				a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[i * C_M + j - 1] / 4 - e_k_mu_arr[(i + 1) * C_M + j] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[(i + 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j + 1)] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
				// v
				a += C_M2;
				a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
				a_arr[a][1] = -2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5;
				a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
				a_arr[a][3] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
				a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4;
				a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[(i + 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j - 1] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[i * C_M + j + 1] / 6 - e_k_mu_arr[(i + 1) * C_M + j] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// u
				int a = i * C_M + j;
				a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
				a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3;
				a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
				a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
				a_arr[a][4] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2;
				a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[i * C_M + j - 1] / 4 - e_k_mu_arr[(i + 1) * C_M + j] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[(i + 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j + 1)] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
				// v
				a += C_M2;
				a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
				a_arr[a][1] = -2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5;
				a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
				a_arr[a][3] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
				a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4;
				a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[(i + 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j - 1] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[i * C_M + j + 1] / 6 - e_k_mu_arr[(i + 1) * C_M + j] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			}
		}
//#pragma omp single nowait
		{
			// u
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
			a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
			a_arr[a][4] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[i * C_M + j - 1] / 4 - e_k_mu_arr[(i + 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[(i + 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j + 1)] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
			a_arr[a][1] = -2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5;
			a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
			a_arr[a][3] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
			a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4;
			a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[(i + 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j - 1] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[i * C_M + j + 1] / 6 - e_k_mu_arr[(i + 1) * C_M + j] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef3;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau + 2 * (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
			a_arr[a][3] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef3;
			a_arr[a][4] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef2;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * C_M + j + 1] / 4 - e_k_mu_arr[(i - 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[i * C_M + j - 1] / 4 - e_k_mu_arr[(i + 1) * C_M + j] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[(i + 1) * C_M + j] / 6 - e_k_mu_arr[i * C_M + (j + 1)] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * C_M + j] + e_k_mu_arr[i * C_M + j]) / c_coef4;
			a_arr[a][1] = -2 * (e_k_mu_arr[i * C_M + (j - 1)] + e_k_mu_arr[i * C_M + j]) / c_coef5;
			a_arr[a][2] = sigma_k1_arr[i * C_M + j] * sigma_k1_arr[i * C_M + j] / C_tau + (e_k_mu_arr[(i - 1) * C_M + j] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu_arr[i * C_M + (j - 1)] + 2 * e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
			a_arr[a][3] = -2 * (e_k_mu_arr[i * C_M + j] + e_k_mu_arr[i * C_M + (j + 1)]) / c_coef5;
			a_arr[a][4] = -(e_k_mu_arr[i * C_M + j] + e_k_mu_arr[(i + 1) * C_M + j]) / c_coef4;
			a_arr[a][5] = (e_k_mu_arr[i * C_M + j - 1] / 6 - e_k_mu_arr[(i - 1) * C_M + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j + 1] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[(i + 1) * C_M + j] / 4 - e_k_mu_arr[i * C_M + j - 1] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[i * C_M + j + 1] / 6 - e_k_mu_arr[(i + 1) * C_M + j] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau - (P(C_gamma, sigma_k_arr[(i + 1) * C_M + j], e_arr[(i + 1) * C_M + j]) - P(C_gamma, sigma_k_arr[(i - 1) * C_M + j], e_arr[(i - 1) * C_M + j])) / (2 * C_hx);
		}
	} // //#pragma omp parallel
}

//Вектор B = A*Xk1
// указывая __restrict мы подсказываем компилятору, что области памяти не пересекаются
// то есть нет явления memory aliasing
// тогда компилятор может применить автоматическую векторизацию цикла
// таким образом уничтожаются зависимости типа ANTI и FLOW 
inline void mtn_calculate_jakobi(cnst_arr_t u_arr, cnst_arr_t v_arr, cnst_ptr_arr_t u2_arr, cnst_ptr_arr_t v2_arr, cnst_arr_t f_arr)
{
//#pragma omp parallel
	{
//#pragma omp for collapse(2) nowait
//#pragma omp for simd nowait
//#pragma omp simd collapse(2)
		__assume((C_M - 1) % 4 == 0);
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * C_M + j;
					float_type b = A[a][0] * u_arr[(i - 1) * C_M + j] + A[a][1] * u_arr[i * C_M + (j - 1)] + A[a][2] * u_arr[i * C_M + j] +
						A[a][3] * u_arr[i * C_M + (j + 1)] + A[a][4] * u_arr[(i + 1) * C_M + j] +
						A[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
						A[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
						A[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
						A[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
					u2_arr[a] = u_arr[a] - 1 / A[a][2] * (b - f_arr[a]);
					// v
					 a = i * C_M + j + C_M2;
					 b = A[a][0] * v_arr[(i - 1) * C_M + j] + A[a][1] * v_arr[i * C_M + j - 1] + A[a][2] * v_arr[i * C_M + j] +
						A[a][3] * v_arr[i * C_M + j + 1] + A[a][4] * v_arr[(i + 1) * C_M + j] +
						A[a][5] * u_arr[(i - 1) * C_M + j - 1] +
						A[a][6] * u_arr[(i - 1) * C_M + j + 1] +
						A[a][7] * u_arr[(i + 1) * C_M + j - 1] +
						A[a][8] * u_arr[(i + 1) * C_M + j + 1];
					v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / A[a][2] * (b - f_arr[a]);
				}
			}
		}
//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				// u
				int a = i * C_M + j;
				float_type b = A[a][0] * u_arr[(i - 1) * C_M + j] + A[a][1] * u_arr[i * C_M + (j - 1)] + A[a][2] * u_arr[i * C_M + j] +
					A[a][3] * u_arr[i * C_M + (j + 1)] + A[a][4] * u_arr[(i + 1) * C_M + j] +
					A[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
				u2_arr[a] = u_arr[a] - 1 / A[a][2] * (b - f_arr[a]);
				// v
				a += C_M2;
				b = A[a][0] * v_arr[(i - 1) * C_M + j] + A[a][1] * v_arr[i * C_M + j - 1] + A[a][2] * v_arr[i * C_M + j] +
					A[a][3] * v_arr[i * C_M + j + 1] + A[a][4] * v_arr[(i + 1) * C_M + j] +
					A[a][5] * u_arr[(i - 1) * C_M + j - 1] +
					A[a][6] * u_arr[(i - 1) * C_M + j + 1] +
					A[a][7] * u_arr[(i + 1) * C_M + j - 1] +
					A[a][8] * u_arr[(i + 1) * C_M + j + 1];
				v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / A[a][2] * (b - f_arr[a]);
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// u
				int a = i * C_M + j;
				float_type b = A[a][0] * u_arr[(i - 1) * C_M + j] + A[a][1] * u_arr[i * C_M + (j - 1)] + A[a][2] * u_arr[i * C_M + j] +
					A[a][3] * u_arr[i * C_M + (j + 1)] + A[a][4] * u_arr[(i + 1) * C_M + j] +
					A[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
				u2_arr[a] = u_arr[a] - 1 / A[a][2] * (b - f_arr[a]);
				// v
				a += C_M2;
				b = A[a][0] * v_arr[(i - 1) * C_M + j] + A[a][1] * v_arr[i * C_M + j - 1] + A[a][2] * v_arr[i * C_M + j] +
					A[a][3] * v_arr[i * C_M + j + 1] + A[a][4] * v_arr[(i + 1) * C_M + j] +
					A[a][5] * u_arr[(i - 1) * C_M + j - 1] +
					A[a][6] * u_arr[(i - 1) * C_M + j + 1] +
					A[a][7] * u_arr[(i + 1) * C_M + j - 1] +
					A[a][8] * u_arr[(i + 1) * C_M + j + 1];
				v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / A[a][2] * (b - f_arr[a]);
			}
		}
//#pragma omp single nowait
		{
			// u
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			float_type b = A[a][0] * u_arr[(i - 1) * C_M + j] + A[a][1] * u_arr[i * C_M + (j - 1)] + A[a][2] * u_arr[i * C_M + j] +
				A[a][3] * u_arr[i * C_M + (j + 1)] + A[a][4] * u_arr[(i + 1) * C_M + j] +
				A[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
			u2_arr[a] = u_arr[a] - 1 / A[a][2] * (b - f_arr[a]);
			// v
			a += C_M2;
			b = A[a][0] * v_arr[(i - 1) * C_M + j] + A[a][1] * v_arr[i * C_M + j - 1] + A[a][2] * v_arr[i * C_M + j] +
				A[a][3] * v_arr[i * C_M + j + 1] + A[a][4] * v_arr[(i + 1) * C_M + j] +
				A[a][5] * u_arr[(i - 1) * C_M + j - 1] +
				A[a][6] * u_arr[(i - 1) * C_M + j + 1] +
				A[a][7] * u_arr[(i + 1) * C_M + j - 1] +
				A[a][8] * u_arr[(i + 1) * C_M + j + 1];
			v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / A[a][2] * (b - f_arr[a]);
			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			b = A[a][0] * u_arr[(i - 1) * C_M + j] + A[a][1] * u_arr[i * C_M + (j - 1)] + A[a][2] * u_arr[i * C_M + j] +
				A[a][3] * u_arr[i * C_M + (j + 1)] + A[a][4] * u_arr[(i + 1) * C_M + j] +
				A[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
			u2_arr[a] = u_arr[a] - 1 / A[a][2] * (b - f_arr[a]);
			// v
			a += C_M2;
			b = A[a][0] * v_arr[(i - 1) * C_M + j] + A[a][1] * v_arr[i * C_M + j - 1] + A[a][2] * v_arr[i * C_M + j] +
				A[a][3] * v_arr[i * C_M + j + 1] + A[a][4] * v_arr[(i + 1) * C_M + j] +
				A[a][5] * u_arr[(i - 1) * C_M + j - 1] +
				A[a][6] * u_arr[(i - 1) * C_M + j + 1] +
				A[a][7] * u_arr[(i + 1) * C_M + j - 1] +
				A[a][8] * u_arr[(i + 1) * C_M + j + 1];
			v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / A[a][2] * (b - f_arr[a]);
		}
	} // //#pragma omp parallel
}

inline int motion(cnst_arr_t sigma_k_arr,
				  cnst_arr_t sigma_k1_arr,
                  cnst_arr_t u_k_arr,
                  cnst_arr_t v_k_arr,
				  cnst_ptr_arr_t f_arr,
                  cnst_ptr_arr_t u_k1_arr,
                  cnst_ptr_arr_t v_k1_arr,
                  cnst_ptr_arr_t u2_arr,
                  cnst_ptr_arr_t v2_arr,
                  cnst_arr_t e_k_arr,
                  cnst_arr_t e_k_mu_arr)
{	
	mtn_calculate_common(A, sigma_k_arr, sigma_k1_arr, e_k_arr, e_k_mu_arr, u_k_arr, v_k_arr, f_arr);

	int c_u;
	int c_v;
	int s_m = 0;
	for (s_m = 0; s_m <= 20; ++s_m)
	{
		mtn_calculate_jakobi(u_k1_arr, v_k1_arr, u2_arr, v2_arr, f_arr);

		c_u = 0;
		c_v = 0;

//#pragma omp parallel for collapse(2) reduction(+:c_u, c_v)
		for (int i = 1; i < C_M1 - 1; i++)
			for (int j = 1; j < C_M - 1; j++)
			{
				if (fabs(u_k1_arr[i * C_M + j] - u2_arr[i * C_M + j]) <= C_epsilon) ++c_u;
				if (fabs(v_k1_arr[i * C_M + j] - v2_arr[i * C_M + j]) <= C_epsilon) ++c_v;
			}

		if (c_u >= C_br && c_v >= C_br) break;

		// Можно копировать с memcpy?
//#pragma omp parallel for collapse(2)
		for (int i = 1; i < C_M1 - 1; i++)
			for (int j = 1; j < C_M - 1; j++)
			{
				u_k1_arr[i * C_M + j] = u2_arr[i * C_M + j];
				v_k1_arr[i * C_M + j] = v2_arr[i * C_M + j];
			}
	}
	return s_m;
}

/* End of motion */

inline int interate_over_nonlinearity(int& s_m, int& s_e, int& s_end)
{
	const int itr = 5;

	int s_itr;
	for (s_itr = 1; s_itr < itr; ++s_itr)
	{
		for (int i = 0; i < C_M2; i++)
			e_k_mu[i] = Mu(e_k[i]);

		continuity(sigma_kk, sigma_k1, u_k, v_k);
		s_m = motion(sigma_k, sigma_k1, u_k, v_k, f, u_k1, v_k1, u2, v2, e_k, e_k_mu);
		s_e = energy(sigma_k, sigma_k1, e_k, e_k1, e2, e_k_mu, u_k, v_k, f);

		if (s_m == 1 && s_e == 1)
		{
			s_end = s_itr;
			s_itr = itr;
		}

		memcpy(sigma_k, sigma_k1, C_M * sizeof *sigma_k);
		memcpy(e_k, e_k1, C_M * sizeof *e_k);
		memcpy(u_k, u_k1, C_M * sizeof *u_k);
		memcpy(v_k, v_k1, C_M * sizeof *v_k);

//#pragma omp parallel for collapse(2) 
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 0; j < C_M; j++)
			{
				if (i < C_qq + 1 || i >= C_qq + C_w - 1)
				{
					int a = i * C_M + j;

					int idx = a;

					if (j == 0) idx++; else if (j == C_N) idx--;

					if (j >= 0 && j <= C_N)
					{
						sigma_k[a] = sigma_k1[idx];
						e_k[a] = e_k1[idx];
						u_k[a] = u_k1[idx];
						v_k[a] = v_k1[idx];
						if (j == 0 || j == C_N)
						{
							sigma_k[a] = sigma_k1[idx];
							sigma_k1[a] = sigma_k1[idx];
							e_k1[a] = e_k1[idx];
							u_k1[a] = u_k1[idx];
							v_k1[a] = v_k1[idx];
							e2[a] = e_k1[idx];
							u2[a] = u_k1[idx];
							v2[a] = v_k1[idx];
						}
					}
				}
			}
		}
//#pragma omp parallel for
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			for (int j = C_cntr + i - C_qq; j < C_M; j++)
			{
				int a = i * C_M + j;
				int idx = j == C_N ? a - 1 : a;

				sigma_k[a] = sigma_k1[idx];
				e_k[a] = e_k1[idx];
				u_k[a] = u_k1[idx];
				v_k[a] = v_k1[idx];

				if (j == C_N)
				{
					sigma_k1[a] = sigma_k1[idx];
					e_k1[a] = e_k1[idx];
					u_k1[a] = u_k1[idx];
					v_k1[a] = v_k1[idx];
					e2[a] = e_k1[idx];
					u2[a] = u_k1[idx];
					v2[a] = v_k1[idx];
				}
			}
			for (int j = C_cntr - i + C_qq; j >= 0; j--)
			{
				int a = i * C_M + j;
				int idx = j == 0 ? a + 1 : a;
				sigma_k[a] = sigma_k1[idx];
				e_k[a] = e_k1[idx];
				u_k[a] = u_k1[idx];
				v_k[a] = v_k1[idx];
				if (j == 0)
				{
					sigma_k1[a] = sigma_k1[idx];
					e_k1[a] = e_k1[idx];
					u_k1[a] = u_k1[idx];
					v_k1[a] = v_k1[idx];
					e2[a] = e_k1[idx];
					u2[a] = u_k1[idx];
					v2[a] = v_k1[idx];
				}
			}
		}
//#pragma omp parallel for
		for (int j = 0; j < C_M; j++)
		{
			int a = C_N1 * C_M + j;
			int indx = (C_N1 - 1) * C_M + j;
			if (j == 0) indx++; else if (j == C_N) indx--;

			sigma_k[a] = sigma_k1[indx + 1];
			sigma_k1[a] = sigma_k1[indx + 1];
			e_k[a] = e_k1[indx + 1];
			e_k1[a] = e_k1[indx + 1];
			u_k[a] = u_k1[indx + 1];
			u_k1[a] = u_k1[indx + 1];
			v_k[a] = v_k1[indx + 1];
			v_k1[a] = v_k1[indx + 1];
			e2[a] = e_k1[indx + 1];
			u2[a] = u_k1[indx + 1];
			v2[a] = v_k1[indx + 1];
		}
	}
	return s_itr;
}

inline void init_arrays_parallel(const int array_element_count, const int param_array_element_count)
{
	int double_size_array = 2 * array_element_count;
	A = new float_type*[double_size_array];
	for (int i = 0; i < double_size_array; ++i)
	{
		A[i] = new float_type[param_array_element_count];
		std::fill_n(A[i], param_array_element_count, 0.);
	}		
	
	f = new float_type[double_size_array];
	std::fill_n(f, double_size_array, 0.);

	sigma_k = new float_type[array_element_count];
	sigma_k1 = new float_type[array_element_count];
	sigma_kk = new float_type[array_element_count];
	std::fill_n(sigma_k, array_element_count, 0.);
	std::fill_n(sigma_k1, array_element_count, 0.);
	std::fill_n(sigma_kk, array_element_count, 0.);

	u_k = new float_type[array_element_count];
	u_k1 = new float_type[array_element_count];
	u_kk = new float_type[array_element_count];
	u2 = new float_type[array_element_count];
	std::fill_n(u_k, array_element_count, 0.);
	std::fill_n(u_k1, array_element_count, 0.);
	std::fill_n(u_kk, array_element_count, 0.);
	std::fill_n(u2, array_element_count, 0.);

	v_k = new float_type[array_element_count];
	v_k1 = new float_type[array_element_count];
	v_kk = new float_type[array_element_count];
	v2 = new float_type[array_element_count];
	std::fill_n(v_k, array_element_count, 0.);
	std::fill_n(v_k1, array_element_count, 0.);
	std::fill_n(v_kk, array_element_count, 0.);
	std::fill_n(v2, array_element_count, 0.);

	e_k = new float_type[array_element_count];
	e_k1 = new float_type[array_element_count];
	e_kk = new float_type[array_element_count];
	e2 = new float_type[array_element_count];
	e_k_mu = new float_type[array_element_count];
	std::fill_n(e_k, array_element_count, 0.);
	std::fill_n(e_k1, array_element_count, 0.);
	std::fill_n(e_kk, array_element_count, 0.);
	std::fill_n(e2, array_element_count, 0.);
	std::fill_n(e_k_mu, array_element_count, 0.);
}
