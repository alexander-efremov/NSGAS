#include "par_calculate.h"
#include "print.h"
#include "par_api.h"
#include "timer.h"

//C_M1 - êîëè÷åñòâî óçëîâ ïî îñè õ
//C_M - êîëè÷åñòâî óçëîâ ïî îñè y
//(2*C_q-1) - êîëè÷åñòâî óçëîâ â îñíîâàíèè êëèíà
//(C_qq,C_cntr) - íîìåð óçëà âåðøèíû êëèíà
//C_tg = C_hx/C_hy

// test case
//static const int C_N = 20;
//static const int C_N1 = 10;
//static const int C_q = 3;
//static const int C_qq = 5;
//static const float_type C_hx = 1.0 / 10; // 1.0 / C_N1
//static const float_type C_hy = 1.0 / 20; // 1.0 / C_N
//static const int C_M = 21; // C_N + 1
//static const int C_M1 = 11; // C_N1 + 1
//static const int C_M2 = 231; // C_M1 * C_M 
//static const int C_w = 3; // C_w = C_q
//static const int C_cntr = 10; // C_N / 2
//static const int C_br = 171; // (C_N1 - 1) * (C_N - 1)
//static const int C_cnt_boundary_nodes = 60; // C_M2 - C_br
//-----------------------
// real test 
static const int C_N = 1200;
static const int C_N1 = 800;
static const int C_q = 101;
static const int C_qq = 20;
static const float_type C_hx = 1.0 / 100;
static const float_type C_hy = 1.0 / 200;
static const int C_M = 1201; // C_N + 1
static const int C_M1 = 801; // C_N1 + 1
static const int C_M2 = 962001; // C_M1 * C_M 
static const int C_w = 101; // C_w = C_q
static const int C_cntr = 600; // C_N / 2
static const int C_br = 958001; // (C_N1 - 1) * (C_N - 1)
static const int C_cnt_boundary_nodes = 4000; // C_M2 - C_br

//static const int time_steps_nbr = 25000; // time_steps_nbr - êîëè÷åñòâî øàãîâ ïî âðåìåíè
static const int time_steps_nbr = 1; // time_steps_nbr - êîëè÷åñòâî øàãîâ ïî âðåìåíè
//-----------------------

static const float_type C_PrRe = 7200; // Pr * C_Re = 0.72 * C_Re
static const float_type C_gamma_Mah2 = 8.96; // C_gamma * (C_gamma - 1) * C_Mah2
static const float_type C_Mah2 = 16; // Mah * Mah
static const float_type C_epsilon = 0.0000000001;
static const float_type C_gamma = 1.4;
static const float_type C_tau = 0.0005;
static const float_type C_tg = 2;
static const float_type C_Re = 10000;


// Â ìàññèâàõ ñ _k õðàíÿòñÿ çíà÷åíèÿ ôóíêöèé ñ ïðåäûäóùåé èòåðàöèè ïî íåëèíåéíîñòè
// Â ìàññèâàõ ñ _kk õðàíÿòñÿ çíà÷åíèÿ ôóíêöèé c (d-1) øàãà ïî âðåìåíè
// Â ìàññèâàõ ñ _k1 õðàíÿòñÿ çíà÷åíèÿ ôóíêöèé íà d-îì øàãå ïî âðåìåíè
// Ìàññèâû ñ "2" èñïîëüçóòñÿ â èòåðàöèÿõ ìåòîäà ßêîáè

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
static float_type* v2;

static float_type* e_k;
static float_type* e_k1;
static float_type* e_kk;
static float_type* e2;
static float_type* e_k_mu;

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

#pragma omp declare simd uniform(tau_d, hx_d, hy_d, u_k_value, v_k_value, m_i) linear(i, j)
inline float_type trajectory(const float_type tau_d, const float_type hx_d, const float_type hy_d, int i, int j, cnst_arr_t arr, const float_type u_k_value, const float_type v_k_value, const int m_i)
{
	int idx = i * m_i + j;
	int idx2 = i * m_i + (j - 1);
	int idx3 = (i - 1) * m_i + j;
	int idx4 = (i + 1) * m_i + j;
	int idx5 = i * m_i + (j + 1);

	if (u_k_value == 0. && v_k_value == 0.) return 0;
	if (u_k_value > 0. && v_k_value > 0.)
		return u_k_value * tau_d * ((arr[idx3] - arr[idx]) / ((i - 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx2] - arr[idx]) / ((j - 1) * hy_d - j * hy_d));
	if (u_k_value < 0. && v_k_value > 0.)
		return u_k_value * tau_d * ((arr[idx4] - arr[idx]) / ((i + 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx2] - arr[idx]) / ((j - 1) * hy_d - j * hy_d));
	if (u_k_value < 0. && v_k_value < 0.)
		return u_k_value * tau_d * ((arr[idx4] - arr[idx]) / ((i + 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx5] - arr[idx]) / ((j + 1) * hy_d - j * hy_d));
	if (u_k_value > 0. && v_k_value < 0.)
		return u_k_value * tau_d * ((arr[idx3] - arr[idx]) / ((i - 1) * hx_d - i * hx_d))
			+ v_k_value * tau_d * ((arr[idx5] - arr[idx]) / ((j + 1) * hy_d - j * hy_d));
	if (u_k_value > 0. && v_k_value == 0.)
		return u_k_value * tau_d * ((arr[idx3] - arr[idx]) / ((i - 1) * hx_d - i * hx_d));
	if (u_k_value == 0. && v_k_value > 0.)
		return v_k_value * tau_d * ((arr[idx2] - arr[idx]) / ((j - 1) * hy_d - j * hy_d));
	if (u_k_value < 0. && v_k_value == 0.)
		return u_k_value * tau_d * ((arr[idx4] - arr[idx]) / ((i + 1) * hx_d - i * hx_d));
	if (u_k_value == 0. && v_k_value < 0.)
		return v_k_value * tau_d * ((arr[idx5] - arr[idx]) / ((j + 1) * hy_d - j * hy_d));
	return 0;
}

inline void continuity(cnst_arr_t sigma_kk_arr, cnst_ptr_arr_t sigma_k1_arr, cnst_arr_t u_k_arr, cnst_arr_t v_k_arr)
{
	//Äëÿ âíóòðåííèõ óçëîâ
	//#pragma omp parallel
	{
		//#pragma omp for nowait
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
			//Äëÿ Ã5. l = C_w-1; m = 1,...,C_q-2;
			int i = C_qq + C_w - 1;
			sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) / (2 * C_tau) / (1 / (2 * C_tau) + (u_k_arr[(i + 1) * C_M + j] - u_k_arr[i * C_M + j]) / (4 * C_hx)
				+ (v_k_arr[i * C_M + j + 1] - v_k_arr[i * C_M + j - 1]) / (8 * C_hy));
		}
		//#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//Äëÿ Ã6.
			int j = C_cntr + i - C_qq;
			sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) * (1 / (4 * C_tau) + 1 / (4 * C_tau)) / (1 / (4 * C_tau) + 1 / (4 * C_tau) + (u_k_arr[i * C_M + j] - u_k_arr[(i - 1) * C_M + j]) / (8 * C_hx)
				- u_k_arr[(i - 1) * C_M + j] / (16 * C_hx) + (v_k_arr[i * C_M + j + 1] - v_k_arr[i * C_M + j]) / (8 * C_hy) + v_k_arr[i * C_M + j + 1] / (16 * C_hy));
			//Äëÿ Ã7.
			j = C_cntr - i + C_qq;
			sigma_k1_arr[i * C_M + j] = (sigma_kk_arr[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[i * C_M + j], v_k_arr[i * C_M + j], C_M)) * (1 / (4 * C_tau) + 1 / (4 * C_tau)) / (1 / (4 * C_tau) + 1 / (4 * C_tau) + (u_k_arr[i * C_M + j] - u_k_arr[(i - 1) * C_M + j]) / (8 * C_hx)
				- u_k_arr[(i - 1) * C_M + j] / (16 * C_hx) + (v_k_arr[i * C_M + j] - v_k_arr[i * C_M + j - 1]) / (8 * C_hy) - v_k_arr[i * C_M + j - 1] / (16 * C_hy));
		}
		//#pragma omp single
		{
			//Äëÿ S_w-1q-1.
			int i = C_qq + C_w - 1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			sigma_k1_arr[a] = (sigma_kk_arr[a] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[a], v_k_arr[a], C_M)) * (3 / (4 * C_tau) + 1 / (8 * C_tau)) / (3 / (4 * C_tau) + 1 / (8 * C_tau) + (2 * u_k_arr[(i + 1) * C_M + j] - u_k_arr[(i - 1) * C_M + j] - u_k_arr[a]) / (8 * C_hx)
				+ (u_k_arr[a] - u_k_arr[(i - 1) * C_M + j]) / (16 * C_hx) + (2 * v_k_arr[a + 1] - v_k_arr[a - 1] - v_k_arr[a]) / (8 * C_hy) + v_k_arr[a] / (16 * C_hy));
			//Äëÿ S.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			sigma_k1_arr[a] = (sigma_kk_arr[a] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[a], v_k_arr[a], C_M)) * (3 / (4 * C_tau) + 1 / (8 * C_tau)) / (3 / (4 * C_tau) + 1 / (8 * C_tau) + (2 * u_k_arr[(i + 1) * C_M + j] - u_k_arr[(i - 1) * C_M + j] - u_k_arr[a]) / (8 * C_hx)
				+ (u_k_arr[a] - u_k_arr[(i - 1) * C_M + j]) / (16 * C_hx) + (v_k_arr[a + 1] - 2 * v_k_arr[a - 1] + v_k_arr[a]) / (8 * C_hy) - v_k_arr[a] / (16 * C_hy));
			//Äëÿ S_qq_0
			i = C_qq;
			j = C_cntr;
			a = i * C_M + j;
			sigma_k1_arr[a] = (sigma_kk_arr[a] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk_arr, u_k_arr[a], v_k_arr[a], C_M)) * (1 / (4 * C_tau) + 1 / (2 * C_tau)) / (1 / (4 * C_tau) + 1 / (2 * C_tau) + (u_k_arr[a] - u_k_arr[(i - 1) * C_M + j]) / (4 * C_hx)
				+ (u_k_arr[(i + 1) * C_M + j] - u_k_arr[a]) / (8 * C_hx) + (v_k_arr[a + 1] - v_k_arr[a - 1]) / (8 * C_hy) + (v_k_arr[a + 1] - v_k_arr[a - 1]) / (16 * C_hy));
		}
	} // //#pragma omp parallel
}

/*Energy*/
/*----- Ôóíêöèÿ çàïîëíÿåò ýëåìåíòû ìàòðèöû äëÿ óðàâíåíèÿ ýíåðãèè----*/
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
		// следующие два цикла разделены на 2, хотя могут выполняться в одном
		// сделано это затем, чтобы убрать предупреждение Intel Adviser'а 
		// High vector register pressure 
		// Recommendation: Split loop into smaller loops 
		// И действительно:
		// До разделения: 218 ms
		// После разделения: 156 ms
		//#pragma omp for nowait
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
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;
					f_arr[a] = (e_kk[a] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k_arr[a], v_k_arr[a], C_M)) * sigma_k1_arr[a] * sigma_k1_arr[a] / C_tau;
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;
					float_type val6 = u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j];
					f_arr[a] -= P(C_gamma, sigma_k_arr[a], e_k_arr[a]) / (4 * e_k_arr[a]) * (val6 / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy);
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;

					float_type val5 = u_k1[a] - u_k1[(i - 1) * C_M + j];
					float_type val7 = u_k1[(i + 1) * C_M + j] - u_k1[a];
					float_type val8 = v_k1[a + 1] - v_k1[a];
					float_type val9 = v_k1[a] - v_k1[a - 1];
					f_arr[a] +=
						e_k_mu_arr[a] / (6 * C_hx * C_hx * C_Re * e_k_arr[a]) * (val7 * val7 + val5 * val5) +
						e_k_mu_arr[a] / (6 * C_hy * C_hy * C_Re * e_k_arr[a]) * (val8 * val8 + val9 * val9) ;
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;

					float_type val1 = v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy;
					float_type val2 = v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy;
					float_type val3 = v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy;
					float_type val4 = v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy;
					f_arr[a] +=
						e_k_mu_arr[a] / (8 * C_Re * e_k_arr[a]) * (val1 * val1 + val2 * val2 + val3 * val3 + val4 * val4);
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;

					float_type val1 = u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy;
					float_type val2 = u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy;
					float_type val3 = u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy;
					float_type val4 = u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy;
					f_arr[a] += e_k_mu_arr[a] / (12 * C_Re * e_k_arr[a]) * (val1 * val1 + val2 * val2 + val3 * val3 + val4 * val4);
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
			//Äëÿ Ã5. l = C_q-1; C_M = 1,...,C_q-1;	
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
			//Äëÿ Ã6. l = 1,...,C_q-1; C_M = C_q-1;
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

			// Äëÿ Ã7.
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

			//Äëÿ S_qq, C_N / 2 + C_q.
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

			//Äëÿ S_qq, C_N/2 - C_q.
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

			//Äëÿ S_qq,C_N/2
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

//Âåêòîð B = A*Xk1
inline void nrg_calculate_jakobi(cnst_arr_t e_k_arr, cnst_ptr_arr_t e2_arr)
{
	//#pragma omp parallel 
	{
		//Äëÿ âíóòðåííèõ óçëîâ
		//#pragma omp for nowait
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
			//Äëÿ Ã6. l = 1,...,C_q-1; C_M = C_q-1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			float_type b = A[a][0] * e_k_arr[(i - 1) * C_M + j] + A[a][2] * e_k_arr[i * C_M + j] + A[a][3] * e_k_arr[a + 1]
				+ A[a][5] * e_k_arr[(i - 1) * C_M + j - 1] + A[a][8] * e_k_arr[(i + 1) * C_M + j + 1];
			e2_arr[a] = e_k_arr[a] - 1 / A[a][2] * (b - f[a]);
			//Äëÿ Ã7.
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

		//#pragma omp parallel for reduction(+:c)
#pragma vector aligned
		for (int i = 0; i < C_M2; i++)
			if (fabs(e_k1_arr[i] - e2_arr[i]) <= C_epsilon) ++c;

		if (c - C_cnt_boundary_nodes >= C_br)
		{
			s_e++;
			return s_e;
		}

		//#pragma omp parallel for 
		for (int i = 1; i < C_M1 - 1; i++)
			for (int j = 1; j < C_M - 1; j++)
				e_k1_arr[i * C_M + j] = e2_arr[i * C_M + j];
	}
	return s_e;
}

/*End of energy*/

/* Motion */
/* ----- Ôóíêöèÿ çàïîëíÿåò ýëåìåíòû ìàòðèöû, ñîñòàâëåííîé äëÿ äâóõ óðàâíåíèè äâèæåíèÿ.---- */
inline void mtn_calculate_common(cnst_ptr_2d_arr_t a_arr, cnst_arr_t sigma_k_arr, cnst_arr_t sigma_k1_arr, cnst_arr_t e_arr, cnst_arr_t e_k_mu_arr, cnst_arr_t u_k_arr, cnst_arr_t v_k_arr, cnst_ptr_arr_t f_arr,
                                 const float_type c_tau_d,
                                 const float_type c_hx_d,
                                 const float_type c_hy_d,
                                 const int c_m_d
)
{
	const float_type c_coef1 = c_hx_d * c_hy_d * C_Re;
	const float_type c_coef2 = 3 * c_hx_d * c_hx_d * C_Re;
	const float_type c_coef3 = 2 * c_hy_d * c_hy_d * C_Re;
	const float_type c_coef4 = 2 * c_hx_d * c_hx_d * C_Re;
	const float_type c_coef5 = 3 * c_hy_d * c_hy_d * C_Re;
	const float_type c_coef6 = 4 * c_hy_d * c_hy_d * C_Re;
	const float_type c_coef7 = 6 * c_hx_d * c_hx_d * C_Re;
	const float_type c_coef8 = 8 * c_hy_d * c_hy_d * C_Re;

	// Óðàâíåíèå äëÿ u	
	//#pragma omp parallel 
	{
		// следующие два цикла разделены на 2, хотя могут выполняться в одном
		// сделано это затем, чтобы убрать предупреждение Intel Adviser'а 
		// High vector register pressure 
		// Recommendation: Split loop into smaller loops 
		// И действительно:
		// До разделения: 440 ms
		// После разделения: 172 ms
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; ++i)
		{
#pragma ivdep
			for (int j = 1; j < c_m_d - 1; ++j)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * c_m_d + j;
					a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
					a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3;
					a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
					a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
					a_arr[a][4] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2;
					a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
					a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
					a_arr[a][7] = (e_k_mu_arr[i * c_m_d + j - 1] / 4 - e_k_mu_arr[(i + 1) * c_m_d + j] / 6) / c_coef1;
					a_arr[a][8] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j + 1)] / 4) / c_coef1;
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; ++i)
		{
#pragma ivdep
			for (int j = 1; j < c_m_d - 1; ++j)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * c_m_d + j;
					f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d ;
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; ++i)
		{
#pragma ivdep
			for (int j = 1; j < c_m_d - 1; ++j)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * c_m_d + j;
					f_arr[a] -=
						(P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
				}
			}
		}
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; ++i)
		{
#pragma ivdep
			for (int j = 1; j < c_m_d - 1; ++j)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// v
					int a = i * c_m_d + j + C_M2;
					a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
					a_arr[a][1] = -2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5;
					a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
					a_arr[a][3] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
					a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4;
					a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
					a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
					a_arr[a][7] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j - 1] / 6) / c_coef1;
					a_arr[a][8] = (e_k_mu_arr[i * c_m_d + j + 1] / 6 - e_k_mu_arr[(i + 1) * c_m_d + j] / 4) / c_coef1;
				}
			}
		}
		//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			// u
			int j = C_cntr + i + 1 - C_qq;
			int a = i * c_m_d + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef6 - (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j]) / c_coef8;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3
				+ (e_k_mu_arr[(i + 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2 + (e_k_mu_arr[(i + 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j]) / c_coef7
				+ (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j - 1)]) / c_coef6 + (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j - 1)]) / c_coef8;
			a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
			a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 - (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef7;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j + 1)] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
			a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5 - (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j]) / (6 * c_hy_d * c_hy_d * C_Re);
			a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5
				+ (e_k_mu_arr[(i + 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / (4 * c_hx_d * c_hx_d * C_Re) + (e_k_mu_arr[(i + 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j]) / (8 * c_hx_d * c_hx_d * C_Re)
				+ (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j - 1)]) / c_coef5 + (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j - 1)]) / (6 * c_hy_d * c_hy_d * C_Re);
			a_arr[a][3] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
			a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / (4 * c_hx_d * c_hx_d * C_Re) - (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / (8 * c_hx_d * c_hx_d * C_Re);
			a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[i * c_m_d + j + 1] / 6 - e_k_mu_arr[(i + 1) * c_m_d + j] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * c_m_d + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3
				+ (e_k_mu_arr[(i + 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2 + (e_k_mu_arr[(i + 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j]) / c_coef7
				+ (e_k_mu_arr[i * c_m_d + (j + 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef6 + (e_k_mu_arr[i * c_m_d + (j + 1)] + 2 * e_k_mu_arr[i * c_m_d + j]) / c_coef8;
			a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef6 - (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef8;
			a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 - (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef7;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[i * c_m_d + j - 1] / 4 - e_k_mu_arr[(i + 1) * c_m_d + j] / 6) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
			a_arr[a][1] = -2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5;
			a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5
				+ (e_k_mu_arr[(i + 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / (4 * c_hx_d * c_hx_d * C_Re) + (e_k_mu_arr[(i + 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j]) / (8 * c_hx_d * c_hx_d * C_Re)
				+ (e_k_mu_arr[i * c_m_d + (j + 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5 + (e_k_mu_arr[i * c_m_d + (j + 1)] + 2 * e_k_mu_arr[i * c_m_d + j]) / (6 * c_hy_d * c_hy_d * C_Re);
			a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5 - (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / (6 * c_hy_d * c_hy_d * C_Re);
			a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / (4 * c_hx_d * c_hx_d * C_Re) - (2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / (8 * c_hx_d * c_hx_d * C_Re);
			a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j - 1] / 6) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
		}
		//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < c_m_d - 1; j++)
			{
				// u
				int a = i * c_m_d + j;
				a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
				a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3;
				a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
				a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
				a_arr[a][4] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2;
				a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[i * c_m_d + j - 1] / 4 - e_k_mu_arr[(i + 1) * c_m_d + j] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j + 1)] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
				// v
				a += C_M2;
				a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
				a_arr[a][1] = -2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5;
				a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
				a_arr[a][3] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
				a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4;
				a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j - 1] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[i * c_m_d + j + 1] / 6 - e_k_mu_arr[(i + 1) * c_m_d + j] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// u
				int a = i * c_m_d + j;
				a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
				a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3;
				a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
				a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
				a_arr[a][4] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2;
				a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[i * c_m_d + j - 1] / 4 - e_k_mu_arr[(i + 1) * c_m_d + j] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j + 1)] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
				// v
				a += C_M2;
				a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
				a_arr[a][1] = -2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5;
				a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
				a_arr[a][3] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
				a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4;
				a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
				a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
				a_arr[a][7] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j - 1] / 6) / c_coef1;
				a_arr[a][8] = (e_k_mu_arr[i * c_m_d + j + 1] / 6 - e_k_mu_arr[(i + 1) * c_m_d + j] / 4) / c_coef1;
				f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			}
		}
		//#pragma omp single nowait
		{
			// u
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * c_m_d + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
			a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
			a_arr[a][4] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[i * c_m_d + j - 1] / 4 - e_k_mu_arr[(i + 1) * c_m_d + j] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j + 1)] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
			a_arr[a][1] = -2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5;
			a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
			a_arr[a][3] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
			a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4;
			a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j - 1] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[i * c_m_d + j + 1] / 6 - e_k_mu_arr[(i + 1) * c_m_d + j] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * c_m_d + j;
			a_arr[a][0] = -2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef2;
			a_arr[a][1] = -(e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef3;
			a_arr[a][2] = sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d + 2 * (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2 + (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
			a_arr[a][3] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef3;
			a_arr[a][4] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef2;
			a_arr[a][5] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j - 1)] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[i * c_m_d + j + 1] / 4 - e_k_mu_arr[(i - 1) * c_m_d + j] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[i * c_m_d + j - 1] / 4 - e_k_mu_arr[(i + 1) * c_m_d + j] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 6 - e_k_mu_arr[i * c_m_d + (j + 1)] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
			// v
			a += C_M2;
			a_arr[a][0] = -(e_k_mu_arr[(i - 1) * c_m_d + j] + e_k_mu_arr[i * c_m_d + j]) / c_coef4;
			a_arr[a][1] = -2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + e_k_mu_arr[i * c_m_d + j]) / c_coef5;
			a_arr[a][2] = sigma_k1_arr[i * c_m_d + j] * sigma_k1_arr[i * c_m_d + j] / c_tau_d + (e_k_mu_arr[(i - 1) * c_m_d + j] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4 + 2 * (e_k_mu_arr[i * c_m_d + (j - 1)] + 2 * e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
			a_arr[a][3] = -2 * (e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[i * c_m_d + (j + 1)]) / c_coef5;
			a_arr[a][4] = -(e_k_mu_arr[i * c_m_d + j] + e_k_mu_arr[(i + 1) * c_m_d + j]) / c_coef4;
			a_arr[a][5] = (e_k_mu_arr[i * c_m_d + j - 1] / 6 - e_k_mu_arr[(i - 1) * c_m_d + j] / 4) / c_coef1;
			a_arr[a][6] = (e_k_mu_arr[(i - 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j + 1] / 6) / c_coef1;
			a_arr[a][7] = (e_k_mu_arr[(i + 1) * c_m_d + j] / 4 - e_k_mu_arr[i * c_m_d + j - 1] / 6) / c_coef1;
			a_arr[a][8] = (e_k_mu_arr[i * c_m_d + j + 1] / 6 - e_k_mu_arr[(i + 1) * c_m_d + j] / 4) / c_coef1;
			f_arr[a] = (u_kk[a] - trajectory(c_tau_d, c_hx_d, c_hy_d, i, j, u_kk, u_k_arr[a], v_k_arr[a], c_m_d)) * sigma_k1_arr[a] * sigma_k1_arr[a] / c_tau_d - (P(C_gamma, sigma_k_arr[(i + 1) * c_m_d + j], e_arr[(i + 1) * c_m_d + j]) - P(C_gamma, sigma_k_arr[(i - 1) * c_m_d + j], e_arr[(i - 1) * c_m_d + j])) / (2 * c_hx_d);
		}
	} // //#pragma omp parallel
}

//Âåêòîð B = A*Xk1 
inline void mtn_calculate_jakobi(cnst_arr_t u_arr, cnst_arr_t v_arr, cnst_ptr_arr_t u2_arr, cnst_ptr_arr_t v2_arr, cnst_arr_t f_arr, cnst_ptr_2d_arr_t a_arr)
{
	//#pragma omp parallel
	{
		//#pragma omp for nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
#pragma ivdep
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * C_M + j;
					float_type b =

						a_arr[a][0] * u_arr[(i - 1) * C_M + j] +
						a_arr[a][1] * u_arr[i * C_M + (j - 1)] +
						a_arr[a][2] * u_arr[i * C_M + j] +
						a_arr[a][3] * u_arr[i * C_M + (j + 1)] +
						a_arr[a][4] * u_arr[(i + 1) * C_M + j] +

						a_arr[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
						a_arr[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
						a_arr[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
						a_arr[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
					u2_arr[a] = u_arr[a] - 1 / a_arr[a][2] * (b - f_arr[a]);

					// v
					a = i * C_M + j + C_M2;
					b = a_arr[a][0] * v_arr[(i - 1) * C_M + j] +
						a_arr[a][1] * v_arr[i * C_M + j - 1] +
						a_arr[a][2] * v_arr[i * C_M + j] +
						a_arr[a][3] * v_arr[i * C_M + j + 1] +
						a_arr[a][4] * v_arr[(i + 1) * C_M + j] +
						a_arr[a][5] * u_arr[(i - 1) * C_M + j - 1] +
						a_arr[a][6] * u_arr[(i - 1) * C_M + j + 1] +
						a_arr[a][7] * u_arr[(i + 1) * C_M + j - 1] +
						a_arr[a][8] * u_arr[(i + 1) * C_M + j + 1];
					v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / a_arr[a][2] * (b - f_arr[a]);
				}
			}
		}
		//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				// u
				u2_arr[i * C_M + j] =
					u_arr[i * C_M + j] -
					(a_arr[i * C_M + j][0] * u_arr[(i - 1) * C_M + j] +
						a_arr[i * C_M + j][1] * u_arr[i * C_M + (j - 1)] +
						a_arr[i * C_M + j][2] * u_arr[i * C_M + j] +
						a_arr[i * C_M + j][3] * u_arr[i * C_M + (j + 1)] +
						a_arr[i * C_M + j][4] * u_arr[(i + 1) * C_M + j] +
						a_arr[i * C_M + j][5] * v_arr[(i - 1) * C_M + (j - 1)] +
						a_arr[i * C_M + j][6] * v_arr[(i - 1) * C_M + (j + 1)] +
						a_arr[i * C_M + j][7] * v_arr[(i + 1) * C_M + (j - 1)] +
						a_arr[i * C_M + j][8] * v_arr[(i + 1) * C_M + (j + 1)] - f_arr[i * C_M + j]) / a_arr[i * C_M + j][2];
			}
		}
		//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				// v
				v2_arr[i * C_M + j] =
					v_arr[i * C_M + j + C_M2] -
					(a_arr[i * C_M + j + C_M2][0] * v_arr[(i - 1) * C_M + j] +
						a_arr[i * C_M + j + C_M2][1] * v_arr[i * C_M + j - 1] +
						a_arr[i * C_M + j + C_M2][2] * v_arr[i * C_M + j] +
						a_arr[i * C_M + j + C_M2][3] * v_arr[i * C_M + j + 1] +
						a_arr[i * C_M + j + C_M2][4] * v_arr[(i + 1) * C_M + j] +
						a_arr[i * C_M + j + C_M2][5] * u_arr[(i - 1) * C_M + j - 1] +
						a_arr[i * C_M + j + C_M2][6] * u_arr[(i - 1) * C_M + j + 1] +
						a_arr[i * C_M + j + C_M2][7] * u_arr[(i + 1) * C_M + j - 1] +
						a_arr[i * C_M + j + C_M2][8] * u_arr[(i + 1) * C_M + j + 1] - f_arr[i * C_M + j + C_M2]) / a_arr[i * C_M + j + C_M2][2];
			}
		}
		//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// u
				int a = i * C_M + j;
				float_type b = a_arr[a][0] * u_arr[(i - 1) * C_M + j] + a_arr[a][1] * u_arr[i * C_M + (j - 1)] + a_arr[a][2] * u_arr[i * C_M + j] +
					a_arr[a][3] * u_arr[i * C_M + (j + 1)] + a_arr[a][4] * u_arr[(i + 1) * C_M + j] +
					a_arr[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
					a_arr[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
					a_arr[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
					a_arr[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
				u2_arr[a] = u_arr[a] - 1 / a_arr[a][2] * (b - f_arr[a]);
			}
		}
		//#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// v
				int a = i * C_M + j + C_M2;
				float_type b = a_arr[a][0] * v_arr[(i - 1) * C_M + j] + a_arr[a][1] * v_arr[i * C_M + j - 1] + a_arr[a][2] * v_arr[i * C_M + j] +
					a_arr[a][3] * v_arr[i * C_M + j + 1] + a_arr[a][4] * v_arr[(i + 1) * C_M + j] +
					a_arr[a][5] * u_arr[(i - 1) * C_M + j - 1] +
					a_arr[a][6] * u_arr[(i - 1) * C_M + j + 1] +
					a_arr[a][7] * u_arr[(i + 1) * C_M + j - 1] +
					a_arr[a][8] * u_arr[(i + 1) * C_M + j + 1];
				v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / a_arr[a][2] * (b - f_arr[a]);
			}
		}
		//#pragma omp single nowait
		{
			// u
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			float_type b = a_arr[a][0] * u_arr[(i - 1) * C_M + j] + a_arr[a][1] * u_arr[i * C_M + (j - 1)] + a_arr[a][2] * u_arr[i * C_M + j] +
				a_arr[a][3] * u_arr[i * C_M + (j + 1)] + a_arr[a][4] * u_arr[(i + 1) * C_M + j] +
				a_arr[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
				a_arr[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
				a_arr[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
				a_arr[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
			u2_arr[a] = u_arr[a] - 1 / a_arr[a][2] * (b - f_arr[a]);
			// v
			a += C_M2;
			b = a_arr[a][0] * v_arr[(i - 1) * C_M + j] + a_arr[a][1] * v_arr[i * C_M + j - 1] + a_arr[a][2] * v_arr[i * C_M + j] +
				a_arr[a][3] * v_arr[i * C_M + j + 1] + a_arr[a][4] * v_arr[(i + 1) * C_M + j] +
				a_arr[a][5] * u_arr[(i - 1) * C_M + j - 1] +
				a_arr[a][6] * u_arr[(i - 1) * C_M + j + 1] +
				a_arr[a][7] * u_arr[(i + 1) * C_M + j - 1] +
				a_arr[a][8] * u_arr[(i + 1) * C_M + j + 1];
			v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / a_arr[a][2] * (b - f_arr[a]);
			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			b = a_arr[a][0] * u_arr[(i - 1) * C_M + j] + a_arr[a][1] * u_arr[i * C_M + (j - 1)] + a_arr[a][2] * u_arr[i * C_M + j] +
				a_arr[a][3] * u_arr[i * C_M + (j + 1)] + a_arr[a][4] * u_arr[(i + 1) * C_M + j] +
				a_arr[a][5] * v_arr[(i - 1) * C_M + (j - 1)] +
				a_arr[a][6] * v_arr[(i - 1) * C_M + (j + 1)] +
				a_arr[a][7] * v_arr[(i + 1) * C_M + (j - 1)] +
				a_arr[a][8] * v_arr[(i + 1) * C_M + (j + 1)];
			u2_arr[a] = u_arr[a] - 1 / a_arr[a][2] * (b - f_arr[a]);
			// v
			a += C_M2;
			b = a_arr[a][0] * v_arr[(i - 1) * C_M + j] + a_arr[a][1] * v_arr[i * C_M + j - 1] + a_arr[a][2] * v_arr[i * C_M + j] +
				a_arr[a][3] * v_arr[i * C_M + j + 1] + a_arr[a][4] * v_arr[(i + 1) * C_M + j] +
				a_arr[a][5] * u_arr[(i - 1) * C_M + j - 1] +
				a_arr[a][6] * u_arr[(i - 1) * C_M + j + 1] +
				a_arr[a][7] * u_arr[(i + 1) * C_M + j - 1] +
				a_arr[a][8] * u_arr[(i + 1) * C_M + j + 1];
			v2_arr[a - C_M2] = v_arr[a - C_M2] - 1 / a_arr[a][2] * (b - f_arr[a]);
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
	mtn_calculate_common(A, sigma_k_arr, sigma_k1_arr, e_k_arr, e_k_mu_arr, u_k_arr, v_k_arr, f_arr, C_tau, C_hx, C_hy, C_M);

	int c_u;
	int c_v;
	int s_m = 0;
	for (s_m = 0; s_m <= 20; ++s_m)
	{
		mtn_calculate_jakobi(u_k1_arr, v_k1_arr, u2_arr, v2_arr, f_arr, A);

		c_u = 0;
		c_v = 0;

		//#pragma omp parallel for reduction(+:c_u, c_v)
#pragma vector aligned
		for (int i = 0; i < C_M2; i++)
		{
			if (fabs(u_k1_arr[i] - u2_arr[i]) <= C_epsilon) ++c_u;
			if (fabs(v_k1_arr[i] - v2_arr[i]) <= C_epsilon) ++c_v;
		}

		if (c_u - C_cnt_boundary_nodes >= C_br && c_v - C_cnt_boundary_nodes >= C_br)
		{
			s_m++;
			return s_m;
		}

		// Ìîæíî êîïèðîâàòü ñ memcpy?
		//#pragma omp parallel for 
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

inline int interate_over_nonlinearity(int* s_m, int* s_e, int* s_end)
{
	const int itr = 5;

	int s_itr;
	for (s_itr = 1; s_itr < itr; s_itr++)
	{
#pragma vector aligned
		for (int i = 0; i < C_M2; i++)
			e_k_mu[i] = Mu(C_gamma_Mah2, e_k[i]);

		continuity(sigma_kk, sigma_k1, u_k, v_k);
		*s_m = motion(sigma_k, sigma_k1, u_k, v_k, f, u_k1, v_k1, u2, v2, e_k, e_k_mu);
		*s_e = energy(sigma_k, sigma_k1, e_k, e_k1, e2, e_k_mu, u_k, v_k, f);

		if (*s_m == 1 && *s_e == 1)
		{
			*s_end = s_itr;
			s_itr = itr;
		}

		memcpy(sigma_k, sigma_k1, C_M * sizeof *sigma_k);
		memcpy(e_k, e_k1, C_M * sizeof *e_k);
		memcpy(u_k, u_k1, C_M * sizeof *u_k);
		memcpy(v_k, v_k1, C_M * sizeof *v_k);

		//#pragma omp parallel for 
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
#pragma ivdep
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
	A = (float_type**)_mm_malloc(double_size_array*sizeof(float_type*), ALIGN);
	for (int i = 0; i < double_size_array; ++i)
	{
		A[i] = (float_type*)_mm_malloc(param_array_element_count*sizeof(float_type*), ALIGN);
		memset(A[i], 0., param_array_element_count * sizeof(float_type*));
	}

	f = (float_type*)_mm_malloc(double_size_array*sizeof(float_type*), ALIGN);
	sigma_k = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	sigma_k1 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	sigma_kk = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	u_k = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	u_k1 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	u_kk = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	u2 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	v_k = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	v_k1 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	v2 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	e_k = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	e_k1 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	e_kk = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	e2 = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);
	e_k_mu = (float_type*)_mm_malloc(array_element_count*sizeof(float_type*), ALIGN);

	memset(f, 0., double_size_array * sizeof(float_type*));
	memset(sigma_k, 0., array_element_count * sizeof(float_type*));
	memset(sigma_k1, 0., array_element_count * sizeof(float_type*));
	memset(sigma_kk, 0., array_element_count * sizeof(float_type*));
	memset(u_k, 0., array_element_count * sizeof(float_type*));
	memset(u_k1, 0., array_element_count * sizeof(float_type*));
	memset(u_kk, 0., array_element_count * sizeof(float_type*));
	memset(u2, 0., array_element_count * sizeof(float_type*));
	memset(v_k, 0., array_element_count * sizeof(float_type*));
	memset(v_k1, 0., array_element_count * sizeof(float_type*));
	memset(v2, 0., array_element_count * sizeof(float_type*));
	memset(e_k, 0., array_element_count * sizeof(float_type*));
	memset(e_k1, 0., array_element_count * sizeof(float_type*));
	memset(e_kk, 0., array_element_count * sizeof(float_type*));
	memset(e2, 0., array_element_count * sizeof(float_type*));
	memset(e_k_mu, 0., array_element_count * sizeof(float_type*));
}

void clear_memory_parallel(const int array_element_count)
{
	for (int i = 0; i < 2 * array_element_count; i++) _mm_free(A[i]);
	_mm_free(A);
	_mm_free(f);
	_mm_free(sigma_k);
	_mm_free(sigma_kk);
	_mm_free(sigma_k1);
	_mm_free(u_k);
	_mm_free(u_k1);
	_mm_free(u_kk);
	_mm_free(u2);
	_mm_free(v_k);
	_mm_free(v_k1);
	_mm_free(v2);
	_mm_free(e_k);
	_mm_free(e_k1);
	_mm_free(e_kk);
	_mm_free(e2);
	_mm_free(e_k_mu);
}

int get_length_parallel()
{
	return C_M2;
}

int get_length_parallel_x()
{
	return C_M1;
}

int get_length_parallel_y()
{
	return C_M;
}

double* get_sigma_parallel()
{
	return sigma_k1;
}

double* get_u_parallel()
{
	return u_k1;
}

double* get_v_parallel()
{
	return v_k1;
}

double* get_e_parallel()
{
	return e_k1;
}

// ReSharper disable once CppParameterNeverUsed
double calculate_parallel(const bool need_print, const int thread_count)
{
	FILE* fout = NULL;
	FILE* fout_itr = NULL;
	FILE* fdensity = NULL;
	FILE* fdensity_new = NULL;
	FILE* fvelocity = NULL;
	FILE* ftemperature = NULL;
	FILE* fpressure = NULL;
	if (need_print)
	{
		fout = fopen("out_p.txt", "w");
		fout_itr = fopen("out_itr_p.txt", "w");
		fdensity = fopen("density_p.dat", "w");
		fdensity_new = fopen("density-new_p.dat", "w");
		fvelocity = fopen("velocity_p.dat", "w");
		ftemperature = fopen("temperature_p.dat", "w");
		fpressure = fopen("pressure_p.dat", "w");
		print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr, C_tau, C_hx, C_hy, C_N);
	}
	double time;

	init_arrays_parallel(C_M2, 9);
	set_initial_boundary_conditions_parallel();

#ifdef _OPENMP
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(thread_count);
	//omp_set_num_threads(omp_get_num_procs());
	printf("OPENMP THREADS COUNT = %d\n", omp_get_max_threads());
	long count = 0;
	// dummy parallel section to get all threads running
	//#pragma omp parallel
	{
		_InterlockedIncrement(&count);
	}
	//printf("OPENMP timer function is used!\n");
	time = omp_get_wtime();
#else
	printf("Standart timer function is used!\n");
	StartTimer();
#endif
	for (int current_time_step = 1; current_time_step <= time_steps_nbr; current_time_step++)
	{
		int s_end = 0;
		int s_m = 0;
		int s_e = 0;
		int s_itr;
		printf("Par TS = %d\n", current_time_step);
		memcpy(sigma_k, sigma_k1, C_M2 * sizeof *sigma_k);
		memcpy(sigma_kk, sigma_k1, C_M2 * sizeof *sigma_kk);
		memcpy(e_k, e_k1, C_M2 * sizeof *e_k);
		memcpy(e_kk, e_k1, C_M2 * sizeof *e_kk);
		memcpy(u_k, u_k1, C_M2 * sizeof *u_k);
		memcpy(u_kk, u_k1, C_M2 * sizeof *u_kk);
		memcpy(v_k, v_k1, C_M2 * sizeof *v_k);
		s_itr = interate_over_nonlinearity(&s_m, &s_e, &s_end);
		if (need_print)
			print_to_file(C_gamma, s_m, s_e, current_time_step, s_itr, s_end, C_tau, C_hx, C_hy, C_M, C_M1, C_N, C_Mah2, sigma_k1, u_k1, v_k1, e_k1, C_gamma_Mah2,
			              C_q, C_w, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
	}

#ifdef _OPENMP 
	time = omp_get_wtime() - time;
#else
	time = GetTimer() / 1000;
#endif

	if (need_print)
	{
		print_new_line(fout, fdensity, fvelocity, ftemperature, fpressure);
		fclose(fout);
		fclose(fdensity);
		fclose(fvelocity);
		fclose(ftemperature);
		fclose(fpressure);
		fclose(fout_itr);
	}
	return time;
}

