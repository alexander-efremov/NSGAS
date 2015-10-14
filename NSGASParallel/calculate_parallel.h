#include <stdio.h>
#include <math.h>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 1
#define omp_get_num_threads() 1
#endif

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
//static const double C_hx = 1.0 / C_N1;
//static const double C_hy = 1.0 / C_N;
//static const int time_steps_nbr = 2; // time_steps_nbr - количество шагов по времени
//-----------------------
// real test 
static const int C_N = 1200;
static const int C_N1 = 800;
static const int C_q = 101;
static const int C_qq = 20;
static const double C_hx = 1.0 / 100;
static const double C_hy = 1.0 / 200;
//static const int time_steps_nbr = 25000; // time_steps_nbr - количество шагов по времени
static const int time_steps_nbr = 1; // time_steps_nbr - количество шагов по времени
//-----------------------
static const int C_M = C_N + 1;
static const int C_M1 = C_N1 + 1;
static const int C_M2 = C_M1 * C_M;
static const int C_w = C_q;
static const int C_cntr = C_N / 2;
static const double C_tau = 0.0005;
static const double C_tg = 2;
static const double C_Re = 10000;
static const double C_PrRe = 0.72 * C_Re; // Pr * C_Re
static const double C_Mah2 = 16; // Mah * Mah
static const double C_epsilon = 0.0000000001;
// В массивах с _k1 хранятся значения функций на d-ом шаге по времени
// В массивах с _k хранятся значения функций с предыдущей итерации по нелинейности
// В массивах с _kk хранятся значения функций c (d-1) шага по времени
// Массивы с "2" использутся в итерациях метода Зейделя
// В массивах с X_k хранятся значения функций, вычисленных методом траекторий
static double** A;
static double* B;
static double* D;
static double* f;
static double* sigma_k;
static double* sigma_k1;
static double* u_k;
static double* u_k1;
static double* v_k;
static double* v_k1;
static double* u2;
static double* v2;
static double* e_k;
static double* e_k1;
static double* e2;
static double* T;
static double* sigma_kk;
static double* u_kk;
static double* v_kk;
static double* e_kk;
static double* sigmaX_k;
static double* uX_k;
static double* vY_k;
static double* eR_k;

inline double Mu(double gamma, double e_k)
{
	const double omega = 0.8;
	return pow(gamma * (gamma - 1) * C_Mah2 * e_k * e_k, omega);
}

inline double P(double gamma, double sigma_k, double e_k)
{
	return (gamma - 1) * sigma_k * sigma_k * e_k * e_k;
}

#include "trajectory_cone.h"
#include "continuity_sigma.h"
#include "motion.h"
#include "energy_epsilon.h"

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

// m = C_M
// m = C_M1
// n = C_N
// mah2 = C_Mah2
inline void print_to_file(const double gamma, int s_m, int s_e,
                          int current_ts, int s_itr, int s_end,
                          const double tau, const double hx,
                          const double hy,
                          const int m, const int m1, const int n,
                          const double mah2,
                          FILE* out, FILE* density, FILE* density_new, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	int i;
	int j;
	int a;
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
		double ts_tau = current_ts * tau;

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

		for (j = 0; j < m; j++)
		{
			for (i = 0; i < m1; i++)
			{
				a = i * m + j;
				double ihx = i * hx;
				double jhy = j * hy;
				fprintf(density_new, "%.3f\t %.4f\t %.15e\n", ihx, jhy, sigma_k1[a] * sigma_k1[a]);

				fprintf(density, "%.3f\t %.4f\t %.10f\n", ihx, jhy, sigma_k1[a] * sigma_k1[a]);
				fprintf(velocity, "%.3f\t %.4f\t%.10f\t %.10f\n", ihx, jhy, u_k1[a], v_k1[a]);
				fprintf(temperature, "%.3f\t %.4f\t %.10f\n", ihx, jhy, e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * mah2));
				fprintf(pressure, "%.3f\t %.4f\t %.10f\n", ihx, jhy, sigma_k1[a] * sigma_k1[a] * e_k1[a] * e_k1[a] * (gamma - 1));
			}
		}

		if (current_ts == 1)
		{
			for (i = 0; i < m1; i++)
			{
				for (j = 0; j < m; j++)
				{
					a = i * m + j;
					fprintf(out, "i=%i j=%i\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\n", i, j,
					        sigma_k1[a] * sigma_k1[a], u_k1[a], v_k1[a], e_k1[a] * e_k1[a], e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * mah2), P(gamma, sigma_k1[a], e_k1[a]), Mu(gamma, e_k1[a]));
				}
			}
		}
		print_new_line(out, density, velocity, temperature, pressure);
		flush_file(out, density, velocity, temperature, pressure, out_itr);
	}
}

// n = C_N
inline void print_file_header(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr, const double tau, const double hx, const double hy, const int n)
{
	fprintf(out, "Cone_2D\n\n");
	fprintf(out, "C_N = %i\t C_hx = %.5f\t C_hy = %.5f\t C_tau = %.5f\n\n", n, hx, hy, tau);
	fprintf(out_itr, "Cone_2D\n\n");
	fprintf(out_itr, "C_N = %i\t C_hx = %.5f\t C_hy = %.5f\t C_tau = %.5f\n\n", n, hx, hy, tau);
	fprintf(density, "TITLE=\"density\"\n\nVARIABLES=\"x\",\"y\",\"Ro\"\n\n");
	fprintf(velocity, "TITLE=\"velocity\"\n\nVARIABLES=\"x\",\"y\",\"u\",\"v\"\n\n");
	fprintf(temperature, "TITLE=\"temperature\"\n\nVARIABLES=\"x\",\"y\",\"T\"\n\n");
	fprintf(pressure, "TITLE=\"pressure\"\n\nVARIABLES=\"x\",\"y\",\"P\"\n\n");
}

inline void close_files(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	fclose(out);
	fclose(density);
	fclose(velocity);
	fclose(temperature);
	fclose(pressure);
	fclose(out_itr);
}

// Initial boundary conditions with t = 0
// qq_i = C_q
// w_i = C_w
// m_i = C_M
// m1_i = C_M1
// m2_i = C_M2
// mah2 = C_Mah2
// cntr_i = C_cntr
inline void set_initial_boundary_conditions(const double gamma, const int qq_i, const int w_i, const int m_i, const int m1_i, const int m2_i, const double mah2, const int cntr_i)
{
	int a;
	int i;
	int j;

	for (i = 0; i < qq_i; i++)
	{
		for (j = 0; j < m_i; j++)
		{
			a = i * m_i + j;
			if (i == 0)
			{
				sigma_k1[a] = 1;
				T[a] = 1;
				e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * mah2));
				u_k1[a] = 1;
				v_k1[a] = 0;
				e2[a] = e_k1[a];
				u2[a] = u_k1[a];
				v2[a] = v_k1[a];
			}
			if (i > 0 && i < qq_i)
			{
				sigma_k1[a] = 1;
				T[a] = 1;
				e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * mah2));
				u_k1[a] = 0;
				v_k1[a] = 0;
				e2[a] = e_k1[a];
				u2[a] = u_k1[a];
				v2[a] = v_k1[a];
			}
		}
	}

	for (i = qq_i; i < qq_i + w_i - 1; i++)
	{
		for (j = cntr_i + i - qq_i; j < m_i; j++)
		{
			a = i * m_i + j;
			sigma_k1[a] = 1;
			T[a] = 1;
			e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * mah2));
			u_k1[a] = 0;
			v_k1[a] = 0;
			e2[a] = e_k1[a];
			u2[a] = u_k1[a];
			v2[a] = v_k1[a];
		}
	}

	for (i = qq_i; i < qq_i + w_i - 1; i++)
	{
		for (j = cntr_i - i + qq_i; j > -1; j--)
		{
			a = i * m_i + j;
			sigma_k1[a] = 1;
			T[a] = 1;
			e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * mah2));
			u_k1[a] = 0;
			v_k1[a] = 0;
			e2[a] = e_k1[a];
			u2[a] = u_k1[a];
			v2[a] = v_k1[a];
		}
	}

	for (i = qq_i + w_i - 1; i < m1_i; i++)
	{
		for (j = 0; j < m_i; j++)
		{
			a = i * m_i + j;
			sigma_k1[a] = 1;
			T[a] = 1;
			e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * mah2));
			u_k1[a] = 0;
			v_k1[a] = 0;
			e2[a] = e_k1[a];
			u2[a] = u_k1[a];
			v2[a] = v_k1[a];
		}
	}

	for (i = 0; i < m2_i; i++)
	{
		sigmaX_k[i] = 0;
		uX_k[i] = 0;
		vY_k[i] = 0;
		eR_k[i] = 0;
	}
}


// qq_i = C_qq
// m_i = C_M
// n = C_N
// m1_i = C_M1
// m2_i = C_M2
// n1_i = C_N1
// C_w = C_w
// cntr_i = C_cntr
// q_i = C_q
inline int interate_over_nonlinearity(const double gamma,
                                      const int qq_i,
                                      const int m_i, 
									  const int m1_i,
									  const int m2_i,
									  const int n1_i,
                                      const int w_i, const int cntr_i,
                                      const int n, const int q_i, int& s_m, int& s_e, int& s_end)
{
	const int itr = 5;
	int i;
	int j;

	int s_itr;
	for (s_itr = 1; s_itr < itr; ++s_itr)
	{
#pragma omp parallel 
		{
#pragma omp for collapse(2) private(i, j) nowait
			for (i = 1; i < qq_i + 1; i++)
			{
				for (j = 1; j < m_i - 1; j++)
				{
					sigmaX_k[i * m_i + j] = sigma_kk[i * m_i + j] - trajectory(i, j, sigma_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					uX_k[i * m_i + j] = u_kk[i * m_i + j] - trajectory(i, j, u_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					vY_k[i * m_i + j] = v_kk[i * m_i + j] - trajectory(i, j, v_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					eR_k[i * m_i + j] = e_kk[i * m_i + j] - trajectory(i, j, e_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i; i < qq_i + w_i - 1; i++)
			{
				for (j = cntr_i + i - qq_i; j < m_i - 1; j++)
				{
					sigmaX_k[i * m_i + j] = sigma_kk[i * m_i + j] - trajectory(i, j, sigma_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					uX_k[i * m_i + j] = u_kk[i * m_i + j] - trajectory(i, j, u_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					vY_k[i * m_i + j] = v_kk[i * m_i + j] - trajectory(i, j, v_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					eR_k[i * m_i + j] = e_kk[i * m_i + j] - trajectory(i, j, e_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i; i < qq_i + w_i - 1; i++)
			{
				for (j = cntr_i - i + qq_i; j > 0; j--)
				{
					sigmaX_k[i * m_i + j] = sigma_kk[i * m_i + j] - trajectory(i, j, sigma_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					uX_k[i * m_i + j] = u_kk[i * m_i + j] - trajectory(i, j, u_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					vY_k[i * m_i + j] = v_kk[i * m_i + j] - trajectory(i, j, v_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					eR_k[i * m_i + j] = e_kk[i * m_i + j] - trajectory(i, j, e_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
				}
			}

#pragma omp for collapse(2) private(i, j) nowait
			for (i = qq_i + w_i - 1; i < m1_i - 1; i++)
			{
				for (j = 1; j < m_i - 1; j++)
				{
					sigmaX_k[i * m_i + j] = sigma_kk[i * m_i + j] - trajectory(i, j, sigma_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					uX_k[i * m_i + j] = u_kk[i * m_i + j] - trajectory(i, j, u_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					vY_k[i * m_i + j] = v_kk[i * m_i + j] - trajectory(i, j, v_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
					eR_k[i * m_i + j] = e_kk[i * m_i + j] - trajectory(i, j, e_kk, u_k[i * m_i + j], v_k[i * m_i + j], m_i);
				}
			}
		} // #pragma omp parallel

		continuity(sigma_k1, u_k, v_k, qq_i, w_i, m_i, cntr_i, m1_i, q_i, C_tau, C_hx, C_hy);
		s_m = motion(gamma, m_i, m1_i, m2_i, qq_i, w_i, cntr_i, sigma_k1, sigma_k, u_k, v_k, u_k1, v_k1, u2, v2, e_k);
		s_e = energy(gamma, sigma_k1, sigma_k, u_k, v_k, u_k1, v_k1, e2, e_k, e_k1, m_i, n, qq_i, w_i, m1_i, q_i, cntr_i);

		if (s_m == 1 && s_e == 1)
		{
			s_end = s_itr;
			s_itr = itr;
		}

#pragma omp parallel for private(j)
		for (j = 0; j < m_i; j++)
		{
			sigma_k[j] = sigma_k1[j];
			e_k[j] = e_k1[j];
			u_k[j] = u_k1[j];
			v_k[j] = v_k1[j];
		}

#pragma omp parallel 
		{
#pragma omp for collapse(2) private(i, j) nowait
			for (i = 1; i < qq_i + 1; i++)
			{
				for (j = 0; j < m_i; j++)
				{
					if (j == 0)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j + 1];
						sigma_k1[i * m_i + j] = sigma_k1[i * m_i + j + 1];
						e_k[i * m_i + j] = e_k1[i * m_i + j + 1];
						e_k1[i * m_i + j] = e_k1[i * m_i + j + 1];
						u_k[i * m_i + j] = u_k1[i * m_i + j + 1];
						u_k1[i * m_i + j] = u_k1[i * m_i + j + 1];
						v_k[i * m_i + j] = v_k1[i * m_i + j + 1];
						v_k1[i * m_i + j] = v_k1[i * m_i + j + 1];
						e2[i * m_i + j] = e_k1[i * m_i + j + 1];
						u2[i * m_i + j] = u_k1[i * m_i + j + 1];
						v2[i * m_i + j] = v_k1[i * m_i + j + 1];
					}
					if (j > 0 && j < n)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j];
						e_k[i * m_i + j] = e_k1[i * m_i + j];
						u_k[i * m_i + j] = u_k1[i * m_i + j];
						v_k[i * m_i + j] = v_k1[i * m_i + j];
					}
					if (j == n)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j - 1];
						sigma_k1[i * m_i + j] = sigma_k1[i * m_i + j - 1];
						e_k[i * m_i + j] = e_k1[i * m_i + j - 1];
						e_k1[i * m_i + j] = e_k1[i * m_i + j - 1];
						u_k[i * m_i + j] = u_k1[i * m_i + j - 1];
						u_k1[i * m_i + j] = u_k1[i * m_i + j - 1];
						v_k[i * m_i + j] = v_k1[i * m_i + j - 1];
						v_k1[i * m_i + j] = v_k1[i * m_i + j - 1];
						e2[i * m_i + j] = e_k1[i * m_i + j - 1];
						u2[i * m_i + j] = u_k1[i * m_i + j - 1];
						v2[i * m_i + j] = v_k1[i * m_i + j - 1];
					}
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i; i < qq_i + w_i - 1; i++)
			{
				for (j = cntr_i + i - qq_i; j < m_i; j++)
				{
					if (j == n)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j - 1];
						sigma_k1[i * m_i + j] = sigma_k1[i * m_i + j - 1];
						e_k[i * m_i + j] = e_k1[i * m_i + j - 1];
						e_k1[i * m_i + j] = e_k1[i * m_i + j - 1];
						u_k[i * m_i + j] = u_k1[i * m_i + j - 1];
						u_k1[i * m_i + j] = u_k1[i * m_i + j - 1];
						v_k[i * m_i + j] = v_k1[i * m_i + j - 1];
						v_k1[i * m_i + j] = v_k1[i * m_i + j - 1];
						e2[i * m_i + j] = e_k1[i * m_i + j - 1];
						u2[i * m_i + j] = u_k1[i * m_i + j - 1];
						v2[i * m_i + j] = v_k1[i * m_i + j - 1];
					}
					else
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j];
						e_k[i * m_i + j] = e_k1[i * m_i + j];
						u_k[i * m_i + j] = u_k1[i * m_i + j];
						v_k[i * m_i + j] = v_k1[i * m_i + j];
					}
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i; i < qq_i + w_i - 1; i++)
			{
				for (j = cntr_i - i + qq_i; j > -1; j--)
				{
					if (j == 0)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j + 1];
						sigma_k1[i * m_i + j] = sigma_k1[i * m_i + j + 1];
						e_k[i * m_i + j] = e_k1[i * m_i + j + 1];
						e_k1[i * m_i + j] = e_k1[i * m_i + j + 1];
						u_k[i * m_i + j] = u_k1[i * m_i + j + 1];
						u_k1[i * m_i + j] = u_k1[i * m_i + j + 1];
						v_k[i * m_i + j] = v_k1[i * m_i + j + 1];
						v_k1[i * m_i + j] = v_k1[i * m_i + j + 1];
						e2[i * m_i + j] = e_k1[i * m_i + j + 1];
						u2[i * m_i + j] = u_k1[i * m_i + j + 1];
						v2[i * m_i + j] = v_k1[i * m_i + j + 1];
					}
					else
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j];
						e_k[i * m_i + j] = e_k1[i * m_i + j];
						u_k[i * m_i + j] = u_k1[i * m_i + j];
						v_k[i * m_i + j] = v_k1[i * m_i + j];
					}
				}
			}

#pragma omp for collapse(2) private(i, j) nowait
			for (i = qq_i + w_i - 1; i < m1_i - 1; i++)
			{
				for (j = 0; j < m_i; j++)
				{
					if (j == 0)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j + 1];
						sigma_k1[i * m_i + j] = sigma_k1[i * m_i + j + 1];
						e_k[i * m_i + j] = e_k1[i * m_i + j + 1];
						e_k1[i * m_i + j] = e_k1[i * m_i + j + 1];
						u_k[i * m_i + j] = u_k1[i * m_i + j + 1];
						u_k1[i * m_i + j] = u_k1[i * m_i + j + 1];
						v_k[i * m_i + j] = v_k1[i * m_i + j + 1];
						v_k1[i * m_i + j] = v_k1[i * m_i + j + 1];
						e2[i * m_i + j] = e_k1[i * m_i + j + 1];
						u2[i * m_i + j] = u_k1[i * m_i + j + 1];
						v2[i * m_i + j] = v_k1[i * m_i + j + 1];
					}
					if (j > 0 && j < n)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j];
						e_k[i * m_i + j] = e_k1[i * m_i + j];
						u_k[i * m_i + j] = u_k1[i * m_i + j];
						v_k[i * m_i + j] = v_k1[i * m_i + j];
					}
					if (j == n)
					{
						sigma_k[i * m_i + j] = sigma_k1[i * m_i + j - 1];
						sigma_k1[i * m_i + j] = sigma_k1[i * m_i + j - 1];
						e_k[i * m_i + j] = e_k1[i * m_i + j - 1];
						e_k1[i * m_i + j] = e_k1[i * m_i + j - 1];
						u_k[i * m_i + j] = u_k1[i * m_i + j - 1];
						u_k1[i * m_i + j] = u_k1[i * m_i + j - 1];
						v_k[i * m_i + j] = v_k1[i * m_i + j - 1];
						v_k1[i * m_i + j] = v_k1[i * m_i + j - 1];
						e2[i * m_i + j] = e_k1[i * m_i + j - 1];
						u2[i * m_i + j] = u_k1[i * m_i + j - 1];
						v2[i * m_i + j] = v_k1[i * m_i + j - 1];
					}
				}
			}
#pragma omp for private(i, j) nowait
			for (j = 0; j < m_i; j++)
			{				
				int indx = (n1_i - 1) * m_i + j;
				if (j == 0)
				{
					sigma_k[n1_i * m_i + j] = sigma_k1[indx + 1];
					sigma_k1[n1_i * m_i + j] = sigma_k1[indx + 1];
					e_k[n1_i * m_i + j] = e_k1[indx + 1];
					e_k1[n1_i * m_i + j] = e_k1[indx + 1];
					u_k[n1_i * m_i + j] = u_k1[indx + 1];
					u_k1[n1_i * m_i + j] = u_k1[indx + 1];
					v_k[n1_i * m_i + j] = v_k1[indx + 1];
					v_k1[n1_i * m_i + j] = v_k1[indx + 1];
					e2[n1_i * m_i + j] = e_k1[indx + 1];
					u2[n1_i * m_i + j] = u_k1[indx + 1];
					v2[n1_i * m_i + j] = v_k1[indx + 1];
				}
				if (j > 0 && j < n)
				{
					sigma_k[n1_i * m_i + j] = sigma_k1[indx];
					sigma_k1[n1_i * m_i + j] = sigma_k1[indx];
					e_k[n1_i * m_i + j] = e_k1[indx];
					e_k1[n1_i * m_i + j] = e_k1[indx];
					u_k[n1_i * m_i + j] = u_k1[indx];
					u_k1[n1_i * m_i + j] = u_k1[indx];
					v_k[n1_i * m_i + j] = v_k1[indx];
					v_k1[n1_i * m_i + j] = v_k1[indx];
					e2[n1_i * m_i + j] = e_k1[indx];
					u2[n1_i * m_i + j] = u_k1[indx];
					v2[n1_i * m_i + j] = v_k1[indx];
				}
				if (j == n)
				{
					sigma_k[n1_i * m_i + j] = sigma_k1[indx - 1];
					sigma_k1[n1_i * m_i + j] = sigma_k1[indx - 1];
					e_k[n1_i * m_i + j] = e_k1[indx - 1];
					e_k1[n1_i * m_i + j] = e_k1[indx - 1];
					u_k[n1_i * m_i + j] = u_k1[indx - 1];
					u_k1[n1_i * m_i + j] = u_k1[indx - 1];
					v_k[n1_i * m_i + j] = v_k1[indx - 1];
					v_k1[n1_i * m_i + j] = v_k1[indx - 1];
					e2[n1_i * m_i + j] = e_k1[indx - 1];
					u2[n1_i * m_i + j] = u_k1[indx - 1];
					v2[n1_i * m_i + j] = v_k1[indx - 1];
				}
			}
		} // #pragma omp parallel		
	}
	return s_itr;
}

// m = C_M
// m1 = C_M1
// qq_i = C_qq
// w_i = C_w
// cntr_i = C_cntr
// Should I try memcpy instead of for loop?
inline void prepare_to_iterate(const int m, const int m1, const int qq_i, const int w_i, const int cntr_i)
{
	int i;
	int j;
#pragma omp parallel
	{
#pragma omp for collapse(2) private(i, j) nowait
		for (i = 0; i < qq_i + 1; i++)
		{
			for (j = 0; j < m; j++)
			{
				sigma_k[i * m + j] = sigma_k1[i * m + j];
				e_k[i * m + j] = e_k1[i * m + j];
				u_k[i * m + j] = u_k1[i * m + j];
				v_k[i * m + j] = v_k1[i * m + j];
				sigma_kk[i * m + j] = sigma_k1[i * m + j];
				u_kk[i * m + j] = u_k1[i * m + j];
				v_kk[i * m + j] = v_k1[i * m + j];
				e_kk[i * m + j] = e_k1[i * m + j];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i; i < qq_i + w_i - 1; i++)
		{
			for (j = cntr_i + i - qq_i; j < m; j++)
			{
				sigma_k[i * m + j] = sigma_k1[i * m + j];
				e_k[i * m + j] = e_k1[i * m + j];
				u_k[i * m + j] = u_k1[i * m + j];
				v_k[i * m + j] = v_k1[i * m + j];
				sigma_kk[i * m + j] = sigma_k1[i * m + j];
				u_kk[i * m + j] = u_k1[i * m + j];
				v_kk[i * m + j] = v_k1[i * m + j];
				e_kk[i * m + j] = e_k1[i * m + j];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i; i < qq_i + w_i - 1; i++)
		{
			for (j = cntr_i - i + qq_i; j > -1; j--)
			{
				sigma_k[i * m + j] = sigma_k1[i * m + j];
				e_k[i * m + j] = e_k1[i * m + j];
				u_k[i * m + j] = u_k1[i * m + j];
				v_k[i * m + j] = v_k1[i * m + j];
				sigma_kk[i * m + j] = sigma_k1[i * m + j];
				u_kk[i * m + j] = u_k1[i * m + j];
				v_kk[i * m + j] = v_k1[i * m + j];
				e_kk[i * m + j] = e_k1[i * m + j];
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = qq_i + w_i - 1; i < m1; i++)
		{
			for (j = 0; j < m; j++)
			{
				sigma_k[i * m + j] = sigma_k1[i * m + j];
				e_k[i * m + j] = e_k1[i * m + j];
				u_k[i * m + j] = u_k1[i * m + j];
				v_k[i * m + j] = v_k1[i * m + j];
				sigma_kk[i * m + j] = sigma_k1[i * m + j];
				u_kk[i * m + j] = u_k1[i * m + j];
				v_kk[i * m + j] = v_k1[i * m + j];
				e_kk[i * m + j] = e_k1[i * m + j];
			}
		}
	} // #pragma omp parallel
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

	sigma_k = new double[array_element_count];
	u_k = new double[array_element_count];
	v_k = new double[array_element_count];

	sigma_k1 = new double[array_element_count];
	u_k1 = new double[array_element_count];
	v_k1 = new double[array_element_count];

	u2 = new double[array_element_count];
	v2 = new double[array_element_count];

	e_k = new double[array_element_count];
	e_k1 = new double[array_element_count];
	e2 = new double[array_element_count];
	T = new double[array_element_count];
	sigma_kk = new double[array_element_count];
	u_kk = new double[array_element_count];
	v_kk = new double[array_element_count];
	e_kk = new double[array_element_count];

	sigmaX_k = new double[array_element_count];
	uX_k = new double[array_element_count];
	vY_k = new double[array_element_count];
	eR_k = new double[array_element_count];

	std::fill_n(B, double_size_array, 0.);
	std::fill_n(D, double_size_array, 0.);
	std::fill_n(f, double_size_array, 0.);

	std::fill_n(sigma_k, array_element_count, 0.);
	std::fill_n(u_k, array_element_count, 0.);
	std::fill_n(v_k, array_element_count, 0.);

	std::fill_n(sigma_k1, array_element_count, 0.);
	std::fill_n(u_k1, array_element_count, 0.);
	std::fill_n(v_k1, array_element_count, 0.);

	std::fill_n(u2, array_element_count, 0.);
	std::fill_n(v2, array_element_count, 0.);

	std::fill_n(e_k, array_element_count, 0.);
	std::fill_n(e_k1, array_element_count, 0.);
	std::fill_n(e2, array_element_count, 0.);
	std::fill_n(T, array_element_count, 0.);

	std::fill_n(sigma_kk, array_element_count, 0.);
	std::fill_n(u_kk, array_element_count, 0.);
	std::fill_n(v_kk, array_element_count, 0.);
	std::fill_n(e_kk, array_element_count, 0.);

	std::fill_n(sigmaX_k, array_element_count, 0.);
	std::fill_n(uX_k, array_element_count, 0.);
	std::fill_n(vY_k, array_element_count, 0.);
	std::fill_n(eR_k, array_element_count, 0.);
}
