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
//#define __pure __attribute__((const))
/*������ �� �� GCC �� �������� � ������ ���������� �������*/
#elif __INTEL_COMPILER
#define __pure __declspec(const)
#elif __NVCC__
#define __pure
#else
#define __pure 
#endif

//C_M1 - ���������� ����� �� ��� �
//C_M - ���������� ����� �� ��� y
//(2*C_q-1) - ���������� ����� � ��������� �����
//(C_qq,C_cntr) - ����� ���� ������� �����
//C_tg = C_hx/C_hy

// test case
//static const int C_N = 20;
//static const int C_N1 = 10;
//static const int C_q = 3;
//static const int C_qq = 5;
//static const double C_hx = 1.0 / C_N1;
//static const double C_hy = 1.0 / C_N;
//static const int time_steps_nbr = 2; // time_steps_nbr - ���������� ����� �� �������
//-----------------------
// real test 
static const int C_N = 1200;
static const int C_N1 = 800;
static const int C_q = 101;
static const int C_qq = 20;
static const double C_hx = 1.0 / 100;
static const double C_hy = 1.0 / 200;

// real test  * 2
//static const int C_N = 2400;
//static const int C_N1 = 1600;
//static const int C_q = 202;
//static const int C_qq = 40;
//static const double C_hx = 1.0 / 200;
//static const double C_hy = 1.0 / 400;

//static const int time_steps_nbr = 25000; // time_steps_nbr - ���������� ����� �� �������
static const int time_steps_nbr = 1; // time_steps_nbr - ���������� ����� �� �������
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
static const double C_gamma = 1.4;
static const double C_gamma_Mah2 = C_gamma * (C_gamma - 1) * C_Mah2;

// � �������� � _k1 �������� �������� ������� �� d-�� ���� �� �������
// � �������� � _k �������� �������� ������� � ���������� �������� �� ������������
// � �������� � _kk �������� �������� ������� c (d-1) ���� �� �������
// ������� � "2" ����������� � ��������� ������ �������
// � �������� � X_k �������� �������� �������, ����������� ������� ����������
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
static double* e_k_mu;
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

inline __pure double Mu(double value)
{
	const double omega = 0.8;
	return pow(C_gamma_Mah2 * value * value, omega);
}

inline __pure double P(double gamma, double sigma_k, double e_k)
{
	return (gamma - 1) * sigma_k * sigma_k * e_k * e_k;
}

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
					        sigma_k1[a] * sigma_k1[a], u_k1[a], v_k1[a], e_k1[a] * e_k1[a], e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * mah2), P(gamma, sigma_k1[a], e_k1[a]), Mu(e_k1[a]));
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
// C_qq = C_q
// C_w = C_w
// C_M = C_M
// C_M1 = C_M1
// m2_i = C_M2
// mah2 = C_Mah2
// C_cntr = C_cntr
inline void set_initial_boundary_conditions(const double gamma, const int C_qq, const int C_w, const int C_M, const int C_M1, const int m2_i, const double mah2, const int C_cntr)
{
	const double sqrt1 = sqrt(1 / (gamma * (gamma - 1) * mah2));

	for (int j = 0; j < C_M; j++)
	{
		for (int i = 0; i < C_qq; i++)
		{
			sigma_k1[i * C_M + j] = 1;
			e_k1[i * C_M + j] = sqrt1;
			u_k1[i * C_M + j] = i == 0 ? 1 : 0;
		}
		for (int i = C_qq + C_w - 1; i < C_M1; i++)
		{
			sigma_k1[i * C_M + j] = 1;
			e_k1[i * C_M + j] = sqrt1;
			u_k1[i * C_M + j] = 0;
		}
	}

	for (int i = C_qq; i < C_qq + C_w - 1; i++)
	{
		for (int j = C_cntr + i - C_qq; j < C_M; j++)
		{
			sigma_k1[i * C_M + j] = 1;
			e_k1[i * C_M + j] = sqrt1;
			u_k1[i * C_M + j] = 0;
		}
		for (int j = C_cntr - i + C_qq; j >= 0; j--)
		{
			sigma_k1[i * C_M + j] = 1;
			e_k1[i * C_M + j] = sqrt1;
			u_k1[i * C_M + j] = 0;
		}
	}

	for (int i = 0; i < m2_i; i++)
	{
		sigmaX_k[i] = 0;
		uX_k[i] = 0;
		vY_k[i] = 0;
		eR_k[i] = 0;
		T[i] = 1;
		v_k1[i] = 0;
		e2[i] = e_k1[i];
		u2[i] = u_k1[i];
		v2[i] = v_k1[i];
	}
}

inline __pure double trajectory(const double C_tau, const double C_hx, const double C_hy, int i, int j, double* arr, double u_k, double v_k, const int m)
{
	int idx = i * m + j;
	int idx2 = i * m + (j - 1);

	if (u_k == 0 && v_k == 0)
	{
		return 0;
	}
	if (u_k > 0 && v_k > 0)
	{
		return u_k * C_tau * ((arr[(i - 1) * m + j] - arr[idx]) / ((i - 1) * C_hx - i * C_hx))
			+ v_k * C_tau * ((arr[idx2] - arr[idx]) / ((j - 1) * C_hy - j * C_hy));
	}
	if (u_k < 0 && v_k > 0)
	{
		return u_k * C_tau * ((arr[(i + 1) * m + j] - arr[idx]) / ((i + 1) * C_hx - i * C_hx))
			+ v_k * C_tau * ((arr[idx2] - arr[idx]) / ((j - 1) * C_hy - j * C_hy));
	}
	if (u_k < 0 && v_k < 0)
	{
		return u_k * C_tau * ((arr[(i + 1) * m + j] - arr[idx]) / ((i + 1) * C_hx - i * C_hx))
			+ v_k * C_tau * ((arr[i * m + (j + 1)] - arr[idx]) / ((j + 1) * C_hy - j * C_hy));
	}
	if (u_k > 0 && v_k < 0)
	{
		return u_k * C_tau * ((arr[(i - 1) * m + j] - arr[idx]) / ((i - 1) * C_hx - i * C_hx))
			+ v_k * C_tau * ((arr[i * m + (j + 1)] - arr[idx]) / ((j + 1) * C_hy - j * C_hy));
	}
	if (u_k > 0 && v_k == 0)
	{
		return u_k * C_tau * ((arr[(i - 1) * m + j] - arr[idx]) / ((i - 1) * C_hx - i * C_hx));
	}
	if (u_k == 0 && v_k > 0)
	{
		return v_k * C_tau * ((arr[idx2] - arr[idx]) / ((j - 1) * C_hy - j * C_hy));
	}
	if (u_k < 0 && v_k == 0)
	{
		return u_k * C_tau * ((arr[(i + 1) * m + j] - arr[idx]) / ((i + 1) * C_hx - i * C_hx));
	}
	if (u_k == 0 && v_k < 0)
	{
		return v_k * C_tau * ((arr[i * m + (j + 1)] - arr[idx]) / ((j + 1) * C_hy - j * C_hy));
	}
	return 0;
}

inline void continuity(double* sigma_k1, double* u_k, double* v_k)
{
	//��� ���������� �����
#pragma omp parallel
	{
#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				sigma_k1[i * C_M + j] = sigmaX_k[i * C_M + j] / C_tau / (1 / C_tau + (u_k[(i + 1) * C_M + j] - u_k[(i - 1) * C_M + j]) / (4 * C_hx)
					+ (v_k[i * C_M + j + 1] - v_k[i * C_M + j - 1]) / (4 * C_hy));
			}
		}
#pragma omp for nowait
		for (int j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
		{
			//��� �5. l = C_w-1; m = 1,...,C_q-2;
			int i = C_qq + C_w - 1;
			sigma_k1[i * C_M + j] = sigmaX_k[i * C_M + j] / (2 * C_tau) / (1 / (2 * C_tau) + (u_k[(i + 1) * C_M + j] - u_k[i * C_M + j]) / (4 * C_hx)
				+ (v_k[i * C_M + j + 1] - v_k[i * C_M + j - 1]) / (8 * C_hy));
		}
#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//��� �6.
			int j = C_cntr + i - C_qq;
			sigma_k1[i * C_M + j] = sigmaX_k[i * C_M + j] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) / (1 / (4 * C_tau) + 1 / (4 * C_tau) + (u_k[i * C_M + j] - u_k[(i - 1) * C_M + j]) / (8 * C_hx)
				- u_k[(i - 1) * C_M + j] / (16 * C_hx) + (v_k[i * C_M + j + 1] - v_k[i * C_M + j]) / (8 * C_hy) + v_k[i * C_M + j + 1] / (16 * C_hy));
			//��� �7.
			j = C_cntr - i + C_qq;
			sigma_k1[i * C_M + j] = sigmaX_k[i * C_M + j] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) / (1 / (4 * C_tau) + 1 / (4 * C_tau) + (u_k[i * C_M + j] - u_k[(i - 1) * C_M + j]) / (8 * C_hx)
				- u_k[(i - 1) * C_M + j] / (16 * C_hx) + (v_k[i * C_M + j] - v_k[i * C_M + j - 1]) / (8 * C_hy) - v_k[i * C_M + j - 1] / (16 * C_hy));
		}
#pragma omp single
		{
			//��� S_w-1q-1.
			int i = C_qq + C_w - 1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			sigma_k1[a] = sigmaX_k[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) / (3 / (4 * C_tau) + 1 / (8 * C_tau) + (2 * u_k[(i + 1) * C_M + j] - u_k[(i - 1) * C_M + j] - u_k[a]) / (8 * C_hx)
				+ (u_k[a] - u_k[(i - 1) * C_M + j]) / (16 * C_hx) + (2 * v_k[a + 1] - v_k[a - 1] - v_k[a]) / (8 * C_hy) + v_k[a] / (16 * C_hy));
			//��� S.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			sigma_k1[a] = sigmaX_k[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) / (3 / (4 * C_tau) + 1 / (8 * C_tau) + (2 * u_k[(i + 1) * C_M + j] - u_k[(i - 1) * C_M + j] - u_k[a]) / (8 * C_hx)
				+ (u_k[a] - u_k[(i - 1) * C_M + j]) / (16 * C_hx) + (v_k[a + 1] - 2 * v_k[a - 1] + v_k[a]) / (8 * C_hy) - v_k[a] / (16 * C_hy));
			//��� S_qq_0
			i = C_qq;
			j = C_cntr;
			a = i * C_M + j;
			sigma_k1[a] = sigmaX_k[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) / (1 / (4 * C_tau) + 1 / (2 * C_tau) + (u_k[a] - u_k[(i - 1) * C_M + j]) / (4 * C_hx)
				+ (u_k[(i + 1) * C_M + j] - u_k[a]) / (8 * C_hx) + (v_k[a + 1] - v_k[a - 1]) / (8 * C_hy) + (v_k[a + 1] - v_k[a - 1]) / (16 * C_hy));
		}
	} // #pragma omp parallel
}

/*Energy*/
/*----- ������� ��������� �������� ������� ��� ��������� �������----*/
inline void nrg_calc_matrix_a(double* sigma_k1, double* e_k, const double* e_k_mu)
{
	const double c_coef1 = 2 * C_hx * C_hx * C_PrRe;
	const double c_coef2 = 2 * C_hy * C_hy * C_PrRe;
	const double c_coef3 = 4 * C_hy * C_hy * C_PrRe;
	const double c_coef4 = 8 * C_hy * C_hy * C_PrRe;
	const double c_coef5 = 2 * C_tg * C_hx * C_hy * C_PrRe;
	const double c_coef6 = 4 * C_hx * C_hx * C_PrRe;
	const double c_coef7 = 8 * C_hx * C_hx * C_PrRe;

#pragma omp parallel 
	{
#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;
					A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
					A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
					A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
						C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
						C_gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j] + e_k_mu[(i - 1) * C_M + j]) +
						C_gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
					A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
					A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));
					D[a] = 1 / A[a][2];
					f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(C_gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +
						e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +
						e_k_mu[a] / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +
						e_k_mu[a] / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
							(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
							(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
							(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
						e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
							(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
							(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
							(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
				}
			}
		}
#pragma omp for nowait 
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
			A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1])) - C_gamma / c_coef3 * (e_k_mu[a - 1] + e_k_mu[a]) - C_gamma / c_coef4 * (e_k_mu[a - 1] + 2 * e_k_mu[a]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				C_gamma / c_coef7 * (8 * e_k_mu[a] + 3 * e_k_mu[(i + 1) * C_M + j] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (8 * e_k_mu[a] + 4 * e_k_mu[a + 1] + 3 * e_k_mu[a - 1]);
			A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a])) - C_gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]) - C_gamma / c_coef7 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]);

			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
			A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				C_gamma / c_coef7 * (8 * e_k_mu[a] + 3 * e_k_mu[(i + 1) * C_M + j] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (8 * e_k_mu[a] + 3 * e_k_mu[a + 1] + 4 * e_k_mu[a - 1]);
			A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a])) - C_gamma / c_coef3 * (e_k_mu[a] + e_k_mu[a + 1]) - C_gamma / c_coef4 * (2 * e_k_mu[a] + e_k_mu[a + 1]);
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a])) - C_gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]) - C_gamma / c_coef7 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]);
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				int a = i * C_M + j;
				A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
				A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
					C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					C_gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j] + e_k_mu[(i - 1) * C_M + j]) +
					C_gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
				A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
				A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				int a = i * C_M + j;
				A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
				A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
					C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					C_gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j] + e_k_mu[(i - 1) * C_M + j]) +
					C_gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
				A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
				A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));
			}
		}
#pragma omp for nowait
		for (int j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
		{
			//��� �5. l = C_q-1; C_M = 1,...,C_q-1;	
			int i = C_qq + C_w - 1;
			int a = i * C_M + j;
			A[a][1] = C_gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - C_gamma / c_coef6 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - 2 * e_k[(i + 1) * C_M + j]) -
				C_gamma / c_coef3 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				C_gamma / c_coef6 * (2 * e_k_mu[a] + 2 * e_k_mu[(i + 1) * C_M + j]) +
				C_gamma / c_coef3 * (2 * e_k_mu[a] + e_k_mu[a - 1] + e_k_mu[a + 1]);
			A[a][3] = -C_gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[(i + 1) * C_M + j] + e_k_mu[a]));
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - P(C_gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * ((2 * u_k1[(i + 1) * C_M + j] - 2 * u_k1[a]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

				e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a])) +

				e_k_mu[a] / (12 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

				e_k_mu[a] / (8 * C_Re * e_k[a]) * (((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) * ((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) +
					((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy) * ((v_k1[(i + 1) * C_M + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy)) +

				e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		}
#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//��� �6. l = 1,...,C_q-1; C_M = C_q-1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			A[a][0] = C_gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i - 1) * C_M + j]) - C_gamma / c_coef7 * (e_k_mu[a] + 2 * e_k_mu[(i - 1) * C_M + j]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[a + 1]) +
				C_gamma / c_coef7 * (4 * e_k_mu[a] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (4 * e_k_mu[a] + 4 * e_k_mu[a + 1]);
			A[a][3] = -C_gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) - C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a])
				- C_gamma / c_coef3 * (e_k_mu[a + 1] + e_k_mu[a]) - C_gamma / c_coef4 * (2 * e_k_mu[a + 1] + e_k_mu[a]);
			A[a][5] = -C_gamma / c_coef5 * e_k_mu[a];
			A[a][8] = C_gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(C_gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a] / C_hy) -
				-P(C_gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a + 1] / C_hy) +

				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (-u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + v_k1[a] * v_k1[a]) +

				e_k_mu[a] / (16 * C_Re * e_k[a]) * ((-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy)) +

				e_k_mu[a] / (24 * C_Re * e_k[a]) * ((-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy));

			// ��� �7.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			A[a][0] = C_gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i - 1) * C_M + j]) - C_gamma / c_coef7 * (e_k_mu[a] + 2 * e_k_mu[(i - 1) * C_M + j]);
			A[a][1] = C_gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) + C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1])
				- C_gamma / c_coef3 * (e_k_mu[a] + e_k_mu[a - 1]) - C_gamma / c_coef4 * (e_k_mu[a] + 2 * e_k_mu[a - 1]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[a - 1]) +
				C_gamma / c_coef7 * (4 * e_k_mu[a] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (4 * e_k_mu[a] + 4 * e_k_mu[a - 1]);
			A[a][6] = -C_gamma / c_coef5 * e_k_mu[a];
			A[a][7] = C_gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(C_gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				-P(C_gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a - 1] / C_hy) +

				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (-u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + -v_k1[a] * -v_k1[a]) +

				e_k_mu[a] / (16 * C_Re * e_k[a]) * ((-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy)) +

				e_k_mu[a] / (24 * C_Re * e_k[a]) * ((-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy));
		}
#pragma omp single nowait
		{
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
			A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				C_gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j] + e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
			A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));

			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
			A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - C_gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				C_gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * C_M + j] + e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
			A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));

			//��� S_qq, C_N / 2 + C_q.
			j = C_cntr + i - C_qq;
			a = i * C_M + j;
			A[a][0] = C_gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]) - C_gamma / c_coef7 * (2 * e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]);
			A[a][1] = C_gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 4 * e_k[(i + 1) * C_M + j] - 3 * e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 4 * e_k[a + 1] - 2 * e_k[a - 1]) +
				C_gamma / c_coef7 * (7 * e_k_mu[a] + 4 * e_k_mu[(i + 1) * C_M + j] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (7 * e_k_mu[a] + 4 * e_k_mu[a + 1] + 2 * e_k_mu[a - 1]) + C_gamma / c_coef5 * e_k_mu[a];
			A[a][3] = -C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));
			A[a][5] = -C_gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(C_gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + 2 * v_k1[a + 1] / C_hy - v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				P(C_gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) +
				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (4 * (u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +
				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (4 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 2 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + v_k1[a] * v_k1[a]) +
				e_k_mu[a] / (16 * C_Re * e_k[a]) * (2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy)) +
				e_k_mu[a] / (24 * C_Re * e_k[a]) * (2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy));

			//��� S_qq, C_N/2 - C_q.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			A[a][0] = C_gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) + C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j])
				- C_gamma / c_coef6 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]) - C_gamma / c_coef7 * (2 * e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]);
			A[a][1] = C_gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 4 * e_k[(i + 1) * C_M + j] - 3 * e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 2 * e_k[a + 1] - 4 * e_k[a - 1]) +
				C_gamma / c_coef7 * (7 * e_k_mu[a] + 4 * e_k_mu[(i + 1) * C_M + j] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (7 * e_k_mu[a] + 2 * e_k_mu[a + 1] + 4 * e_k_mu[a - 1]) + C_gamma / c_coef5 * e_k_mu[a];
			A[a][3] = -C_gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * C_M + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * C_M + j]));
			A[a][6] = -C_gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(C_gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx + v_k1[a + 1] / C_hy + v_k1[a] / C_hy - u_k1[(i - 1) * C_M + j] / C_hx - 2 * v_k1[a - 1] / C_hy) -
				P(C_gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy) +
				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (4 * (u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +
				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (4 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + 2 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + -v_k1[a] * -v_k1[a]) +
				e_k_mu[a] / (16 * C_Re * e_k[a]) * (2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[a] / C_hy) +
					2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
				e_k_mu[a] / (24 * C_Re * e_k[a]) * (2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a] / C_hy) +
					2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));

			//��� S_qq,C_N/2
			i = C_qq;
			j = C_cntr;
			a = i * C_M + j;
			A[a][0] = C_gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * C_M + j]) - (e_k_mu[(i - 1) * C_M + j] + e_k_mu[a]));
			A[a][1] = C_gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) + C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1])
				- C_gamma / c_coef3 * (e_k_mu[a - 1] + e_k_mu[a]) - C_gamma / c_coef4 * (2 * e_k_mu[a - 1] + e_k_mu[a]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - C_gamma / c_coef7 * e_k_mu[a] / e_k[a] * (6 * e_k[a] - 4 * e_k[(i - 1) * C_M + j]) -
				C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (6 * e_k[a] - 3 * e_k[a + 1] - 3 * e_k[a - 1]) +
				C_gamma / c_coef7 * (6 * e_k_mu[a] + 4 * e_k_mu[(i - 1) * C_M + j]) +
				C_gamma / c_coef4 * (6 * e_k_mu[a] + 4 * e_k_mu[a + 1] + 4 * e_k_mu[a - 1]) - C_gamma / (C_tg * C_hx * C_hy * C_PrRe) * e_k_mu[a];
			A[a][3] = -C_gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) - C_gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a])
				- C_gamma / c_coef3 * (e_k_mu[a] + e_k_mu[a + 1]) - C_gamma / c_coef4 * (e_k_mu[a] + 2 * e_k_mu[a + 1]);
			A[a][7] = C_gamma / c_coef5 * e_k_mu[a];
			A[a][8] = C_gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - P(C_gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[a] / C_hx - 2 * u_k1[(i - 1) * C_M + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) -
				P(C_gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-2 * u_k1[a] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) +
				e_k_mu[a] / (12 * C_hx * C_hx * C_Re * e_k[a]) * (2 * (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j]) + -u_k1[a] * -u_k1[a]) +
				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +
				e_k_mu[a] / (16 * C_Re * e_k[a]) * (2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					(-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
				e_k_mu[a] / (24 * C_Re * e_k[a]) * (2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		} // #pragma omp single		
	} // #pragma omp parallel

#pragma omp parallel 
	{
		//�������� ������������ ������� ��� ������� �. ������������ � ���� ������� �� ��������� �������� ��������� ������� ��������� ������� �
#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			for (int j = C_cntr + i - C_qq; j < C_M - 1; j++)
			{
				D[i * C_M + j] = 1 / A[i * C_M + j][2];
			}
			for (int j = C_cntr - i + C_qq; j > 0; j--)
			{
				D[i * C_M + j] = 1 / A[i * C_M + j][2];
			}
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				int a = i * C_M + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(C_gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

					e_k_mu[a] / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					e_k_mu[a] / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
			for (int j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				int a = i * C_M + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(C_gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] - u_k1[a]) * (u_k1[(i + 1) * C_M + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * C_M + j]) * (u_k1[a] - u_k1[(i - 1) * C_M + j])) +

					e_k_mu[a] / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					e_k_mu[a] / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}
	}
}

//������ B = A*Xk1
inline void nrg_calc_energy(double* e_k1)
{
#pragma omp parallel 
	{
		//��� ���������� �����
#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					int a = i * C_M + j;
					B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + (j - 1)] + A[a][2] * e_k1[a] +
						A[a][3] * e_k1[i * C_M + (j + 1)] + A[a][4] * e_k1[(i + 1) * C_M + j];
					e2[a] = e_k1[a] - D[a] * (B[a] - f[a]);
				}
			}
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				int a = i * C_M + j;
				B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * C_M + (j + 1)] + A[a][4] * e_k1[(i + 1) * C_M + j];
			}
			for (int j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				int a = i * C_M + j;
				B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * C_M + (j + 1)] + A[a][4] * e_k1[(i + 1) * C_M + j];
			}
		}
#pragma omp for nowait
		for (int j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
		{
			//��� �5. l = C_q-1; C_M = 1,...,C_q-1;
			int i = C_qq + C_w - 1;
			int a = i * C_M + j;
			B[a] = A[a][1] * e_k1[a - 1] + A[a][2] * e_k1[i * C_M + j] + A[a][3] * e_k1[a + 1] +
				A[a][4] * e_k1[(i + 1) * C_M + j];
		}
#pragma omp for nowait
		for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			//��� �6. l = 1,...,C_q-1; C_M = C_q-1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][2] * e_k1[i * C_M + j] + A[a][3] * e_k1[a + 1]
				+ A[a][5] * e_k1[(i - 1) * C_M + j - 1] + A[a][8] * e_k1[(i + 1) * C_M + j + 1];
			//��� �7.
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + j - 1] + A[a][2] * e_k1[a]
				+ A[a][6] * e_k1[(i - 1) * C_M + j + 1] + A[a][7] * e_k1[(i + 1) * C_M + j - 1];
		}
#pragma omp single
		{
			//��� S_qq,C_N/2+C_q.
			int i = C_qq + C_w - 1;
			int j = C_cntr + i - C_qq;
			int a = i * C_M + j;
			B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * C_M + (j + 1)] + A[a][4] * e_k1[(i + 1) * C_M + j]
				+ A[a][5] * e_k1[(i - 1) * C_M + j - 1];

			//��� S_qq,C_N/2-C_q.
			i = C_qq + C_w - 1;
			j = C_cntr - i + C_qq;
			a = i * C_M + j;
			B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * C_M + (j + 1)] + A[a][4] * e_k1[(i + 1) * C_M + j]
				+ A[a][6] * e_k1[(i - 1) * C_M + j + 1];

			//��� S_qq,C_N/2
			i = C_qq;
			j = C_cntr;
			a = i * C_M + j;
			B[a] = A[a][0] * e_k1[(i - 1) * C_M + j] + A[a][1] * e_k1[i * C_M + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * C_M + (j + 1)]
				+ A[a][7] * e_k1[(i + 1) * C_M + j - 1] + A[a][8] * e_k1[(i + 1) * C_M + j + 1];
		}
	} // #pragma omp parallel 
#pragma omp parallel for
	//����� �����
	for (int i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		for (int j = C_cntr + i - C_qq; j < C_M - 1; j++)
		{
			e2[i * C_M + j] = e_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
		}
		for (int j = C_cntr - i + C_qq; j > 0; j--)
		{
			e2[i * C_M + j] = e_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
		}
	}
}

inline int energy(double* sigma_k1, double* e2, double* e_k, double* e_k1, const double* e_k_mu)
{
	int c;
	const int c_br = (C_N1 - 1) * (C_N - 1);
	nrg_calc_matrix_a(sigma_k1, e_k, e_k_mu);
	int s_e = 0;
	for (s_e = 0; s_e <= 20; ++s_e)
	{
		nrg_calc_energy(e_k1);
		c = 0;

#pragma omp parallel for collapse(2) reduction(+:c)
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (fabs(e_k1[i * C_M + j] - e2[i * C_M + j]) <= C_epsilon)
					++c;
			}
		}

		if (c == c_br)
			break;

#pragma omp parallel for collapse(2)
		for (int i = 1; i < C_M1 - 1; i++)
			for (int j = 1; j < C_M - 1; j++)
				e_k1[i * C_M + j] = e2[i * C_M + j];
	}

	return s_e;
}

/*End of energy*/

/* Motion */
/* ----- ������� ��������� �������� �������, ������������ ��� ���� ��������� ��������.---- */
inline void mtn_calculate_common(double* sigma_k1, double* e_k, const double* e_k_mu)
{
	const double c_coef = C_hx * C_hy * C_Re;
	const double c_coef2 = 3 * C_hx * C_hx * C_Re;
	const double c_coef3 = 2 * C_hy * C_hy * C_Re;
	const double c_coef4 = 2 * C_hx * C_hx * C_Re;
	const double c_coef5 = 3 * C_hy * C_hy * C_Re;
	const double c_coef6 = 4 * C_hy * C_hy * C_Re;
	const double c_coef7 = 6 * C_hx * C_hx * C_Re;
	const double c_coef8 = 8 * C_hy * C_hy * C_Re;

	// ��������� ��� u	
#pragma omp parallel 
	{
#pragma omp for collapse(2) nowait

		for (int i = 1; i < C_M1 - 1; ++i)
		{
			for (int j = 1; j < C_M - 1; ++j)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * C_M + j;
					A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
					A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3;
					A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
					A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
					A[a][4] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2;
					A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
					A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
					A[a][7] = (e_k_mu[i * C_M + j - 1] / 4 - e_k_mu[(i + 1) * C_M + j] / 6) / c_coef;
					A[a][8] = (e_k_mu[(i + 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j + 1)] / 4) / c_coef;
					D[a] = 1 / A[a][2];
					f[a] = uX_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - (P(C_gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(C_gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);

					// v
					a += C_M2;
					A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
					A[a][1] = -2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5;
					A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
					A[a][3] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
					A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4;
					A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
					A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
					A[a][7] = (e_k_mu[(i + 1) * C_M + j] / 4 - e_k_mu[i * C_M + j - 1] / 6) / c_coef;
					A[a][8] = (e_k_mu[i * C_M + j + 1] / 6 - e_k_mu[(i + 1) * C_M + j] / 4) / c_coef;
					D[a] = 1 / A[a][2];
				}
			}
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			// u
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef6 - (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j]) / c_coef8;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3
				+ (e_k_mu[(i + 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2 + (e_k_mu[(i + 1) * C_M + j] + 2 * e_k_mu[i * C_M + j]) / c_coef7
				+ (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j - 1)]) / c_coef6 + (2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j - 1)]) / c_coef8;
			A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
			A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 - (2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef7;
			A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
			A[a][8] = (e_k_mu[(i + 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j + 1)] / 4) / c_coef;
			A[a][9] = (1. / 4. - 1. / 8.) * e_k_mu[i * C_M + (j - 1)] / c_coef;
			A[a][11] = (1. / 12. - 1. / 6.) * e_k_mu[(i + 1) * C_M + j] / c_coef;

			// v
			a += C_M2;
			A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
			A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5 - (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j]) / (6 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5
				+ (e_k_mu[(i + 1) * C_M + j] + e_k_mu[i * C_M + j]) / (4 * C_hx * C_hx * C_Re) + (e_k_mu[(i + 1) * C_M + j] + 2 * e_k_mu[i * C_M + j]) / (8 * C_hx * C_hx * C_Re)
				+ (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j - 1)]) / c_coef5 + (2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j - 1)]) / (6 * C_hy * C_hy * C_Re);
			A[a][3] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
			A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / (4 * C_hx * C_hx * C_Re) - (2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / (8 * C_hx * C_hx * C_Re);
			A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
			A[a][8] = (e_k_mu[i * C_M + j + 1] / 6 - e_k_mu[(i + 1) * C_M + j] / 4) / c_coef;
			A[a][9] = (1. / 12. - 1. / 6.) * e_k_mu[i * C_M + (j - 1)] / c_coef;
			A[a][11] = (1. / 4. - 1. / 8.) * e_k_mu[(i + 1) * C_M + j] / c_coef;

			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3
				+ (e_k_mu[(i + 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2 + (e_k_mu[(i + 1) * C_M + j] + 2 * e_k_mu[i * C_M + j]) / c_coef7
				+ (e_k_mu[i * C_M + (j + 1)] + e_k_mu[i * C_M + j]) / c_coef6 + (e_k_mu[i * C_M + (j + 1)] + 2 * e_k_mu[i * C_M + j]) / c_coef8;
			A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef6 - (2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef8;
			A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 - (2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef7;
			A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
			A[a][7] = (e_k_mu[i * C_M + j - 1] / 4 - e_k_mu[(i + 1) * C_M + j] / 6) / c_coef;
			A[a][10] = (1. / 8. - 1. / 4.) * e_k_mu[i * C_M + (j + 1)] / c_coef;
			A[a][11] = (1. / 6. - 1. / 12.) * e_k_mu[(i + 1) * C_M + j] / c_coef;

			// v
			a += C_M2;
			A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
			A[a][1] = -2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5;
			A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5
				+ (e_k_mu[(i + 1) * C_M + j] + e_k_mu[i * C_M + j]) / (4 * C_hx * C_hx * C_Re) + (e_k_mu[(i + 1) * C_M + j] + 2 * e_k_mu[i * C_M + j]) / (8 * C_hx * C_hx * C_Re)
				+ (e_k_mu[i * C_M + (j + 1)] + e_k_mu[i * C_M + j]) / c_coef5 + (e_k_mu[i * C_M + (j + 1)] + 2 * e_k_mu[i * C_M + j]) / (6 * C_hy * C_hy * C_Re);
			A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5 - (2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / (6 * C_hy * C_hy * C_Re);
			A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / (4 * C_hx * C_hx * C_Re) - (2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / (8 * C_hx * C_hx * C_Re);
			A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
			A[a][7] = (e_k_mu[(i + 1) * C_M + j] / 4 - e_k_mu[i * C_M + j - 1] / 6) / c_coef;
			A[a][10] = (1. / 6. - 1. / 12.) * e_k_mu[i * C_M + (j + 1)] / c_coef;
			A[a][11] = (1. / 8. - 1. / 4.) * e_k_mu[(i + 1) * C_M + j] / c_coef;
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				// u
				int a = i * C_M + j;
				A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
				A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3;
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
				A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
				A[a][4] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2;
				A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
				A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
				A[a][7] = (e_k_mu[i * C_M + j - 1] / 4 - e_k_mu[(i + 1) * C_M + j] / 6) / c_coef;
				A[a][8] = (e_k_mu[(i + 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j + 1)] / 4) / c_coef;

				// v
				a += C_M2;
				A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
				A[a][1] = -2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5;
				A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
				A[a][3] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
				A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4;
				A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
				A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
				A[a][7] = (e_k_mu[(i + 1) * C_M + j] / 4 - e_k_mu[i * C_M + j - 1] / 6) / c_coef;
				A[a][8] = (e_k_mu[i * C_M + j + 1] / 6 - e_k_mu[(i + 1) * C_M + j] / 4) / c_coef;
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// u
				int a = i * C_M + j;
				A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
				A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3;
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
				A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
				A[a][4] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2;
				A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
				A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
				A[a][7] = (e_k_mu[i * C_M + j - 1] / 4 - e_k_mu[(i + 1) * C_M + j] / 6) / c_coef;
				A[a][8] = (e_k_mu[(i + 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j + 1)] / 4) / c_coef;

				// v
				a += C_M2;
				A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
				A[a][1] = -2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5;
				A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
				A[a][3] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
				A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4;
				A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
				A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
				A[a][7] = (e_k_mu[(i + 1) * C_M + j] / 4 - e_k_mu[i * C_M + j - 1] / 6) / c_coef;
				A[a][8] = (e_k_mu[i * C_M + j + 1] / 6 - e_k_mu[(i + 1) * C_M + j] / 4) / c_coef;
			}
		}
#pragma omp single nowait
		{
			// u
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
			A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
			A[a][4] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2;
			A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
			A[a][7] = (e_k_mu[i * C_M + j - 1] / 4 - e_k_mu[(i + 1) * C_M + j] / 6) / c_coef;
			A[a][8] = (e_k_mu[(i + 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j + 1)] / 4) / c_coef;
			// v
			a += C_M2;
			A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
			A[a][1] = -2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5;
			A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
			A[a][3] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
			A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4;
			A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
			A[a][7] = (e_k_mu[(i + 1) * C_M + j] / 4 - e_k_mu[i * C_M + j - 1] / 6) / c_coef;
			A[a][8] = (e_k_mu[i * C_M + j + 1] / 6 - e_k_mu[(i + 1) * C_M + j] / 4) / c_coef;

			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef3;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2 + (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
			A[a][3] = -(e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef3;
			A[a][4] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef2;
			A[a][5] = (e_k_mu[(i - 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * C_M + j + 1] / 4 - e_k_mu[(i - 1) * C_M + j] / 6) / c_coef;
			A[a][7] = (e_k_mu[i * C_M + j - 1] / 4 - e_k_mu[(i + 1) * C_M + j] / 6) / c_coef;
			A[a][8] = (e_k_mu[(i + 1) * C_M + j] / 6 - e_k_mu[i * C_M + (j + 1)] / 4) / c_coef;
			// v
			a += C_M2;
			A[a][0] = -(e_k_mu[(i - 1) * C_M + j] + e_k_mu[i * C_M + j]) / c_coef4;
			A[a][1] = -2 * (e_k_mu[i * C_M + (j - 1)] + e_k_mu[i * C_M + j]) / c_coef5;
			A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (e_k_mu[(i - 1) * C_M + j] + 2 * e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4 + 2 * (e_k_mu[i * C_M + (j - 1)] + 2 * e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
			A[a][3] = -2 * (e_k_mu[i * C_M + j] + e_k_mu[i * C_M + (j + 1)]) / c_coef5;
			A[a][4] = -(e_k_mu[i * C_M + j] + e_k_mu[(i + 1) * C_M + j]) / c_coef4;
			A[a][5] = (e_k_mu[i * C_M + j - 1] / 6 - e_k_mu[(i - 1) * C_M + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * C_M + j] / 4 - e_k_mu[i * C_M + j + 1] / 6) / c_coef;
			A[a][7] = (e_k_mu[(i + 1) * C_M + j] / 4 - e_k_mu[i * C_M + j - 1] / 6) / c_coef;
			A[a][8] = (e_k_mu[i * C_M + j + 1] / 6 - e_k_mu[(i + 1) * C_M + j] / 4) / c_coef;
		}
	} // #pragma omp parallel

	// �������� ������������ ������� ��� ������� �. ������������ � ���� ������� �� ��������� �������� ��������� ������� ��������� ������� �
#pragma omp parallel for 
	for (int i = C_qq; i < C_qq + C_w; i++)
	{
		for (int j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			D[i * C_M + j] = 1 / A[i * C_M + j][2];
			D[C_M2 + i * C_M + j] = 1 / A[C_M2 + i * C_M + j][2];
			//������ ������ ������ ������� ���������
			f[i * C_M + j] = uX_k[i * C_M + j] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(C_gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(C_gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
		}
		for (int j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			D[i * C_M + j] = 1 / A[i * C_M + j][2];
			D[C_M2 + i * C_M + j] = 1 / A[C_M2 + i * C_M + j][2];
			//������ ������ ������ ������� ���������
			f[i * C_M + j] = uX_k[i * C_M + j] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(C_gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(C_gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
		}
	} // #pragma omp parallel
}

//������ B = A*Xk1
inline void mtn_calculate_jakobi(double* u_k1, double* v_k1)
{
#pragma omp parallel
	{
#pragma omp for collapse(2) nowait
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (i < C_qq || i >= C_qq + C_w)
				{
					// u
					int a = i * C_M + j;
					B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
						A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
						A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
						A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
						A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
						A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

					u2[a] = u_k1[a] - D[a] * (B[a] - f[a]);

					// v
					a += C_M2;
					B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
						A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
						A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
						A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
						A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
						A[a][8] * u_k1[(i + 1) * C_M + j + 1];

					v2[a - C_M2] = v_k1[a - C_M2] - D[a] * (B[a] - f[a]);
				}
			}
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w - 1; i++)
		{
			// u
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)] +
				A[a][9] * v_k1[i * C_M + (j - 1)] +
				A[a][11] * v_k1[(i + 1) * C_M + j];
			// v
			a += C_M2;
			B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
				A[a][8] * u_k1[(i + 1) * C_M + j + 1] +
				A[a][9] * u_k1[i * C_M + j - 1] +
				A[a][11] * u_k1[(i + 1) * C_M + j];

			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][10] * v_k1[i * C_M + (j + 1)] +
				A[a][11] * v_k1[(i + 1) * C_M + j];
			// v
			a += C_M2;
			B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
				A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
				A[a][10] * u_k1[i * C_M + j + 1] +
				A[a][11] * u_k1[(i + 1) * C_M + j];
		}
#pragma omp for nowait
		for (int i = C_qq; i < C_qq + C_w; i++)
		{
			for (int j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				// u
				int a = i * C_M + j;
				B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
				// v
				a += C_M2;
				B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * u_k1[(i + 1) * C_M + j + 1];
			}
			for (int j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				// u
				int a = i * C_M + j;
				B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
				// v
				a += C_M2;
				B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * u_k1[(i + 1) * C_M + j + 1];
			}
		}
#pragma omp single nowait
		{
			// u
			int i = C_qq + C_w - 1;
			int j = C_cntr + i + 1 - C_qq;
			int a = i * C_M + j;
			B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

			// v
			a += C_M2;
			B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
				A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
				A[a][8] * u_k1[(i + 1) * C_M + j + 1];

			// u
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
			// v
			a += C_M2;
			B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
				A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
				A[a][8] * u_k1[(i + 1) * C_M + j + 1];
		}
	} // #pragma omp parallel
	// ����� ����� �������

#pragma omp parallel for 
	for (int i = C_qq; i < C_qq + C_w; i++)
	{
		for (int j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			u2[i * C_M + j] = u_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
			v2[i * C_M + j] = v_k1[i * C_M + j] - D[C_M2 + i * C_M + j] * (B[C_M2 + i * C_M + j] - f[C_M2 + i * C_M + j]);
		}
		for (int j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			u2[i * C_M + j] = u_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
			v2[i * C_M + j] = v_k1[i * C_M + j] - D[C_M2 + i * C_M + j] * (B[C_M2 + i * C_M + j] - f[C_M2 + i * C_M + j]);
		}
	}
}

inline int motion(double* sigma_k1,
                  double* u_k1, double* v_k1, double* u2, double* v2, double* e_k, const double* e_k_mu)
{
	int c_u;
	int c_v;
	const int break_value = (C_N1 - 1) * (C_N - 1);

	mtn_calculate_common(sigma_k1, e_k, e_k_mu);

	int s_m = 0;
	for (s_m = 0; s_m <= 20; ++s_m)
	{
		mtn_calculate_jakobi(u_k1, v_k1);

		c_u = 0;
		c_v = 0;

#pragma omp parallel for collapse(2) reduction(+:c_u, c_v)
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				if (fabs(u_k1[i * C_M + j] - u2[i * C_M + j]) <= C_epsilon)
				{
					++c_u;
				}
				if (fabs(v_k1[i * C_M + j] - v2[i * C_M + j]) <= C_epsilon)
				{
					++c_v;
				}
			}
		}

		if (c_u >= break_value && c_v >= break_value)
			break;

#pragma omp parallel for collapse(2)
		for (int i = 1; i < C_M1 - 1; i++)
		{
			for (int j = 1; j < C_M - 1; j++)
			{
				u_k1[i * C_M + j] = u2[i * C_M + j];
				v_k1[i * C_M + j] = v2[i * C_M + j];
			}
		}
	}
	return s_m;
}

/* End of motion */

inline int interate_over_nonlinearity(const double gamma,
                                      const int C_qq,
                                      const int C_M,
                                      const int C_M1,
                                      const int m2_i,
                                      const int n1_i,
                                      const int C_w,
                                      const int C_cntr,
                                      const int n_i, const int C_q, int& s_m, int& s_e, int& s_end)
{
	const int itr = 5;
	int i;
	int j;

	int s_itr;
	for (s_itr = 1; s_itr < itr; ++s_itr)
	{
#pragma omp parallel 
		{
#pragma omp for collapse(2) nowait
			for (int i = 1; i < C_M1 - 1; i++)
			{
				for (int j = 1; j < C_M - 1; j++)
				{
					sigmaX_k[i * C_M + j] = sigma_kk[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, sigma_kk, u_k[i * C_M + j], v_k[i * C_M + j], C_M);
					uX_k[i * C_M + j] = u_kk[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, u_kk, u_k[i * C_M + j], v_k[i * C_M + j], C_M);
					vY_k[i * C_M + j] = v_kk[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, v_kk, u_k[i * C_M + j], v_k[i * C_M + j], C_M);
					eR_k[i * C_M + j] = e_kk[i * C_M + j] - trajectory(C_tau, C_hx, C_hy, i, j, e_kk, u_k[i * C_M + j], v_k[i * C_M + j], C_M);
				}
			}
#pragma omp for nowait
			for (int i = 0; i < m2_i; ++i)
			{
				e_k_mu[i] = Mu(e_k[i]);
			}
		} // #pragma omp parallel

		continuity(sigma_k1, u_k, v_k);
		s_m = motion(sigma_k1, u_k1, v_k1, u2, v2, e_k, e_k_mu);
		s_e = energy(sigma_k1, e2, e_k, e_k1, e_k_mu);

		if (s_m == 1 && s_e == 1)
		{
			s_end = s_itr;
			s_itr = itr;
		}

#pragma omp parallel for private(j)
		for (j = 0; j < C_M; j++)
		{
			sigma_k[j] = sigma_k1[j];
			e_k[j] = e_k1[j];
			u_k[j] = u_k1[j];
			v_k[j] = v_k1[j];
		}

#pragma omp parallel for private(i,j)
		for (i = 1; i < C_qq + 1; i++)
		{
			for (j = 0; j < C_M; j++)
			{
				if (j == 0)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j + 1];
					sigma_k1[i * C_M + j] = sigma_k1[i * C_M + j + 1];
					e_k[i * C_M + j] = e_k1[i * C_M + j + 1];
					e_k1[i * C_M + j] = e_k1[i * C_M + j + 1];
					u_k[i * C_M + j] = u_k1[i * C_M + j + 1];
					u_k1[i * C_M + j] = u_k1[i * C_M + j + 1];
					v_k[i * C_M + j] = v_k1[i * C_M + j + 1];
					v_k1[i * C_M + j] = v_k1[i * C_M + j + 1];
					e2[i * C_M + j] = e_k1[i * C_M + j + 1];
					u2[i * C_M + j] = u_k1[i * C_M + j + 1];
					v2[i * C_M + j] = v_k1[i * C_M + j + 1];
				}
				if (j > 0 && j < n_i)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j];
					e_k[i * C_M + j] = e_k1[i * C_M + j];
					u_k[i * C_M + j] = u_k1[i * C_M + j];
					v_k[i * C_M + j] = v_k1[i * C_M + j];
				}
				if (j == n_i)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j - 1];
					sigma_k1[i * C_M + j] = sigma_k1[i * C_M + j - 1];
					e_k[i * C_M + j] = e_k1[i * C_M + j - 1];
					e_k1[i * C_M + j] = e_k1[i * C_M + j - 1];
					u_k[i * C_M + j] = u_k1[i * C_M + j - 1];
					u_k1[i * C_M + j] = u_k1[i * C_M + j - 1];
					v_k[i * C_M + j] = v_k1[i * C_M + j - 1];
					v_k1[i * C_M + j] = v_k1[i * C_M + j - 1];
					e2[i * C_M + j] = e_k1[i * C_M + j - 1];
					u2[i * C_M + j] = u_k1[i * C_M + j - 1];
					v2[i * C_M + j] = v_k1[i * C_M + j - 1];
				}
			}
		}
#pragma omp parallel for private(i,j)
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			for (j = C_cntr + i - C_qq; j < C_M; j++)
			{
				if (j == n_i)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j - 1];
					sigma_k1[i * C_M + j] = sigma_k1[i * C_M + j - 1];
					e_k[i * C_M + j] = e_k1[i * C_M + j - 1];
					e_k1[i * C_M + j] = e_k1[i * C_M + j - 1];
					u_k[i * C_M + j] = u_k1[i * C_M + j - 1];
					u_k1[i * C_M + j] = u_k1[i * C_M + j - 1];
					v_k[i * C_M + j] = v_k1[i * C_M + j - 1];
					v_k1[i * C_M + j] = v_k1[i * C_M + j - 1];
					e2[i * C_M + j] = e_k1[i * C_M + j - 1];
					u2[i * C_M + j] = u_k1[i * C_M + j - 1];
					v2[i * C_M + j] = v_k1[i * C_M + j - 1];
				}
				else
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j];
					e_k[i * C_M + j] = e_k1[i * C_M + j];
					u_k[i * C_M + j] = u_k1[i * C_M + j];
					v_k[i * C_M + j] = v_k1[i * C_M + j];
				}
			}
			for (j = C_cntr - i + C_qq; j > -1; j--)
			{
				if (j == 0)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j + 1];
					sigma_k1[i * C_M + j] = sigma_k1[i * C_M + j + 1];
					e_k[i * C_M + j] = e_k1[i * C_M + j + 1];
					e_k1[i * C_M + j] = e_k1[i * C_M + j + 1];
					u_k[i * C_M + j] = u_k1[i * C_M + j + 1];
					u_k1[i * C_M + j] = u_k1[i * C_M + j + 1];
					v_k[i * C_M + j] = v_k1[i * C_M + j + 1];
					v_k1[i * C_M + j] = v_k1[i * C_M + j + 1];
					e2[i * C_M + j] = e_k1[i * C_M + j + 1];
					u2[i * C_M + j] = u_k1[i * C_M + j + 1];
					v2[i * C_M + j] = v_k1[i * C_M + j + 1];
				}
				else
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j];
					e_k[i * C_M + j] = e_k1[i * C_M + j];
					u_k[i * C_M + j] = u_k1[i * C_M + j];
					v_k[i * C_M + j] = v_k1[i * C_M + j];
				}
			}
		}
#pragma omp parallel for private(i,j)
		for (i = C_qq + C_w - 1; i < C_M1 - 1; i++)
		{
			for (j = 0; j < C_M; j++)
			{
				if (j == 0)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j + 1];
					sigma_k1[i * C_M + j] = sigma_k1[i * C_M + j + 1];
					e_k[i * C_M + j] = e_k1[i * C_M + j + 1];
					e_k1[i * C_M + j] = e_k1[i * C_M + j + 1];
					u_k[i * C_M + j] = u_k1[i * C_M + j + 1];
					u_k1[i * C_M + j] = u_k1[i * C_M + j + 1];
					v_k[i * C_M + j] = v_k1[i * C_M + j + 1];
					v_k1[i * C_M + j] = v_k1[i * C_M + j + 1];
					e2[i * C_M + j] = e_k1[i * C_M + j + 1];
					u2[i * C_M + j] = u_k1[i * C_M + j + 1];
					v2[i * C_M + j] = v_k1[i * C_M + j + 1];
				}
				if (j > 0 && j < n_i)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j];
					e_k[i * C_M + j] = e_k1[i * C_M + j];
					u_k[i * C_M + j] = u_k1[i * C_M + j];
					v_k[i * C_M + j] = v_k1[i * C_M + j];
				}
				if (j == n_i)
				{
					sigma_k[i * C_M + j] = sigma_k1[i * C_M + j - 1];
					sigma_k1[i * C_M + j] = sigma_k1[i * C_M + j - 1];
					e_k[i * C_M + j] = e_k1[i * C_M + j - 1];
					e_k1[i * C_M + j] = e_k1[i * C_M + j - 1];
					u_k[i * C_M + j] = u_k1[i * C_M + j - 1];
					u_k1[i * C_M + j] = u_k1[i * C_M + j - 1];
					v_k[i * C_M + j] = v_k1[i * C_M + j - 1];
					v_k1[i * C_M + j] = v_k1[i * C_M + j - 1];
					e2[i * C_M + j] = e_k1[i * C_M + j - 1];
					u2[i * C_M + j] = u_k1[i * C_M + j - 1];
					v2[i * C_M + j] = v_k1[i * C_M + j - 1];
				}
			}
		}
#pragma omp parallel for private(j)
		for (j = 0; j < C_M; j++)
		{
			int indx = (n1_i - 1) * C_M + j;
			if (j == 0)
			{
				sigma_k[n1_i * C_M + j] = sigma_k1[indx + 1];
				sigma_k1[n1_i * C_M + j] = sigma_k1[indx + 1];
				e_k[n1_i * C_M + j] = e_k1[indx + 1];
				e_k1[n1_i * C_M + j] = e_k1[indx + 1];
				u_k[n1_i * C_M + j] = u_k1[indx + 1];
				u_k1[n1_i * C_M + j] = u_k1[indx + 1];
				v_k[n1_i * C_M + j] = v_k1[indx + 1];
				v_k1[n1_i * C_M + j] = v_k1[indx + 1];
				e2[n1_i * C_M + j] = e_k1[indx + 1];
				u2[n1_i * C_M + j] = u_k1[indx + 1];
				v2[n1_i * C_M + j] = v_k1[indx + 1];
			}
			if (j > 0 && j < n_i)
			{
				sigma_k[n1_i * C_M + j] = sigma_k1[indx];
				sigma_k1[n1_i * C_M + j] = sigma_k1[indx];
				e_k[n1_i * C_M + j] = e_k1[indx];
				e_k1[n1_i * C_M + j] = e_k1[indx];
				u_k[n1_i * C_M + j] = u_k1[indx];
				u_k1[n1_i * C_M + j] = u_k1[indx];
				v_k[n1_i * C_M + j] = v_k1[indx];
				v_k1[n1_i * C_M + j] = v_k1[indx];
				e2[n1_i * C_M + j] = e_k1[indx];
				u2[n1_i * C_M + j] = u_k1[indx];
				v2[n1_i * C_M + j] = v_k1[indx];
			}
			if (j == n_i)
			{
				sigma_k[n1_i * C_M + j] = sigma_k1[indx - 1];
				sigma_k1[n1_i * C_M + j] = sigma_k1[indx - 1];
				e_k[n1_i * C_M + j] = e_k1[indx - 1];
				e_k1[n1_i * C_M + j] = e_k1[indx - 1];
				u_k[n1_i * C_M + j] = u_k1[indx - 1];
				u_k1[n1_i * C_M + j] = u_k1[indx - 1];
				v_k[n1_i * C_M + j] = v_k1[indx - 1];
				v_k1[n1_i * C_M + j] = v_k1[indx - 1];
				e2[n1_i * C_M + j] = e_k1[indx - 1];
				u2[n1_i * C_M + j] = u_k1[indx - 1];
				v2[n1_i * C_M + j] = v_k1[indx - 1];
			}
		}
	}
	return s_itr;
}

// Should I try memcpy instead of for loop?
inline void prepare_to_iterate(const int m2_i)
{
#pragma omp parallel for
	for (int i = 0; i < m2_i; i++)
	{
		sigma_k[i] = sigma_k1[i];
		sigma_kk[i] = sigma_k1[i];
		e_k[i] = e_k1[i];
		u_k[i] = u_k1[i];
		v_k[i] = v_k1[i];
		u_kk[i] = u_k1[i];
		v_kk[i] = v_k1[i];
		e_kk[i] = e_k1[i];
	}
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
	e_k_mu = new double[array_element_count];
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
	std::fill_n(e_k_mu, array_element_count, 0.);
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
