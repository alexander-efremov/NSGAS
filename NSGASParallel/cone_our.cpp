#include "math.h"
#include <stdio.h>

//M1 - количество узлов по оси х
//M - количество узлов по оси y
//(2*q-1) - количество узлов в основании клина
//(qq,cntr) - номер узла вершины клина
//tg = hx/hy
const int N = 20;
const int N1 = 10;
const int M = N + 1;
const int M1 = N1 + 1;
const int M2 = M1 * M;
const int q = 3;
const int qq = 5;
const int w = q;
const int cntr = N / 2;
const double hx = 1.0 / N1;
const double hy = 1.0 / N;
const double tau = 0.0005;
const double tg = 2;
const double Re = 10000;
const double PrRe = 0.72 * Re; // Pr * Re
const double Mah2 = 16; // Mah * Mah
const double epsilon = 0.0000000001;
// В массивах с _k1 хранятся значения функций на d-ом шаге по времени
// В массивах с _k хранятся значения функций с предыдущей итерации по нелинейности
// В массивах с _kk хранятся значения функций c (d-1) шага по времени
// Массивы с "2" использутся в итерациях метода Зейделя
// В массивах с X_k хранятся значения функций, вычисленных методом траекторий
double A[2 * M2][12];
double D[2 * M2];
double f[2 * M2];
double sigma_k[M2];
double sigma_k1[M2];
double u_k[M2];
double u_k1[M2];
double v_k[M2];
double v_k1[M2];
double B[2 * M2];
double u2[M2];
double v2[M2];
double e_k[M2];
double e_k1[M2];
double e2[M2];
double T[M2];
double sigma_kk[M2];
double u_kk[M2];
double v_kk[M2];
double e_kk[M2];
double sigmaX_k[M2];
double uX_k[M2];
double vY_k[M2];
double eR_k[M2];

double Mu(double gamma, double e_k)
{
	const double omega = 0.8;
	return pow(gamma * (gamma - 1) * Mah2 * e_k * e_k, omega);
}

double P(double gamma, double sigma_k, double e_k)
{
	return (gamma - 1) * sigma_k * sigma_k * e_k * e_k;
}

#include "trajectory_cone.h"
#include "continuity_sigma.h"
#include "motion.h"
#include "energy_epsilon.h"

void print_new_line(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure)
{
	fprintf(out, "\n\n");
	fprintf(density, "\n\n");
	fprintf(velocity, "\n\n");
	fprintf(temperature, "\n\n");
	fprintf(pressure, "\n\n");
}

void flush_file(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	fflush(out);
	fflush(density);
	fflush(velocity);
	fflush(temperature);
	fflush(pressure);
	fflush(out_itr);
}

// m = M
// m = M1
// n = N
// mah2 = Mah2
void print_to_file(const double gamma, int s_m, int s_e, int current_ts, int s_itr, int s_end, const double tau, const double hx, const double hy,
	               const int m, const int m1, const int n, const double mah2,
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

		fprintf(out, "\t\t d = %i\t d*tau = %.5f\n\n", current_ts, ts_tau);
		fprintf(out, "q = %i\t w = %i\n\n", q, w);
		fprintf(out_itr, "\t\t d = %i\t d*tau = %.5f\n\n", current_ts, ts_tau);
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

// n = N
void print_file_header(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr, const double tau, const double hx, const double hy, const int n)
{
	fprintf(out, "Cone_2D\n\n");
	fprintf(out, "N = %i\t hx = %.5f\t hy = %.5f\t tau = %.5f\n\n", n, hx, hy, tau);
	fprintf(out_itr, "Cone_2D\n\n");
	fprintf(out_itr, "N = %i\t hx = %.5f\t hy = %.5f\t tau = %.5f\n\n", n, hx, hy, tau);
	fprintf(density, "TITLE=\"density\"\n\nVARIABLES=\"x\",\"y\",\"Ro\"\n\n");
	fprintf(velocity, "TITLE=\"velocity\"\n\nVARIABLES=\"x\",\"y\",\"u\",\"v\"\n\n");
	fprintf(temperature, "TITLE=\"temperature\"\n\nVARIABLES=\"x\",\"y\",\"T\"\n\n");
	fprintf(pressure, "TITLE=\"pressure\"\n\nVARIABLES=\"x\",\"y\",\"P\"\n\n");
}

void close_files(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	fclose(out);
	fclose(density);
	fclose(velocity);
	fclose(temperature);
	fclose(pressure);
	fclose(out_itr);
}

// Initial boundary conditions with t = 0
// qq_i = q
// w_i = w
// m = M
// m1 = M1
// m2 = M2
// mah2 = Mah2
void set_initial_boundary_conditions(const double gamma, const int qq_i, const int w_i, const int m, const int m1, const int m2, const double mah2)
{
	int a;
	for (int i = 0; i < qq_i; i++)
	{
		for (int j = 0; j < m; j++)
		{
			a = i * m + j;
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

	for (int i = qq_i; i < qq_i + w_i - 1; i++)
	{
		for (int j = cntr + i - qq_i; j < m; j++)
		{
			a = i * m + j;
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

	for (int i = qq_i; i < qq_i + w_i - 1; i++)
	{
		for (int j = cntr - i + qq_i; j > -1; j--)
		{
			a = i * m + j;
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

	for (int i = qq_i + w_i - 1; i < m1; i++)
	{
		for (int j = 0; j < m; j++)
		{
			a = i * m + j;
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

	for (int i = 0; i < m2; i++)
	{
		sigmaX_k[i] = 0;
		uX_k[i] = 0;
		vY_k[i] = 0;
		eR_k[i] = 0;
	}
}

// Set array values to zero
// n = 2 * M2
// m = 12
void zeroed_arrays(int n, int m)
{
	for (int i = 0; i < 2 * n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			A[i][j] = 0;
		}
		D[i] = 0;
		B[i] = 0;
		f[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		sigma_k[i] = 0;
		sigma_k1[i] = 0;
		u_k[i] = 0;
		v_k[i] = 0;
		u_k1[i] = 0;
		v_k1[i] = 0;
		u2[i] = 0;
		v2[i] = 0;
		e_k[i] = 0;
		e_k1[i] = 0;
		e2[i] = 0;
		T[i] = 0;
		sigma_kk[i] = 0;
		u_kk[i] = 0;
		v_kk[i] = 0;
		e_kk[i] = 0;
	}
}

// qq_i = qq
// m = M
// m1 = M1
// w = w
// cntr = cntr
// n = N
int interate_over_nonlinearity(const double gamma, const int qq_i, const int m, const int m1, const int w, const int cntr, const int n, int& s_m, int& s_e, int& s_end)
{
	const int itr = 5;
	int i;
	int j;
	int a;
	int s_itr;
	for (s_itr = 1; s_itr < itr; ++s_itr)
	{
		for (i = 1; i < qq_i + 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		for (i = qq_i; i < qq_i + w - 1; i++)
		{
			for (j = cntr + i - qq_i; j < m - 1; j++)
			{
				a = i * m + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		for (i = qq_i; i < qq_i + w - 1; i++)
		{
			for (j = cntr - i + qq_i; j > 0; j--)
			{
				a = i * m + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		for (i = qq_i + w - 1; i < m1 - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		continuity(sigma_k1, u_k, v_k);
		s_m = motion(gamma, sigma_k1, sigma_k, u_k, v_k, u_k1, v_k1, u2, v2, e_k);
		s_e = energy(gamma, sigma_k1, sigma_k, u_k, v_k, u_k1, v_k1, e2, e_k, e_k1);

		if (s_m == 1 && s_e == 1)
		{
			s_end = s_itr;
			s_itr = itr;
		}

		for (j = 0; j < m; j++)
		{
			sigma_k[j] = sigma_k1[j];
			e_k[j] = e_k1[j];
			u_k[j] = u_k1[j];
			v_k[j] = v_k1[j];
		}

		for (i = 1; i < qq_i + 1; i++)
		{
			for (j = 0; j < m; j++)
			{
				a = i * m + j;
				if (j == 0)
				{
					sigma_k[a] = sigma_k1[a + 1];
					sigma_k1[a] = sigma_k1[a + 1];
					e_k[a] = e_k1[a + 1];
					e_k1[a] = e_k1[a + 1];
					u_k[a] = u_k1[a + 1];
					u_k1[a] = u_k1[a + 1];
					v_k[a] = v_k1[a + 1];
					v_k1[a] = v_k1[a + 1];
					e2[a] = e_k1[a + 1];
					u2[a] = u_k1[a + 1];
					v2[a] = v_k1[a + 1];
				}
				if (j > 0 && j < n)
				{
					sigma_k[a] = sigma_k1[a];
					e_k[a] = e_k1[a];
					u_k[a] = u_k1[a];
					v_k[a] = v_k1[a];
				}
				if (j == n)
				{
					sigma_k[a] = sigma_k1[a - 1];
					sigma_k1[a] = sigma_k1[a - 1];
					e_k[a] = e_k1[a - 1];
					e_k1[a] = e_k1[a - 1];
					u_k[a] = u_k1[a - 1];
					u_k1[a] = u_k1[a - 1];
					v_k[a] = v_k1[a - 1];
					v_k1[a] = v_k1[a - 1];
					e2[a] = e_k1[a - 1];
					u2[a] = u_k1[a - 1];
					v2[a] = v_k1[a - 1];
				}
			}
		}

		for (i = qq_i; i < qq_i + w - 1; i++)
		{
			for (j = cntr + i - qq_i; j < m; j++)
			{
				a = i * m + j;
				if (j == n)
				{
					sigma_k[a] = sigma_k1[a - 1];
					sigma_k1[a] = sigma_k1[a - 1];
					e_k[a] = e_k1[a - 1];
					e_k1[a] = e_k1[a - 1];
					u_k[a] = u_k1[a - 1];
					u_k1[a] = u_k1[a - 1];
					v_k[a] = v_k1[a - 1];
					v_k1[a] = v_k1[a - 1];
					e2[a] = e_k1[a - 1];
					u2[a] = u_k1[a - 1];
					v2[a] = v_k1[a - 1];
				}
				else
				{
					sigma_k[a] = sigma_k1[a];
					e_k[a] = e_k1[a];
					u_k[a] = u_k1[a];
					v_k[a] = v_k1[a];
				}
			}
		}

		for (i = qq_i; i < qq_i + w - 1; i++)
		{
			for (j = cntr - i + qq_i; j > -1; j--)
			{
				a = i * m + j;
				if (j == 0)
				{
					sigma_k[a] = sigma_k1[a + 1];
					sigma_k1[a] = sigma_k1[a + 1];
					e_k[a] = e_k1[a + 1];
					e_k1[a] = e_k1[a + 1];
					u_k[a] = u_k1[a + 1];
					u_k1[a] = u_k1[a + 1];
					v_k[a] = v_k1[a + 1];
					v_k1[a] = v_k1[a + 1];
					e2[a] = e_k1[a + 1];
					u2[a] = u_k1[a + 1];
					v2[a] = v_k1[a + 1];
				}
				else
				{
					sigma_k[a] = sigma_k1[a];
					e_k[a] = e_k1[a];
					u_k[a] = u_k1[a];
					v_k[a] = v_k1[a];
				}
			}
		}

		for (i = qq_i + w - 1; i < m1 - 1; i++)
		{
			for (j = 0; j < m; j++)
			{
				a = i * m + j;
				if (j == 0)
				{
					sigma_k[a] = sigma_k1[a + 1];
					sigma_k1[a] = sigma_k1[a + 1];
					e_k[a] = e_k1[a + 1];
					e_k1[a] = e_k1[a + 1];
					u_k[a] = u_k1[a + 1];
					u_k1[a] = u_k1[a + 1];
					v_k[a] = v_k1[a + 1];
					v_k1[a] = v_k1[a + 1];
					e2[a] = e_k1[a + 1];
					u2[a] = u_k1[a + 1];
					v2[a] = v_k1[a + 1];
				}
				if (j > 0 && j < n)
				{
					sigma_k[a] = sigma_k1[a];
					e_k[a] = e_k1[a];
					u_k[a] = u_k1[a];
					v_k[a] = v_k1[a];
				}
				if (j == n)
				{
					sigma_k[a] = sigma_k1[a - 1];
					sigma_k1[a] = sigma_k1[a - 1];
					e_k[a] = e_k1[a - 1];
					e_k1[a] = e_k1[a - 1];
					u_k[a] = u_k1[a - 1];
					u_k1[a] = u_k1[a - 1];
					v_k[a] = v_k1[a - 1];
					v_k1[a] = v_k1[a - 1];
					e2[a] = e_k1[a - 1];
					u2[a] = u_k1[a - 1];
					v2[a] = v_k1[a - 1];
				}
			}
		}

		for (j = 0; j < m; j++)
		{
			a = N1 * m + j;
			int indx = (N1 - 1) * m + j;

			if (j == 0)
			{
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
			if (j > 0 && j < n)
			{
				sigma_k[a] = sigma_k1[indx];
				sigma_k1[a] = sigma_k1[indx];
				e_k[a] = e_k1[indx];
				e_k1[a] = e_k1[indx];
				u_k[a] = u_k1[indx];
				u_k1[a] = u_k1[indx];
				v_k[a] = v_k1[indx];
				v_k1[a] = v_k1[indx];
				e2[a] = e_k1[indx];
				u2[a] = u_k1[indx];
				v2[a] = v_k1[indx];
			}
			if (j == n)
			{
				sigma_k[a] = sigma_k1[indx - 1];
				sigma_k1[a] = sigma_k1[indx - 1];
				e_k[a] = e_k1[indx - 1];
				e_k1[a] = e_k1[indx - 1];
				u_k[a] = u_k1[indx - 1];
				u_k1[a] = u_k1[indx - 1];
				v_k[a] = v_k1[indx - 1];
				v_k1[a] = v_k1[indx - 1];
				e2[a] = e_k1[indx - 1];
				u2[a] = u_k1[indx - 1];
				v2[a] = v_k1[indx - 1];
			}
		}
	}
	return s_itr;
}

// m = M
// m1 = M1
// qq_i = qq
// w_i = w
// cntr_i = cntr
void prepare_to_iterate(const int m, const int m1, const int qq_i, const int w_i, const int cntr_i)
{
	int i;
	int j;
	int a;
	for (i = 0; i < qq_i + 1; i++)
	{
		for (j = 0; j < m; j++)
		{
			a = i * m + j;
			sigma_k[a] = sigma_k1[a];
			e_k[a] = e_k1[a];
			u_k[a] = u_k1[a];
			v_k[a] = v_k1[a];
			sigma_kk[a] = sigma_k1[a];
			u_kk[a] = u_k1[a];
			v_kk[a] = v_k1[a];
			e_kk[a] = e_k1[a];
		}
	}
	for (i = qq_i; i < qq_i + w_i - 1; i++)
	{
		for (j = cntr_i + i - qq_i; j < m; j++)
		{
			a = i * m + j;
			sigma_k[a] = sigma_k1[a];
			e_k[a] = e_k1[a];
			u_k[a] = u_k1[a];
			v_k[a] = v_k1[a];
			sigma_kk[a] = sigma_k1[a];
			u_kk[a] = u_k1[a];
			v_kk[a] = v_k1[a];
			e_kk[a] = e_k1[a];
		}
	}
	for (i = qq_i; i < qq_i + w_i - 1; i++)
	{
		for (j = cntr_i - i + qq_i; j > -1; j--)
		{
			a = i * m + j;
			sigma_k[a] = sigma_k1[a];
			e_k[a] = e_k1[a];
			u_k[a] = u_k1[a];
			v_k[a] = v_k1[a];
			sigma_kk[a] = sigma_k1[a];
			u_kk[a] = u_k1[a];
			v_kk[a] = v_k1[a];
			e_kk[a] = e_k1[a];
		}
	}
	for (i = qq_i + w_i - 1; i < m1; i++)
	{
		for (j = 0; j < m; j++)
		{
			a = i * m + j;
			sigma_k[a] = sigma_k1[a];
			e_k[a] = e_k1[a];
			u_k[a] = u_k1[a];
			v_k[a] = v_k1[a];
			sigma_kk[a] = sigma_k1[a];
			u_kk[a] = u_k1[a];
			v_kk[a] = v_k1[a];
			e_kk[a] = e_k1[a];
		}
	}
}

int main()
{
	FILE* fout = fopen("out.txt", "w");
	FILE* fout_itr = fopen("out_itr.txt", "w");
	FILE* fdensity = fopen("density.dat", "w");
	FILE* fdensity_new = fopen("density-new.dat", "w");
	FILE* fvelocity = fopen("velocity.dat", "w");
	FILE* ftemperature = fopen("temperature.dat", "w");
	FILE* fpressure = fopen("pressure.dat", "w");
	print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr, tau, hx, hy, N);
	const double gamma = 1.4;
	const int time_steps_nbr = 1; // time_steps_nbr - количество шагов по времени
	zeroed_arrays(M2, 12);
	set_initial_boundary_conditions(gamma, qq, w, M, M1, M2, Mah2);
	for (int current_time_step = 1; current_time_step <= time_steps_nbr; current_time_step++)
	{
		int s_end = 0;
		int s_m = 0;
		int s_e = 0;
		int s_itr;
		prepare_to_iterate(M, M1, qq, w, cntr);
		s_itr = interate_over_nonlinearity(gamma, qq, M, M1, w, cntr, N, s_m, s_e, s_end);
		print_to_file(gamma, s_m, s_e, current_time_step, s_itr, s_end, tau, hx, hy, M, M1, N, Mah2, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
	}
	print_new_line(fout, fdensity, fvelocity, ftemperature, fpressure);
	close_files(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr);
	return 0;
}
