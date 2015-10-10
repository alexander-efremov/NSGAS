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

void print_to_file(const double gamma, int s_m, int s_e, int d, int s_itr, int s_end, FILE* out, FILE* density, FILE* density_new, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	int i;
	int j;
	int a;
	if (d / 1. == 1 || /*(d/2. == 1)||(d/3. == 1)||(d/4. == 1)||(d/5. == 1)||(d/10. == 1)||*/d / 100. == 1/*||(d/200. == 1)||(d/300. == 1)||(d/400. == 1)*/ || d / 500. == 1
		/*||(d/600. == 1)||(d/800. == 1)*/ || d / 1000. == 1/*||(d/1200. == 1)||(d/1500. == 1)||(d/1700. == 1)*/
		|| d / 2000. == 1/*||(d/2200. == 1)||(d/2500. == 1)||(d/2600. == 1)||(d/2700. == 1)||(d/2800. == 1)||(d/2900. == 1)*/
		|| d / 3000. == 1/*|| (d/3100. == 1)||(d/3200. == 1)|| (d/3300. == 1)|| (d/3400. == 1)*/
		|| d / 3500. == 1/*||(d/3600. == 1)||(d/3700. == 1)||(d/3800. == 1)||(d/3900. == 1)*/ || d / 4000. == 1/*||(d/4100. == 1)||(d/4200. == 1)||(d/4300. == 1)||(d/4400. == 1)*/ || d / 4500. == 1/*||(d/4600. == 1)||(d/4700. == 1)||(d/4800. == 1)||(d/4900. == 1)*/ || d / 5000. == 1/*|| (d/5100. == 1) || (d/5200. == 1)
																																																																						   || (d/5300. == 1)|| (d/5400. == 1)*/ || d / 5500. == 1/*|| (d/5600. == 1)|| (d/5700. == 1)|| (d/5800. == 1)|| (d/5900. == 1)*/ || d / 6000. == 1/*|| (d/6100. == 1)|| (d/6200. == 1)|| (d/6300. == 1)|| (d/6400. == 1)*/ || d / 6500. == 1/*|| (d/6700. == 1)*/ || d / 7000. == 1 || d / 7500. == 1 ||
		d / 8000. == 1 || d / 8500. == 1 || d / 9000. == 1 || d / 9500. == 1
		|| d / 10000. == 1 || d / 11000. == 1 || d / 12000. == 1
		|| d / 13000. == 1 || d / 14000. == 1 || d / 15000. == 1 || d / 16000. == 1 || d / 17000. == 1 || d / 18000. == 1 || d / 19000. == 1
		|| d / 20000. == 1/*|| (d/21000. == 1)|| (d/22000. == 1)*/ || d / 23000. == 1/*|| (d/24000. == 1)*/ || d / 25000. == 1
		/*|| (d/26000. == 1)*/ || d / 28000. == 1 || d / 30000. == 1 || d / 33000. == 1 || d / 35000. == 1 || d / 38000. == 1
		|| d / 40000. == 1 || d / 43000. == 1 || d / 45000. == 1 || d / 48000. == 1 || d / 50000. == 1 || d / 53000. == 1 || d / 55000. == 1
		|| d / 58000. == 1 || d / 60000. == 1 || d / 63000. == 1 || d / 65000. == 1 || d / 68000. == 1 || d / 70000. == 1 || d / 73000. == 1
		|| d / 75000. == 1 || d / 78000. == 1 || d / 80000. == 1 || d / 83000. == 1 || d / 85000. == 1 || d / 88000. == 1 || d / 90000. == 1
		|| d / 93000. == 1 || d / 95000. == 1 || d / 98000. == 1 || d / 100000. == 1
		|| d / 103000. == 1 || d / 105000. == 1 || d / 108000. == 1 || d / 110000. == 1 || d / 113000. == 1 || d / 115000. == 1 || d / 118000. == 1
		|| d / 120000. == 1 || d / 123000. == 1 || d / 125000. == 1 || d / 128000. == 1 || d / 130000. == 1 || d / 133000. == 1 || d / 135000. == 1
		|| d / 138000. == 1 || d / 140000. == 1 || d / 143000. == 1 || d / 145000. == 1 || d / 148000. == 1 || d / 150000. == 1
		|| d / 153000. == 1 || d / 155000. == 1 || d / 158000. == 1 || d / 160000. == 1 || d / 163000. == 1 || d / 165000. == 1 || d / 168000. == 1
		|| d / 170000. == 1 || d / 173000. == 1 || d / 175000. == 1 || d / 178000. == 1 || d / 180000. == 1 || d / 183000. == 1 || d / 185000. == 1
		|| d / 188000. == 1 || d / 190000. == 1 || d / 193000. == 1 || d / 195000. == 1 || d / 198000. == 1 || d / 200000. == 1)
	{
		fprintf(out, "\t\t d = %i\t d*tau = %.5f\n\n", d, d * tau);
		fprintf(out, "q = %i\t w = %i\n\n", q, w);
		fprintf(out_itr, "\t\t d = %i\t d*tau = %.5f\n\n", d, d * tau);
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
		fprintf(density, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", N, d * tau, M1, M);
		fprintf(velocity, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", N, d * tau, M1, M);
		fprintf(temperature, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", N, d * tau, M1, M);
		fprintf(pressure, "ZONE T=\"n = %i t = %.4f\",I=%i,J=%i,ZONETYPE=ORDERED,DATAPACKING=POINT\n\n", N, d * tau, M1, M);

		for (j = 0; j < M; j++)
		{
			for (i = 0; i < M1; i++)
			{
				a = i * M + j;
				double ihx = i * hx;
				double jhy = j * hy;
				fprintf(density, "%.3f\t %.4f\t %.10f\n", ihx, jhy, sigma_k1[a] * sigma_k1[a]);
				fprintf(density_new, "%.3f\t %.4f\t %.15e\n", ihx, jhy, sigma_k1[a] * sigma_k1[a]);
				fprintf(velocity, "%.3f\t %.4f\t%.10f\t %.10f\n", ihx, jhy, u_k1[a], v_k1[a]);
				fprintf(temperature, "%.3f\t %.4f\t %.10f\n", ihx, jhy, e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * Mah2));
				fprintf(pressure, "%.3f\t %.4f\t %.10f\n", ihx, jhy, sigma_k1[a] * sigma_k1[a] * e_k1[a] * e_k1[a] * (gamma - 1));
			}
		}

		if (d == 1)
		{
			for (i = 0; i < M1; i++)
			{
				for (j = 0; j < M; j++)
				{
					a = i * M + j;
					fprintf(out, "i=%i j=%i\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\n", i, j, sigma_k1[a] * sigma_k1[a], u_k1[a], v_k1[a], e_k1[a] * e_k1[a], e_k1[a] * e_k1[a] * (gamma * (gamma - 1) * Mah2), P(gamma, sigma_k1[a], e_k1[a]), Mu(gamma, e_k1[a]));
				}
			}
		}
		print_new_line(out, density, velocity, temperature, pressure);
		flush_file(out, density, velocity, temperature, pressure, out_itr);
	}
}

void print_file_header(FILE* out, FILE* density, FILE* velocity, FILE* temperature, FILE* pressure, FILE* out_itr)
{
	fprintf(out, "Cone_2D\n\n");
	fprintf(out, "N = %i\t hx = %.5f\t hy = %.5f\t tau = %.5f\n\n", N, hx, hy, tau);
	fprintf(out_itr, "Cone_2D\n\n");
	fprintf(out_itr, "N = %i\t hx = %.5f\t hy = %.5f\t tau = %.5f\n\n", N, hx, hy, tau);
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
void set_initial_boundary_conditions(const double gamma, const int qq_i, const int w_i, const int m, const int m1, const int m2)
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
				e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * Mah2));
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
				e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * Mah2));
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
			e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * Mah2));
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
			e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * Mah2));
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
			e_k1[a] = sqrt(T[a] / (gamma * (gamma - 1) * Mah2));
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
	for (int i = 0; i < 2*n; i++)
	{
		for ( int j = 0; j < m; j++)
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

int interate_over_nonlinearity(const double gamma, int& s_m, int& s_e, int& s_end)
{
	const int itr = 5;
	int i;
	int j;
	int a;
	int s_itr;
	for (s_itr = 1; s_itr < itr; ++s_itr)
	{
		for (i = 1; i < qq + 1; i++)
		{
			for (j = 1; j < M - 1; j++)
			{
				a = i * M + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		for (i = qq; i < qq + w - 1; i++)
		{
			for (j = cntr + i - qq; j < M - 1; j++)
			{
				a = i * M + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		for (i = qq; i < qq + w - 1; i++)
		{
			for (j = cntr - i + qq; j > 0; j--)
			{
				a = i * M + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		for (i = qq + w - 1; i < M1 - 1; i++)
		{
			for (j = 1; j < M - 1; j++)
			{
				a = i * M + j;
				sigmaX_k[a] = sigma_kk[a] - trajectory(i, j, sigma_kk, u_k[a], v_k[a]);
				uX_k[a] = u_kk[a] - trajectory(i, j, u_kk, u_k[a], v_k[a]);
				vY_k[a] = v_kk[a] - trajectory(i, j, v_kk, u_k[a], v_k[a]);
				eR_k[a] = e_kk[a] - trajectory(i, j, e_kk, u_k[a], v_k[a]);
			}
		}

		continuity(sigma_k, sigma_k1, u_k, v_k);
		s_m = motion(gamma, sigma_k1, sigma_k, u_k, v_k, u_k1, v_k1, u2, v2, e_k);
		s_e = energy(gamma, sigma_k1, sigma_k, u_k, v_k, u_k1, v_k1, e2, e_k, e_k1);

		if (s_m == 1 && s_e == 1)
		{
			s_end = s_itr;
			s_itr = itr;
		}

		for (j = 0; j < M; j++)
		{
			sigma_k[j] = sigma_k1[j];
			e_k[j] = e_k1[j];
			u_k[j] = u_k1[j];
			v_k[j] = v_k1[j];
		}

		for (i = 1; i < qq + 1; i++)
		{
			for (j = 0; j < M; j++)
			{
				a = i * M + j;
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
				if (j > 0 && j < N)
				{
					sigma_k[a] = sigma_k1[a];
					e_k[a] = e_k1[a];
					u_k[a] = u_k1[a];
					v_k[a] = v_k1[a];
				}
				if (j == N)
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

		for (i = qq; i < qq + w - 1; i++)
		{
			for (j = cntr + i - qq; j < M; j++)
			{
				a = i * M + j;
				if (j == N)
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

		for (i = qq; i < qq + w - 1; i++)
		{
			for (j = cntr - i + qq; j > -1; j--)
			{
				a = i * M + j;
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

		for (i = qq + w - 1; i < M1 - 1; i++)
		{
			for (j = 0; j < M; j++)
			{
				a = i * M + j;
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
				if (j > 0 && j < N)
				{
					sigma_k[a] = sigma_k1[a];
					e_k[a] = e_k1[a];
					u_k[a] = u_k1[a];
					v_k[a] = v_k1[a];
				}
				if (j == N)
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

		for (j = 0; j < M; j++)
		{
			a = N1 * M + j;
			int indx = (N1 - 1) * M + j;

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
			if (j > 0 && j < N)
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
			if (j == N)
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

void prepare_to_iterate()
{
	int i;
	int j;
	int a;
	for (i = 0; i < qq + 1; i++)
	{
		for (j = 0; j < M; j++)
		{
			a = i * M + j;
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
	for (i = qq; i < qq + w - 1; i++)
	{
		for (j = cntr + i - qq; j < M; j++)
		{
			a = i * M + j;
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
	for (i = qq; i < qq + w - 1; i++)
	{
		for (j = cntr - i + qq; j > -1; j--)
		{
			a = i * M + j;
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
	for (i = qq + w - 1; i < M1; i++)
	{
		for (j = 0; j < M; j++)
		{
			a = i * M + j;
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
	print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr);
	const double gamma = 1.4;
	const int time_steps_nbr = 1; // time_steps_nbr - количество шагов по времени
	zeroed_arrays(M2, 12);
	set_initial_boundary_conditions(gamma, qq, w, M, M1, M2);
	for (int current_time_step = 1; current_time_step <= time_steps_nbr; current_time_step++)
	{
		int s_end = 0;
		int s_m = 0;
		int s_e = 0;
		int s_itr;
		prepare_to_iterate();		
		s_itr = interate_over_nonlinearity(gamma, s_m, s_e, s_end);
		print_to_file(gamma, s_m, s_e, current_time_step, s_itr, s_end, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
	}	
	print_new_line(fout, fdensity, fvelocity, ftemperature, fpressure);
	close_files(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr);
	return 0;
}
