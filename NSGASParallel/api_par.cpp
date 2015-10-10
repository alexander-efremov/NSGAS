#include "calculate_parallel.h"
#include "api_par.h"

int get_length_parallel()
{
	return M2;
}

int get_length_parallel_x()
{
	return M1;
}

int get_length_parallel_y()
{
	return M;
}

double* get_sigma_parallel()
{
	double* r = new double[M2];
	for (int i = 0; i < M2; ++i)
	{
		r[i] = sigma_k1[i];
	}
	return r;
}

double* get_u_parallel()
{
	double* r = new double[M2];
	for (int i = 0; i < M2; ++i)
	{
		r[i] = u_k1[i];
	}
	return r;
}

double* get_v_parallel()
{
	double* r = new double[M2];
	for (int i = 0; i < M2; ++i)
	{
		r[i] = v_k1[i];
	}
	return r;
}

double* get_e_parallel()
{
	double* r = new double[M2];
	for (int i = 0; i < M2; ++i)
	{
		r[i] = e_k1[i];
	}
	return r;
}

void calculate_parallel()
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
}