#include "calculate_parallel.h"
#include "api_par.h"
#include "timer.h"

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

double calculate_parallel(bool need_print)
{
	FILE* fout = nullptr;
	FILE* fout_itr = nullptr;
	FILE* fdensity = nullptr;
	FILE* fdensity_new = nullptr;
	FILE* fvelocity = nullptr;
	FILE* ftemperature = nullptr;
	FILE* fpressure = nullptr;
	if (need_print)
	{
		fout = fopen("out_p.txt", "w");
		fout_itr = fopen("out_itr_p.txt", "w");
		fdensity = fopen("density_p.dat", "w");
		fdensity_new = fopen("density-new_p.dat", "w");
		fvelocity = fopen("velocity_p.dat", "w");
		ftemperature = fopen("temperature_p.dat", "w");
		fpressure = fopen("pressure_p.dat", "w");
		print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr, tau, hx, hy, N);
	}
	const double gamma = 1.4;
	const int time_steps_nbr = 1; // time_steps_nbr - количество шагов по времени
	zeroed_arrays(M2, 12);
	set_initial_boundary_conditions(gamma, qq, w, M, M1, M2, Mah2);
	StartTimer();
	for (int current_time_step = 1; current_time_step <= time_steps_nbr; current_time_step++)
	{
		int s_end = 0;
		int s_m = 0;
		int s_e = 0;
		int s_itr;
		prepare_to_iterate(M, M1, qq, w, cntr);
		s_itr = interate_over_nonlinearity(gamma, qq, M, M1, w, cntr, N, q, s_m, s_e, s_end);
		if (need_print)
			print_to_file(gamma, s_m, s_e, current_time_step, s_itr, s_end, tau, hx, hy, M, M1, N, Mah2, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
	}
	double time = GetTimer();
	if (need_print)
	{
		print_new_line(fout, fdensity, fvelocity, ftemperature, fpressure);
		close_files(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr);
	}
	return time;
}
