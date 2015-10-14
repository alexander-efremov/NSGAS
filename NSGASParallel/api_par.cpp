#include "calculate_parallel.h"
#include "api_par.h"
#include "timer.h"


void clear_memory_parallel(const int array_element_count)
{
	for (int i = 0; i < 2 * array_element_count; i++)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] B;
	delete[] D;
	delete[] f;
	delete[] sigma_k;
	delete[] e_k;
	delete[] e_k1;
	delete[] e_kk;
	delete[] v_kk;
	delete[] e2;
	delete[] T;
	delete[] sigma_kk;
	delete[] u_k;
	delete[] u_kk;
	delete[] v_k;
	delete[] sigma_k1;
	delete[] u_k1;
	delete[] v_k1;
	delete[] u2;
	delete[] v2;
	delete[] sigmaX_k;
	delete[] uX_k;
	delete[] vY_k;
	delete[] eR_k;
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
	double* r = new double[C_M2];
	for (int i = 0; i < C_M2; ++i)
	{
		r[i] = sigma_k1[i];
	}
	return r;
}

double* get_u_parallel()
{
	double* r = new double[C_M2];
	for (int i = 0; i < C_M2; ++i)
	{
		r[i] = u_k1[i];
	}
	return r;
}

double* get_v_parallel()
{
	double* r = new double[C_M2];
	for (int i = 0; i < C_M2; ++i)
	{
		r[i] = v_k1[i];
	}
	return r;
}

double* get_e_parallel()
{
	double* r = new double[C_M2];
	for (int i = 0; i < C_M2; ++i)
	{
		r[i] = e_k1[i];
	}
	return r;
}

double calculate_parallel(const bool need_print, const int thread_count)
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
		fout = fopen("out_p.txt", "C_w");
		fout_itr = fopen("out_itr_p.txt", "C_w");
		fdensity = fopen("density_p.dat", "C_w");
		fdensity_new = fopen("density-new_p.dat", "C_w");
		fvelocity = fopen("velocity_p.dat", "C_w");
		ftemperature = fopen("temperature_p.dat", "C_w");
		fpressure = fopen("pressure_p.dat", "C_w");
		print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr, C_tau, C_hx, C_hy, C_N);
	}
	const double gamma = 1.4;	
	double time;

	init_arrays(C_M2, 12);	
	set_initial_boundary_conditions(gamma, C_qq, C_w, C_M, C_M1, C_M2, C_Mah2, C_cntr);

#ifdef _OPENMP
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(omp_get_num_procs());
	printf("OPENMP THREADS COUNT = %d\n", omp_get_max_threads());
	long count = 0;
	// dummy parallel section to get all threads running
#pragma omp parallel
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
		printf("Seq TS = %d\n", current_time_step);
		prepare_to_iterate(C_M, C_M1, C_qq, C_w, C_cntr);
		s_itr = interate_over_nonlinearity(gamma, C_qq, C_M, C_M1, C_M2, C_N1, C_w, C_cntr, C_N, C_q, s_m, s_e, s_end);
		if (need_print)
			print_to_file(gamma, s_m, s_e, current_time_step, s_itr, s_end, C_tau, C_hx, C_hy, C_M, C_M1, C_N, C_Mah2, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
	}

#ifdef _OPENMP	
	time = omp_get_wtime() - time;
#else
	time = GetTimer() / 1000;
#endif

	if (need_print)
	{
		print_new_line(fout, fdensity, fvelocity, ftemperature, fpressure);
		close_files(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr);
	}
	return time;
}
