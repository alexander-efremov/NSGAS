#include "par_calculate.h"
#include "par_api.h"
#include "timer.h"

void clear_memory_parallel(const int array_element_count)
{
	for (int i = 0; i < 2 * array_element_count; i++)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] f;
	delete[] sigma_k;
	delete[] e_k;
	delete[] e_k_mu;
	delete[] e_k1;
	delete[] e_kk;
	delete[] v_kk;
	delete[] e2;
	delete[] sigma_kk;
	delete[] u_k;
	delete[] u_kk;
	delete[] v_k;
	delete[] sigma_k1;
	delete[] u_k1;
	delete[] v_k1;
	delete[] u2;
	delete[] v2;
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
		print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr, C_tau, C_hx, C_hy);
	}
	double time;

	init_arrays_parallel(C_M2, 12);	
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
	// use FastLibC from here https://software.intel.com/en-us/articles/optimizing-without-breaking-a-sweat
		memcpy(v_kk, v_k1, C_M2 * sizeof *v_kk);
		s_itr = interate_over_nonlinearity(s_m, s_e, s_end);
		if (need_print)
			print_to_file(C_gamma, s_m, s_e, current_time_step, s_itr, s_end, C_tau, C_hx, C_hy, C_M, C_M1, C_N, C_Mah2, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
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
