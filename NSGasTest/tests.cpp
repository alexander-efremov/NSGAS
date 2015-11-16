#include <gtest/gtest.h>
#include <seq_api.h>

extern "C"
{
#include <par_api.h>
}


using namespace std;
using namespace ::testing;

inline void _print_matrix(double* a, int n, int m, int precision = 8)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			int k = i * n + j;
			switch (precision)
			{
			case 1:
				printf("%.1f ", a[k]);
				break;
			case 2:
				printf("%.2f ", a[k]);
				break;
			case 3:
				printf("%.3f ", a[k]);
				break;
			case 4:
				printf("%.4f ", a[k]);
				break;
			case 5:
				printf("%.5f ", a[k]);
				break;
			case 6:
				printf("%.6f ", a[k]);
				break;
			case 7:
				printf("%.7f ", a[k]);
				break;
			case 8:
				printf("%.8f ", a[k]);
				break;
			}
		}
		printf("\n");
	}
}

int main(int ac, char* av [])
{
	InitGoogleTest(&ac, av);
	return RUN_ALL_TESTS();
}

TEST(nsgas, main_test)
{
	const bool need_print = false;
	const bool need_seq = true;
	
	const int thread_count = 4;
	bool need_out = get_length_parallel() < 300;
	double abs_error = 1e-12;
	double time;
	if (need_seq)
	{
		printf("Start sequential execution\n");
		time = calculate(need_print);
		printf("Seq time = %f s.\n", time);
		printf("Finish sequential execution\n");
	}	
	printf("Start parallel execution\n");
	double time_p = calculate_parallel(need_print, thread_count);
	printf("Par time = %f s.\n", time_p);
	printf("Finish parallel execution\n");
	if (need_seq)
		printf("Seq time / par time = %f\n", time / time_p);
	double* sigma_seq;
	double* u_seq;
	double* v_seq;
	double* e_seq;

	if (need_seq)
	{
		sigma_seq = get_sigma();
		u_seq = get_u();
		v_seq = get_v();
		e_seq = get_e();
	}
	double* sigma_par = get_sigma_parallel();
	double* u_par = get_u_parallel();
	double* v_par = get_v_parallel();
	double* e_par = get_e_parallel();
	if (need_out && need_seq)
	{
		printf("Sigma Seq\n");
		_print_matrix(sigma_seq, get_length_x(), get_length_y());
		printf("Sigma Par\n");
		_print_matrix(sigma_par, get_length_parallel_x(), get_length_parallel_y());
		printf("U Seq\n");
		_print_matrix(u_seq, get_length_x(), get_length_y());
		printf("U Par\n");
		_print_matrix(u_par, get_length_parallel_x(), get_length_parallel_y());
		printf("V Seq\n");
		_print_matrix(v_seq, get_length_x(), get_length_y());
		printf("V Par\n");
		_print_matrix(v_par, get_length_parallel_x(), get_length_parallel_y());
		printf("E Seq\n");
		_print_matrix(e_seq, get_length_x(), get_length_y());
		printf("E Par\n");
		_print_matrix(e_par, get_length_parallel_x(), get_length_parallel_y());
	}

	
	if (need_seq)
	for (int i = 0; i < get_length_parallel(); i++)
	{
		ASSERT_NEAR(v_seq[i], v_par[i], abs_error);
		ASSERT_NEAR(sigma_seq[i], sigma_par[i], abs_error);
		ASSERT_NEAR(u_seq[i], u_par[i], abs_error);
		ASSERT_NEAR(e_seq[i], e_par[i], abs_error);
	}

	clear_memory_parallel(get_length_parallel());
	if (need_seq)
		clear_memory(get_length());
}

