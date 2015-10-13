#include <gtest/gtest.h>
#include <api_seq.h>
#include <api_par.h>

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
	int a = RUN_ALL_TESTS();
	//getchar();
	return a;
}

TEST(nsgas, main_test)
{
	const bool need_print = false;
	const int thread_count = 4;
	bool need_out = false;
	double abs_error = 1e-12;

	printf("Start sequential execution\n");
	double time = calculate(need_print);
	printf("Seq time = %f s.\n", time);
	printf("Finish sequential execution\n");
	printf("Start parallel execution\n");
	double time_p = calculate_parallel(need_print, thread_count);
	printf("Par time = %f s.\n", time_p);
	printf("Finish parallel execution\n");

	double* sigma_seq = get_sigma();
	double* u_seq = get_u();
	double* v_seq = get_v();
	double* e_seq = get_e();
	double* sigma_par = get_sigma_parallel();
	double* u_par = get_u_parallel();
	double* v_par = get_v_parallel();
	double* e_par = get_e_parallel();
		
	for (int i = 0; i < get_length_parallel(); i++)
	{
		ASSERT_NEAR(sigma_seq[i], sigma_par[i], abs_error);
		ASSERT_NEAR(u_seq[i], u_par[i], abs_error);
		ASSERT_NEAR(v_seq[i], v_par[i], abs_error);
		ASSERT_NEAR(e_seq[i], e_par[i], abs_error);		
	}
	
	if (need_out)
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
	printf("clear test mem\n");
	delete[] sigma_seq;
	delete[] v_seq;
	delete[] e_seq;
	delete[] sigma_par;
	delete[] u_par;
	delete[] v_par;
	delete[] e_par;
	delete[] u_seq;
	printf("run clear_memory_parallel\n");
	clear_memory_parallel(get_length_parallel());
}
