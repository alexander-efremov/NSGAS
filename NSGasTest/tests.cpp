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
	const int thread_count = 8;
	bool need_out = get_length_parallel() < 300;
	double abs_error12 = 1e-12;
	double abs_error11 = 1e-11;

	printf("Start sequential execution\n");
	double time = calculate(need_print);
printf("Seq time = %f s.\n", time);
	printf("Finish sequential execution\n");
	printf("Start parallel execution\n");
	double time_p = calculate_parallel(need_print, thread_count);
	printf("Par time = %f s.\n", time_p);
	printf("Finish parallel execution\n");
	printf("Seq time / par time = %f\n", time / time_p);

	double* sigma_seq = get_sigma();
	double* u_seq = get_u();
	double* v_seq = get_v();
	double* e_seq = get_e();
	double* sigma_par = get_sigma_parallel();
	double* u_par = get_u_parallel();
	double* v_par = get_v_parallel();
	double* e_par = get_e_parallel();
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

	/*
	Comparing floating-point numbers is tricky. Due to round-off errors, it is very unlikely that two floating-points will match exactly.
	Therefore, ASSERT_EQ 's naive comparison usually doesn't work.
	And since floating-points can have a wide value range, no single fixed error bound works.
	It's better to compare by a fixed relative error bound, except for values close to 0 due to the loss of precision there.
	*/
	// пока что v не считается... И что - то GTF не умеет работать с нулями double нормально

	for (int i = 0; i < get_length_parallel(); i++)
		ASSERT_NEAR(v_seq[i], v_par[i], abs_error12);

	for (int i = 0; i < get_length_parallel(); i++)
	{
		ASSERT_NEAR(sigma_seq[i], sigma_par[i], abs_error12);
		ASSERT_NEAR(u_seq[i], u_par[i], abs_error12);
		ASSERT_NEAR(e_seq[i], e_par[i], abs_error12);
	}

	printf("run clear_memory_parallel\n");
	clear_memory_parallel(get_length_parallel());
	clear_memory(get_length());
}

