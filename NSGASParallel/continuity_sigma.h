// w_i = C_w
// qq_i = C_qq
// m = C_M
// cntr_i = C_cntr
// tau_d = C_tau
// hx_d = C_hx
// hy_d = C_hy
inline void continuity(double* sigma_k1, double* u_k, double* v_k,
                         const int qq_i, const int w_i, const int m,
                         const int cntr_i, const double tau_d,
                         const double hx_d,
                         const double hy_d)
{
	int i = 0;
	int j = 0;
	int a;
#pragma omp parallel for collapse(2) private(i, j, a)
	//Для внутренних узлов
	for (i = 1; i < qq_i; i++)
	{
		for (j = 1; j < m - 1; j++)
		{
			a = i * m + j;
			sigma_k1[a] = sigmaX_k[a] / tau_d / (1 / tau_d + (u_k[(i + 1) * m + j] - u_k[(i - 1) * m + j]) / (4 * hx_d)
				+ (v_k[a + 1] - v_k[a - 1]) / (4 * hy_d));
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = qq_i; i < qq_i + w_i; i++)
	{
		for (j = cntr_i + i + 1 - qq_i; j < m - 1; j++)
		{
			a = i * m + j;
			sigma_k1[a] = sigmaX_k[a] / tau_d / (1 / tau_d + (u_k[(i + 1) * m + j] - u_k[(i - 1) * m + j]) / (4 * hx_d)
				+ (v_k[a + 1] - v_k[a - 1]) / (4 * hy_d));
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = qq_i; i < qq_i + w_i; i++)
	{
		for (j = cntr_i - i - 1 + qq_i; j > 0; j--)
		{
			a = i * m + j;
			sigma_k1[a] = sigmaX_k[a] / tau_d / (1 / tau_d + (u_k[(i + 1) * m + j] - u_k[(i - 1) * m + j]) / (4 * hx_d)
				+ (v_k[a + 1] - v_k[a - 1]) / (4 * hy_d));
		}
	}

#pragma omp parallel for collapse(2) private(i, j, a)

	for (i = qq_i + w_i; i < C_M1 - 1; i++)
	{
		for (j = 1; j < m - 1; j++)
		{
			a = i * m + j;
			sigma_k1[a] = sigmaX_k[a] / tau_d / (1 / tau_d + (u_k[(i + 1) * m + j] - u_k[(i - 1) * m + j]) / (4 * hx_d)
				+ (v_k[a + 1] - v_k[a - 1]) / (4 * hy_d));
		}
	}

	//Для Г5. l = C_w-1; m = 1,...,C_q-2;
	i = qq_i + w_i - 1;
#pragma omp parallel for private(j, a)
	for (j = cntr_i - C_q + 2; j < cntr_i + C_q - 1; j++)
	{
		a = i * m + j;
		sigma_k1[a] = sigmaX_k[a] / (2 * tau_d) / (1 / (2 * tau_d) + (u_k[(i + 1) * m + j] - u_k[a]) / (4 * hx_d)
			+ (v_k[a + 1] - v_k[a - 1]) / (8 * hy_d));
	}

	//Для Г6.

#pragma omp parallel for private(i, j, a)
	for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
	{
		j = cntr_i + i - qq_i;
		a = i * m + j;
		sigma_k1[a] = sigmaX_k[a] * (1 / (4 * tau_d) + 1 / (4 * tau_d)) / (1 / (4 * tau_d) + 1 / (4 * tau_d) + (u_k[a] - u_k[(i - 1) * m + j]) / (8 * hx_d)
			- u_k[(i - 1) * m + j] / (16 * hx_d) + (v_k[a + 1] - v_k[a]) / (8 * hy_d) + v_k[a + 1] / (16 * hy_d));
	}

	//Для Г7.

#pragma omp parallel for private(i, j, a)
	for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
	{
		j = cntr_i - i + qq_i;
		a = i * m + j;
		sigma_k1[a] = sigmaX_k[a] * (1 / (4 * tau_d) + 1 / (4 * tau_d)) / (1 / (4 * tau_d) + 1 / (4 * tau_d) + (u_k[a] - u_k[(i - 1) * m + j]) / (8 * hx_d)
			- u_k[(i - 1) * m + j] / (16 * hx_d) + (v_k[a] - v_k[a - 1]) / (8 * hy_d) - v_k[a - 1] / (16 * hy_d));
	}

	//Для S_w-1q-1.
	i = qq_i + w_i - 1;
	j = cntr_i + i - qq_i;

	a = i * m + j;
	sigma_k1[a] = sigmaX_k[a] * (3 / (4 * tau_d) + 1 / (8 * tau_d)) / (3 / (4 * tau_d) + 1 / (8 * tau_d) + (2 * u_k[(i + 1) * m + j] - u_k[(i - 1) * m + j] - u_k[a]) / (8 * hx_d)
		+ (u_k[a] - u_k[(i - 1) * m + j]) / (16 * hx_d) + (2 * v_k[a + 1] - v_k[a - 1] - v_k[a]) / (8 * hy_d) + v_k[a] / (16 * hy_d));

	//Для S.
	i = qq_i + w_i - 1;
	j = cntr_i - i + qq_i;

	a = i * m + j;
	sigma_k1[a] = sigmaX_k[a] * (3 / (4 * tau_d) + 1 / (8 * tau_d)) / (3 / (4 * tau_d) + 1 / (8 * tau_d) + (2 * u_k[(i + 1) * m + j] - u_k[(i - 1) * m + j] - u_k[a]) / (8 * hx_d)
		+ (u_k[a] - u_k[(i - 1) * m + j]) / (16 * hx_d) + (v_k[a + 1] - 2 * v_k[a - 1] + v_k[a]) / (8 * hy_d) - v_k[a] / (16 * hy_d));

	//Для S_qq_0
	i = qq_i;
	j = cntr_i;

	a = i * m + j;
	sigma_k1[a] = sigmaX_k[a] * (1 / (4 * tau_d) + 1 / (2 * tau_d)) / (1 / (4 * tau_d) + 1 / (2 * tau_d) + (u_k[a] - u_k[(i - 1) * m + j]) / (4 * hx_d)
		+ (u_k[(i + 1) * m + j] - u_k[a]) / (8 * hx_d) + (v_k[a + 1] - v_k[a - 1]) / (8 * hy_d) + (v_k[a + 1] - v_k[a - 1]) / (16 * hy_d));
}
