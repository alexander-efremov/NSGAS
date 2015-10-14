/*----- Функция заполняет элементы матрицы для уравнения энергии----*/
// This method requires refactoring: change i * m + j to a
// m = C_M
// qq = C_qq
// w_i = C_w
// cntr_i = C_cntr
inline void energy_a(const double gamma, const int m, const int qq_i, const int w_i, const int cntr_i, const int q_i, double* sigma_k1, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;
	const double c_coef1 = 2 * C_hx * C_hx * C_PrRe;
	const double c_coef2 = 2 * C_hy * C_hy * C_PrRe;	
	const double c_coef3 = 4 * C_hy * C_hy * C_PrRe;
	const double c_coef4 = 8 * C_hy * C_hy * C_PrRe;
	const double c_coef5 = 2 * C_tg * C_hx * C_hy * C_PrRe;
	const double c_coef6 = 4 * C_hx * C_hx * C_PrRe;

// begin inner nodes
#pragma omp parallel					
	{
#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = 1; i < qq_i; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
				A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
					gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					gamma / c_coef1 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[(i - 1) * m + j])) +
					gamma / c_coef2 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a - 1]));
				A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
				A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + w_i - 1; i++)
		{
			j = cntr_i + i + 1 - qq_i;
			a = i * m + j;
			A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
			A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1])) - gamma / c_coef3 * (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])) - gamma / c_coef4 * (Mu(gamma, e_k[a - 1]) + 2 * Mu(gamma, e_k[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
				gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / (8 * C_hx * C_hx * C_PrRe) * (8 * Mu(gamma, e_k[a]) + 3 * Mu(gamma, e_k[(i + 1) * m + j]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
				gamma / c_coef4 * (8 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[a + 1]) + 3 * Mu(gamma, e_k[a - 1]));
			A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
			A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a])) - gamma / c_coef6 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]));
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + w_i - 1; i++)
		{
			j = cntr_i - i - 1 + qq_i;
			a = i * m + j;
			A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
			A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
				gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / (8 * C_hx * C_hx * C_PrRe) * (8 * Mu(gamma, e_k[a]) + 3 * Mu(gamma, e_k[(i + 1) * m + j]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
				gamma / c_coef4 * (8 * Mu(gamma, e_k[a]) + 3 * Mu(gamma, e_k[a + 1]) + 4 * Mu(gamma, e_k[a - 1]));
			A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a])) - gamma / c_coef3 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])) - gamma / c_coef4 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]));
			A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a])) - gamma / c_coef6 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]));
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i + i + 2 - qq_i; j < m - 1; j++)
			{
				a = i * m + j;
				A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
				A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
					gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					gamma / c_coef1 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[(i - 1) * m + j])) +
					gamma / c_coef2 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a - 1]));
				A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
				A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i - i - 2 + qq_i; j > 0; j--)
			{
				a = i * m + j;
				A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
				A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
					gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					gamma / c_coef1 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[(i - 1) * m + j])) +
					gamma / c_coef2 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a - 1]));
				A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
				A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = qq_i + w_i; i < C_M1 - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
				A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
					gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					gamma / c_coef1 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[(i - 1) * m + j])) +
					gamma / c_coef2 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a - 1]));
				A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
				A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));
			}
		}
	} // #pragma omp parallel

	i = qq_i + w_i - 1;
	j = cntr_i + i + 1 - qq_i;
	a = i * m + j;
	A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
	A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
		gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
		gamma / c_coef1 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[(i - 1) * m + j])) +
		gamma / c_coef2 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a - 1]));
	A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
	A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));

	i = qq_i + w_i - 1;
	j = cntr_i - i - 1 + qq_i;
	a = i * m + j;
	A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
	A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m + j] - e_k[(i - 1) * m + j]) -
		gamma / c_coef2 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
		gamma / c_coef1 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[(i - 1) * m + j])) +
		gamma / c_coef2 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a - 1]));
	A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
	A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));

	// end inner nodes


	//Для Г5. l = C_q-1; m = 1,...,C_q-1;	
	i = qq_i + w_i - 1;	
#pragma omp parallel 	
	{
#pragma omp for private(j, a) nowait
		for (j = cntr_i - q_i + 2; j < cntr_i + q_i - 1; j++)
		{
			a = i * m + j;
			A[a][1] = gamma / c_coef3 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - gamma / c_coef6 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - 2 * e_k[(i + 1) * m + j]) -
				gamma / c_coef3 * Mu(gamma, e_k[a]) / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / c_coef6 * (2 * Mu(gamma, e_k[a]) + 2 * Mu(gamma, e_k[(i + 1) * m + j])) +
				gamma / c_coef3 * (2 * Mu(gamma, e_k[a]) + Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a + 1]));
			A[a][3] = -gamma / c_coef3 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
			A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[(i + 1) * m + j]) + Mu(gamma, e_k[a])));
		}

		//Для Г6. l = 1,...,C_q-1; m = C_q-1;
#pragma omp for private(i, j, a) nowait
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			j = cntr_i + i - qq_i;
			a = i * m + j;
			A[a][0] = gamma / c_coef6 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j])
				- gamma / c_coef6 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i - 1) * m + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[a]) + 2 * Mu(gamma, e_k[(i - 1) * m + j]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (4 * e_k[a] - 3 * e_k[(i - 1) * m + j]) -
				gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (4 * e_k[a] - 3 * e_k[a + 1]) +
				gamma / (8 * C_hx * C_hx * C_PrRe) * (4 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
				gamma / c_coef4 * (4 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[a + 1]));
			A[a][3] = -gamma / c_coef3 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) - gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a])
				- gamma / c_coef3 * (Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a])) - gamma / c_coef4 * (2 * Mu(gamma, e_k[a + 1]) + Mu(gamma, e_k[a]));
			A[a][5] = -gamma / c_coef5 * Mu(gamma, e_k[a]);
			A[a][8] = gamma / c_coef5 * Mu(gamma, e_k[a]);
		}

		//Для Г7.
#pragma omp for private(i, j, a) nowait
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			j = cntr_i - i + qq_i;
			a = i * m + j;
			A[a][0] = gamma / c_coef6 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j])
				- gamma / c_coef6 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i - 1) * m + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[a]) + 2 * Mu(gamma, e_k[(i - 1) * m + j]));
			A[a][1] = gamma / c_coef3 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) + gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1])
				- gamma / c_coef3 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a - 1])) - gamma / c_coef4 * (Mu(gamma, e_k[a]) + 2 * Mu(gamma, e_k[a - 1]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (4 * e_k[a] - 3 * e_k[(i - 1) * m + j]) -
				gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (4 * e_k[a] - 3 * e_k[a - 1]) +
				gamma / (8 * C_hx * C_hx * C_PrRe) * (4 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
				gamma / c_coef4 * (4 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[a - 1]));
			A[a][6] = -gamma / c_coef5 * Mu(gamma, e_k[a]);
			A[a][7] = gamma / c_coef5 * Mu(gamma, e_k[a]);
		}
	} // #pragma omp parallel

	//Для S_qq,C_N/2+C_q.
	i = qq_i + w_i - 1;
	j = cntr_i + i - qq_i;
	a = i * m + j;
	A[a][0] = gamma / c_coef6 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j])
		- gamma / c_coef6 * (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a]));
	A[a][1] = gamma / c_coef3 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (7 * e_k[a] - 4 * e_k[(i + 1) * m + j] - 3 * e_k[(i - 1) * m + j]) -
		gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (7 * e_k[a] - 4 * e_k[a + 1] - 2 * e_k[a - 1]) +
		gamma / (8 * C_hx * C_hx * C_PrRe) * (7 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[(i + 1) * m + j]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
		gamma / c_coef4 * (7 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[a + 1]) + 2 * Mu(gamma, e_k[a - 1])) + gamma / c_coef5 * Mu(gamma, e_k[a]);
	A[a][3] = -gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
	A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));
	A[a][5] = -gamma / c_coef5 * Mu(gamma, e_k[a]);

	//Для S_qq,C_N/2-C_q.
	i = qq_i + w_i - 1;
	j = cntr_i - i + qq_i;
	a = i * m + j;
	A[a][0] = gamma / c_coef6 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j])
		- gamma / c_coef6 * (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a]));
	A[a][1] = gamma / c_coef2 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) - (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (7 * e_k[a] - 4 * e_k[(i + 1) * m + j] - 3 * e_k[(i - 1) * m + j]) -
		gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (7 * e_k[a] - 2 * e_k[a + 1] - 4 * e_k[a - 1]) +
		gamma / (8 * C_hx * C_hx * C_PrRe) * (7 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[(i + 1) * m + j]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
		gamma / c_coef4 * (7 * Mu(gamma, e_k[a]) + 2 * Mu(gamma, e_k[a + 1]) + 4 * Mu(gamma, e_k[a - 1])) + gamma / c_coef5 * Mu(gamma, e_k[a]);
	A[a][3] = -gamma / c_coef3 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])));
	A[a][4] = -gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[(i + 1) * m + j] - e_k[a]) + (Mu(gamma, e_k[a]) + Mu(gamma, e_k[(i + 1) * m + j])));
	A[a][6] = -gamma / c_coef5 * Mu(gamma, e_k[a]);

	//Для S_qq,C_N/2
	i = qq_i;
	j = cntr_i;
	a = i * m + j;
	A[a][0] = gamma / c_coef1 * (Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[(i - 1) * m + j]) - (Mu(gamma, e_k[(i - 1) * m + j]) + Mu(gamma, e_k[a])));
	A[a][1] = gamma / c_coef3 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1]) + gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a] - e_k[a - 1])
		- gamma / c_coef3 * (Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a])) - gamma / c_coef4 * (2 * Mu(gamma, e_k[a - 1]) + Mu(gamma, e_k[a]));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[a]) / e_k[a] * (6 * e_k[a] - 4 * e_k[(i - 1) * m + j]) -
		gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (6 * e_k[a] - 3 * e_k[a + 1] - 3 * e_k[a - 1]) +
		gamma / (8 * C_hx * C_hx * C_PrRe) * (6 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[(i - 1) * m + j])) +
		gamma / c_coef4 * (6 * Mu(gamma, e_k[a]) + 4 * Mu(gamma, e_k[a + 1]) + 4 * Mu(gamma, e_k[a - 1])) - gamma / (C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[a]);
	A[a][3] = -gamma / c_coef3 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a]) - gamma / c_coef4 * Mu(gamma, e_k[a]) / e_k[a] * (e_k[a + 1] - e_k[a])
		- gamma / c_coef3 * (Mu(gamma, e_k[a]) + Mu(gamma, e_k[a + 1])) - gamma / c_coef4 * (Mu(gamma, e_k[a]) + 2 * Mu(gamma, e_k[a + 1]));
	A[a][7] = gamma / c_coef5 * Mu(gamma, e_k[a]);
	A[a][8] = gamma / c_coef5 * Mu(gamma, e_k[a]);
}

//Вектор правых частей системы линейных уравнений
// qq_i = C_qq
// m = C_M
inline void energy_f(const double gamma, const int qq_i, const int m, double* sigma_k, double* sigma_k1, double* u_k, double* v_k, double* u_k1, double* v_k1, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;

	//Для внутренних узлов
#pragma omp parallel 
	{
#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = 1; i < qq_i; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[(i - 1) * m + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					Mu(gamma, e_k[a]) / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

					Mu(gamma, e_k[a]) / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					Mu(gamma, e_k[a]) / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					Mu(gamma, e_k[a]) / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + C_w; i++)
		{
			for (j = C_cntr + i + 1 - qq_i; j < m - 1; j++)
			{
				a = i * m + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[(i - 1) * m + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					Mu(gamma, e_k[a]) / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

					Mu(gamma, e_k[a]) / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					Mu(gamma, e_k[a]) / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					Mu(gamma, e_k[a]) / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + C_w; i++)
		{
			for (j = C_cntr - i - 1 + qq_i; j > 0; j--)
			{
				a = i * m + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[(i - 1) * m + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					Mu(gamma, e_k[a]) / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

					Mu(gamma, e_k[a]) / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					Mu(gamma, e_k[a]) / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					Mu(gamma, e_k[a]) / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = qq_i + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[(i - 1) * m + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					Mu(gamma, e_k[a]) / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

					Mu(gamma, e_k[a]) / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					Mu(gamma, e_k[a]) / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					Mu(gamma, e_k[a]) / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}
	} // #pragma omp parallel

	//Для Г5. l = C_q-1; m = 1,...,C_q-1;
	i = qq_i + C_w - 1;
#pragma omp parallel
	{
#pragma omp for private(j, a) nowait
		for (j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
		{
			a = i * m + j;
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * ((2 * u_k1[(i + 1) * m + j] - 2 * u_k1[a]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

				Mu(gamma, e_k[a]) / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a])) +

				Mu(gamma, e_k[a]) / (12 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

				Mu(gamma, e_k[a]) / (8 * C_Re * e_k[a]) * (((v_k1[(i + 1) * m + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) * ((v_k1[(i + 1) * m + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) +
					((v_k1[(i + 1) * m + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy) * ((v_k1[(i + 1) * m + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy)) +

				Mu(gamma, e_k[a]) / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		}

		//Для Г6. l = 1,...,C_q-1; m = C_q-1;

#pragma omp for private(i, j, a) nowait
		for (i = qq_i + 1; i < qq_i + C_w - 1; i++)
		{
			j = C_cntr + i - qq_i;
			a = i * m + j;
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a] / C_hy) -
				-P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-u_k1[(i - 1) * m + j] / C_hx + v_k1[a + 1] / C_hy) +

				Mu(gamma, e_k[a]) / (24 * C_hx * C_hx * C_Re * e_k[a]) * (1 * -u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

				Mu(gamma, e_k[a]) / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 1 * v_k1[a] * v_k1[a]) +

				Mu(gamma, e_k[a]) / (16 * C_Re * e_k[a]) * (1 * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy)) +

				Mu(gamma, e_k[a]) / (24 * C_Re * e_k[a]) * (1 * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy));
		}

		//Для Г7.
#pragma omp for private(i, j, a) nowait
		for (i = qq_i + 1; i < qq_i + C_w - 1; i++)
		{
			j = C_cntr - i + qq_i;
			a = i * m + j;
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				-P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-u_k1[(i - 1) * m + j] / C_hx - v_k1[a - 1] / C_hy) +

				Mu(gamma, e_k[a]) / (24 * C_hx * C_hx * C_Re * e_k[a]) * (1 * -u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

				Mu(gamma, e_k[a]) / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + 1 * -v_k1[a] * -v_k1[a]) +

				Mu(gamma, e_k[a]) / (16 * C_Re * e_k[a]) * (1 * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx - u_k1[a] / C_hy)) +

				Mu(gamma, e_k[a]) / (24 * C_Re * e_k[a]) * (1 * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a] / C_hy));
		}
	} // #pragma omp parallel
	//Для S_qq,C_N/2+C_q.
	i = qq_i + C_w - 1;
	j = C_cntr + i - qq_i;
	a = i * m + j;
	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + 2 * v_k1[a + 1] / C_hy - v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
		P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a] / C_hy) +

		Mu(gamma, e_k[a]) / (24 * C_hx * C_hx * C_Re * e_k[a]) * (4 * (u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

		Mu(gamma, e_k[a]) / (24 * C_hy * C_hy * C_Re * e_k[a]) * (4 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 2 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + v_k1[a] * v_k1[a]) +

		Mu(gamma, e_k[a]) / (16 * C_Re * e_k[a]) * (2 * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
			2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
			2 * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
			1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy)) +

		Mu(gamma, e_k[a]) / (24 * C_Re * e_k[a]) * (2 * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
			2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
			2 * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
			1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy));


	//Для S_qq,C_N/2-C_q.
	i = qq_i + C_w - 1;
	j = C_cntr - i + qq_i;
	a = i * m + j;
	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx + v_k1[a + 1] / C_hy + v_k1[a] / C_hy - u_k1[(i - 1) * m + j] / C_hx - 2 * v_k1[a - 1] / C_hy) -
		P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy) +

		Mu(gamma, e_k[a]) / (24 * C_hx * C_hx * C_Re * e_k[a]) * (4 * (u_k1[(i + 1) * m + j] - u_k1[a]) * (u_k1[(i + 1) * m + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j])) +

		Mu(gamma, e_k[a]) / (24 * C_hy * C_hy * C_Re * e_k[a]) * (4 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + 2 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + -v_k1[a] * -v_k1[a]) +

		Mu(gamma, e_k[a]) / (16 * C_Re * e_k[a]) * (2 * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
			1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx - u_k1[a] / C_hy) +
			2 * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
			2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

		Mu(gamma, e_k[a]) / (24 * C_Re * e_k[a]) * (2 * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
			1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx + v_k1[a] / C_hy) +
			2 * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
			2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));


	//Для S_qq,C_N/2
	i = qq_i;
	j = C_cntr;
	a = i * m + j;
	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[a] / C_hx - 2 * u_k1[(i - 1) * m + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) -
		P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-2 * u_k1[a] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) +

		Mu(gamma, e_k[a]) / (12 * C_hx * C_hx * C_Re * e_k[a]) * (2 * (u_k1[a] - u_k1[(i - 1) * m + j]) * (u_k1[a] - u_k1[(i - 1) * m + j]) + -u_k1[a] * -u_k1[a]) +

		Mu(gamma, e_k[a]) / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

		Mu(gamma, e_k[a]) / (16 * C_Re * e_k[a]) * (2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
			2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
			1 * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
			1 * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

		Mu(gamma, e_k[a]) / (24 * C_Re * e_k[a]) * (2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
			2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
			1 * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
			1 * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
inline void energy_d(const int m)
{
	int i = 0;
	int j = 0;
#pragma omp parallel
	{
#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < C_qq + 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				D[i * m + j] = 1 / A[i * m + j][2];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			for (j = C_cntr + i - C_qq; j < m - 1; j++)
			{
				D[i * m + j] = 1 / A[i * m + j][2];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
		{
			for (j = C_cntr - i + C_qq; j > 0; j--)
			{
				D[i * m + j] = 1 / A[i * m + j][2];
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = C_qq + C_w - 1; i < C_M1 - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				D[i * m + j] = 1 / A[i * m + j][2];
			}
		}
	} // #pragma omp parallel
}

//Вектор B = A*Xk1
// m = C_M
// qq_i = C_qq
// w_i = C_w
// cntr_i = C_cntr
// m1_i = C_M1
inline void energy_b(double* e_k1, const int m, const int qq_i, const int w_i, const int cntr_i, const int m1_i)
{
	int i = 0;
	int j = 0;
	int a;
#pragma omp parallel 
	{
		//Для внутренних узлов
#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = 1; i < qq_i; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i + i + 1 - qq_i; j < m - 1; j++)
			{
				a = i * m + j;

				B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i - i - 1 + qq_i; j > 0; j--)
			{
				a = i * m + j;
				B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = qq_i + w_i; i < m1_i - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
			}
		}
	} // #pragma omp parallel 

	//Для Г5. l = C_q-1; m = 1,...,C_q-1;
	i = qq_i + w_i - 1;
#pragma omp parallel 
	{
#pragma omp for private(j, a) nowait
		for (j = cntr_i - C_q + 2; j < cntr_i + C_q - 1; j++)
		{
			a = i * m + j;
			B[a] = A[a][1] * e_k1[a - 1] + A[a][2] * e_k1[i * m + j] + A[a][3] * e_k1[a + 1] +
				A[a][4] * e_k1[(i + 1) * m + j];
		}

		//Для Г6. l = 1,...,C_q-1; m = C_q-1;
#pragma omp for private(i, j, a) nowait
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			j = cntr_i + i - qq_i;
			a = i * m + j;
			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][2] * e_k1[i * m + j] + A[a][3] * e_k1[a + 1]
				+ A[a][5] * e_k1[(i - 1) * m + j - 1] + A[a][8] * e_k1[(i + 1) * m + j + 1];
		}

		//Для Г7.
#pragma omp for private(i, j, a) nowait
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			j = cntr_i - i + qq_i;
			a = i * m + j;
			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + j - 1] + A[a][2] * e_k1[a]
				+ A[a][6] * e_k1[(i - 1) * m + j + 1] + A[a][7] * e_k1[(i + 1) * m + j - 1];
		}
	} // #pragma omp parallel 

	//Для S_qq,C_N/2+C_q.
	i = qq_i + w_i - 1;
	j = cntr_i + i - qq_i;
	a = i * m + j;
	B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
		A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j]
		+ A[a][5] * e_k1[(i - 1) * m + j - 1];

	//Для S_qq,C_N/2-C_q.
	i = qq_i + w_i - 1;
	j = cntr_i - i + qq_i;
	a = i * m + j;
	B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
		A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j]
		+ A[a][6] * e_k1[(i - 1) * m + j + 1];

	//Для S_qq,C_N/2
	i = qq_i;
	j = cntr_i;
	a = i * m + j;
	B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
		A[a][3] * e_k1[i * m + (j + 1)]
		+ A[a][7] * e_k1[(i + 1) * m + j - 1] + A[a][8] * e_k1[(i + 1) * m + j + 1];
}

//Метод Якоби
// m = C_M
// qq_i = C_qq
// w_i = C_w
// m1 = C_M1
inline void energy_jakobi(double* e_k1, double* e2, const int m,
                          const int qq_i,
                          const int w_i, const int m1)
{
	int i = 0;
	int j = 0;
#pragma omp parallel 
	{
#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < qq_i + 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			for (j = C_cntr + i - qq_i; j < m - 1; j++)
			{
				e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			for (j = C_cntr - i + qq_i; j > 0; j--)
			{
				e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = qq_i + w_i - 1; i < m1 - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
			}
		}
	} // #pragma omp parallel 
}

// m = C_M
// n = C_N
// qq_i = C_qq
// w_i = C_w
// m1 = C_M1
// q_i = C_q
inline int energy(const double gamma,
                  double* sigma_k1, double* sigma_k,
                  double* u_k, double* v_k, double* u_k1,
                  double* v_k1, double* e2, double* e_k, double* e_k1,
                  const int m,
                  const int n, const int qq_i,
                  const int w_i, const int m1_i,
                  const int q_i, const int cntr_i)
{
	int i = 0;
	int j = 0;
	int c1;
	int c2;
	int c3;
	int c4;

	energy_a(gamma, m, qq_i, w_i, cntr_i, q_i, sigma_k1, e_k);
	energy_d(m);
	energy_f(gamma, qq_i, m, sigma_k, sigma_k1, u_k, v_k, u_k1, v_k1, e_k);
	int s_e = 0;
	for (s_e = 0; s_e <= 20; ++s_e)
	{
		energy_b(e_k1, m, qq_i, w_i, cntr_i, m1_i);
		energy_jakobi(e_k1, e2, m, qq_i, w_i, m1_i);

		c1 = 0;
		c2 = 0;
		c3 = 0;
		c4 = 0;
#pragma omp parallel
		{
#pragma omp for collapse(2) private(i, j) reduction(+:c1) nowait
			for (i = 1; i < qq_i + 1; i++)
			{
				for (j = 1; j < m - 1; j++)
				{
					if (fabs(e_k1[i * m + j] - e2[i * m + j]) <= C_epsilon)
					{
						++c1;
					}
				}
			}

#pragma omp for private(i, j) reduction(+:c2) nowait
			for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
			{
				for (j = C_cntr + i - qq_i; j < m - 1; j++)
				{
					if (fabs(e_k1[i * m + j] - e2[i * m + j]) <= C_epsilon)
					{
						++c2;
					}
				}
			}

#pragma omp for private(i, j) reduction(+:c3) nowait
			for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
			{
				for (j = C_cntr - i + qq_i; j > 0; j--)
				{
					if (fabs(e_k1[i * m + j] - e2[i * m + j]) <= C_epsilon)
					{
						++c3;
					}
				}
			}

#pragma omp for collapse(2) private(i, j) reduction(+:c4) nowait
			for (i = qq_i + w_i - 1; i < m1_i - 1; i++)
			{
				for (j = 1; j < m - 1; j++)
				{
					if (fabs(e_k1[i * m + j] - e2[i * m + j]) <= C_epsilon)
					{
						++c4;
					}
				}
			}
		} //#pragma omp parallel

		if (c1 + c2 + c3 + c4 == (C_N1 - 1) * (n - 1) - (2 + (q_i - 2 - 1) * 2) / 2 * (q_i - 2))
		{
			break;
		}

#pragma omp parallel
		{
#pragma omp for collapse(2) private(i, j) nowait
			for (i = 1; i < qq_i + 1; i++)
			{
				for (j = 1; j < m - 1; j++)
				{
					e_k1[i * m + j] = e2[i * m + j];
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
			{
				for (j = C_cntr + i - qq_i; j < m - 1; j++)
				{
					e_k1[i * m + j] = e2[i * m + j];
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
			{
				for (j = C_cntr - i + qq_i; j > 0; j--)
				{
					e_k1[i * m + j] = e2[i * m + j];
				}
			}

#pragma omp for collapse(2) private(i, j) nowait
			for (i = qq_i + w_i - 1; i < m1_i - 1; i++)
			{
				for (j = 1; j < m - 1; j++)
				{
					e_k1[i * m + j] = e2[i * m + j];
				}
			}
		} // #pragma omp parallel
	}
	return s_e;
}
