/*----- Функция заполняет элементы матрицы для уравнения энергии----*/
inline void energy_a(double gamma, double* sigma_k1, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;
	//Для внутренних узлов
#pragma omp parallel for collapse(2) private(i, j, a)	
	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
				gamma / (2 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) +
				gamma / (2 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j - 1]));
			A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
			A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr + i + 1 - C_qq;
		a = i * C_M + j;
		A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
		A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1])) - gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])) - gamma / (8 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j - 1]) + 2 * Mu(gamma, e_k[i * C_M + j]));
		A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
			gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
			gamma / (8 * C_hx * C_hx * C_PrRe) * (8 * Mu(gamma, e_k[i * C_M + j]) + 3 * Mu(gamma, e_k[(i + 1) * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
			gamma / (8 * C_hy * C_hy * C_PrRe) * (8 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[i * C_M + j + 1]) + 3 * Mu(gamma, e_k[i * C_M + j - 1]));
		A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
		A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j])) - gamma / (4 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]));
	}

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;
	a = i * C_M + j;
	A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
		gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
		gamma / (2 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) +
		gamma / (2 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j - 1]));
	A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
	A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));


#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr - i - 1 + C_qq;
		a = i * C_M + j;
		A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
		A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
		A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
			gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
			gamma / (8 * C_hx * C_hx * C_PrRe) * (8 * Mu(gamma, e_k[i * C_M + j]) + 3 * Mu(gamma, e_k[(i + 1) * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
			gamma / (8 * C_hy * C_hy * C_PrRe) * (8 * Mu(gamma, e_k[i * C_M + j]) + 3 * Mu(gamma, e_k[i * C_M + j + 1]) + 4 * Mu(gamma, e_k[i * C_M + j - 1]));
		A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j])) - gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])) - gamma / (8 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]));
		A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j])) - gamma / (4 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]));
	}

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;
	a = i * C_M + j;
	A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
		gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
		gamma / (2 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) +
		gamma / (2 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j - 1]));
	A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
	A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
				gamma / (2 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) +
				gamma / (2 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j - 1]));
			A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
			A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 2 + C_qq; j > 0; j--)
		{
			a = i * C_M + j;

			A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
				gamma / (2 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) +
				gamma / (2 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j - 1]));
			A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
			A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));
		}
	}

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / (2 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[(i + 1) * C_M + j] - e_k[(i - 1) * C_M + j]) -
				gamma / (2 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
				gamma / (2 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) +
				gamma / (2 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j - 1]));
			A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
			A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));
		}
	}

	//Для Г5. l = C_q-1; m = 1,...,C_q-1;
	i = C_qq + C_w - 1;

#pragma omp parallel for private(j, a)
	for (j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
	{
		a = i * C_M + j;
		A[a][1] = gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
		A[a][2] = sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - gamma / (4 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - 2 * e_k[(i + 1) * C_M + j]) -
			gamma / (4 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (2 * e_k[i * C_M + j] - e_k[i * C_M + j + 1] - e_k[i * C_M + j - 1]) +
			gamma / (4 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + 2 * Mu(gamma, e_k[(i + 1) * C_M + j])) +
			gamma / (4 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j + 1]));
		A[a][3] = -gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
		A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
	}

	//Для Г6. l = 1,...,C_q-1; m = C_q-1;

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr + i - C_qq;
		a = i * C_M + j;
		A[a][0] = gamma / (4 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j])
			- gamma / (4 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + 2 * Mu(gamma, e_k[(i - 1) * C_M + j]));
		A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (4 * e_k[i * C_M + j] - 3 * e_k[(i - 1) * C_M + j]) -
			gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (4 * e_k[i * C_M + j] - 3 * e_k[i * C_M + j + 1]) +
			gamma / (8 * C_hx * C_hx * C_PrRe) * (4 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
			gamma / (8 * C_hy * C_hy * C_PrRe) * (4 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[i * C_M + j + 1]));
		A[a][3] = -gamma / (4 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) - gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j])
			- gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j])) - gamma / (8 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j + 1]) + Mu(gamma, e_k[i * C_M + j]));
		A[a][5] = -gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
		A[a][8] = gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
	}

	//Для Г7.

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr - i + C_qq;
		a = i * C_M + j;
		A[a][0] = gamma / (4 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j])
			- gamma / (4 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i - 1) * C_M + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + 2 * Mu(gamma, e_k[(i - 1) * C_M + j]));
		A[a][1] = gamma / (4 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) + gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1])
			- gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j - 1])) - gamma / (8 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j - 1]));
		A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (4 * e_k[i * C_M + j] - 3 * e_k[(i - 1) * C_M + j]) -
			gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (4 * e_k[i * C_M + j] - 3 * e_k[i * C_M + j - 1]) +
			gamma / (8 * C_hx * C_hx * C_PrRe) * (4 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
			gamma / (8 * C_hy * C_hy * C_PrRe) * (4 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[i * C_M + j - 1]));
		A[a][6] = -gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
		A[a][7] = gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
	}

	//Для S_qq,C_N/2+C_q.
	i = C_qq + C_w - 1;
	j = C_cntr + i - C_qq;
	a = i * C_M + j;
	A[a][0] = gamma / (4 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j])
		- gamma / (4 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j]));
	A[a][1] = gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (7 * e_k[i * C_M + j] - 4 * e_k[(i + 1) * C_M + j] - 3 * e_k[(i - 1) * C_M + j]) -
		gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (7 * e_k[i * C_M + j] - 4 * e_k[i * C_M + j + 1] - 2 * e_k[i * C_M + j - 1]) +
		gamma / (8 * C_hx * C_hx * C_PrRe) * (7 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[(i + 1) * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
		gamma / (8 * C_hy * C_hy * C_PrRe) * (7 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[i * C_M + j + 1]) + 2 * Mu(gamma, e_k[i * C_M + j - 1])) + gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
	A[a][3] = -gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
	A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));
	A[a][5] = -gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);


	//Для S_qq,C_N/2-C_q.
	i = C_qq + C_w - 1;
	j = C_cntr - i + C_qq;
	a = i * C_M + j;
	A[a][0] = gamma / (4 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) + gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j])
		- gamma / (4 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) - gamma / (8 * C_hx * C_hx * C_PrRe) * (2 * Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j]));
	A[a][1] = gamma / (2 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) - (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (7 * e_k[i * C_M + j] - 4 * e_k[(i + 1) * C_M + j] - 3 * e_k[(i - 1) * C_M + j]) -
		gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (7 * e_k[i * C_M + j] - 2 * e_k[i * C_M + j + 1] - 4 * e_k[i * C_M + j - 1]) +
		gamma / (8 * C_hx * C_hx * C_PrRe) * (7 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[(i + 1) * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
		gamma / (8 * C_hy * C_hy * C_PrRe) * (7 * Mu(gamma, e_k[i * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j + 1]) + 4 * Mu(gamma, e_k[i * C_M + j - 1])) + gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
	A[a][3] = -gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])));
	A[a][4] = -gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[(i + 1) * C_M + j] - e_k[i * C_M + j]) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])));
	A[a][6] = -gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);

	//Для S_qq,C_N/2
	i = C_qq;
	j = C_cntr;
	a = i * C_M + j;
	A[a][0] = gamma / (2 * C_hx * C_hx * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[(i - 1) * C_M + j]) - (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])));
	A[a][1] = gamma / (4 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1]) + gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j] - e_k[i * C_M + j - 1])
		- gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j])) - gamma / (8 * C_hy * C_hy * C_PrRe) * (2 * Mu(gamma, e_k[i * C_M + j - 1]) + Mu(gamma, e_k[i * C_M + j]));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - gamma / (8 * C_hx * C_hx * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (6 * e_k[i * C_M + j] - 4 * e_k[(i - 1) * C_M + j]) -
		gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (6 * e_k[i * C_M + j] - 3 * e_k[i * C_M + j + 1] - 3 * e_k[i * C_M + j - 1]) +
		gamma / (8 * C_hx * C_hx * C_PrRe) * (6 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[(i - 1) * C_M + j])) +
		gamma / (8 * C_hy * C_hy * C_PrRe) * (6 * Mu(gamma, e_k[i * C_M + j]) + 4 * Mu(gamma, e_k[i * C_M + j + 1]) + 4 * Mu(gamma, e_k[i * C_M + j - 1])) - gamma / (C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
	A[a][3] = -gamma / (4 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j]) - gamma / (8 * C_hy * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]) / e_k[i * C_M + j] * (e_k[i * C_M + j + 1] - e_k[i * C_M + j])
		- gamma / (4 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + j + 1])) - gamma / (8 * C_hy * C_hy * C_PrRe) * (Mu(gamma, e_k[i * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j + 1]));
	A[a][7] = gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
	A[a][8] = gamma / (2 * C_tg * C_hx * C_hy * C_PrRe) * Mu(gamma, e_k[i * C_M + j]);
}

//Вектор правых частей системы линейных уравнений
inline void energy_f(double gamma, double* sigma_k, double* sigma_k1, double* u_k, double* v_k, double* u_k1, double* v_k1, double* e_k)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (4 * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j - 1]) / C_hy) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) + (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * ((v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1])) +

				Mu(gamma, e_k[i * C_M + j]) / (8 * C_Re * e_k[i * C_M + j]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy)) +

				Mu(gamma, e_k[i * C_M + j]) / (12 * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (4 * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j - 1]) / C_hy) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) + (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * ((v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1])) +

				Mu(gamma, e_k[i * C_M + j]) / (8 * C_Re * e_k[i * C_M + j]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy)) +

				Mu(gamma, e_k[i * C_M + j]) / (12 * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = i * C_M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (4 * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j - 1]) / C_hy) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) + (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * ((v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1])) +

				Mu(gamma, e_k[i * C_M + j]) / (8 * C_Re * e_k[i * C_M + j]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy)) +

				Mu(gamma, e_k[i * C_M + j]) / (12 * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));
		}
	}

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (4 * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[(i - 1) * C_M + j]) / C_hx + (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j - 1]) / C_hy) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) + (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

				Mu(gamma, e_k[i * C_M + j]) / (6 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * ((v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1])) +

				Mu(gamma, e_k[i * C_M + j]) / (8 * C_Re * e_k[i * C_M + j]) * ((v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
					(v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
					(v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy)) +

				Mu(gamma, e_k[i * C_M + j]) / (12 * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
					(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
					(u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));
		}
	}

	//Для Г5. l = C_q-1; m = 1,...,C_q-1;
	i = C_qq + C_w - 1;
#pragma omp parallel for private(j, a)
	for (j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
	{
		a = i * C_M + j;
		f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (8 * e_k[i * C_M + j]) * ((2 * u_k1[(i + 1) * C_M + j] - 2 * u_k1[i * C_M + j]) / C_hx + (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j - 1]) / C_hy) +

			Mu(gamma, e_k[i * C_M + j]) / (6 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j])) +

			Mu(gamma, e_k[i * C_M + j]) / (12 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * ((v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1])) +

			Mu(gamma, e_k[i * C_M + j]) / (8 * C_Re * e_k[i * C_M + j]) * (((v_k1[(i + 1) * C_M + j] - v_k1[i * C_M + j]) / C_hx + (u_k1[i * C_M + j + 1] - u_k1[i * C_M + j]) / C_hy) * ((v_k1[(i + 1) * C_M + j] - v_k1[i * C_M + j]) / C_hx + (u_k1[i * C_M + j + 1] - u_k1[i * C_M + j]) / C_hy) +
				((v_k1[(i + 1) * C_M + j] - v_k1[i * C_M + j]) / C_hx + (u_k1[i * C_M + j] - u_k1[i * C_M + j - 1]) / C_hy) * ((v_k1[(i + 1) * C_M + j] - v_k1[i * C_M + j]) / C_hx + (u_k1[i * C_M + j] - u_k1[i * C_M + j - 1]) / C_hy)) +

			Mu(gamma, e_k[i * C_M + j]) / (12 * C_Re * e_k[i * C_M + j]) * ((u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
				(u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));
	}

	//Для Г6. l = 1,...,C_q-1; m = C_q-1;

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr + i - C_qq;
		a = i * C_M + j;
		f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (8 * e_k[i * C_M + j]) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j + 1] / C_hy - v_k1[i * C_M + j] / C_hy) -
			-P(gamma, sigma_k[a], e_k[i * C_M + j]) / (16 * e_k[i * C_M + j]) * (-u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j + 1] / C_hy) +

			Mu(gamma, e_k[i * C_M + j]) / (24 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * (1 * -u_k1[i * C_M + j] * -u_k1[i * C_M + j] + 3 * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

			Mu(gamma, e_k[i * C_M + j]) / (24 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * (3 * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + 1 * v_k1[i * C_M + j] * v_k1[i * C_M + j]) +

			Mu(gamma, e_k[i * C_M + j]) / (16 * C_Re * e_k[i * C_M + j]) * (1 * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
				2 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
				1 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy)) +

			Mu(gamma, e_k[i * C_M + j]) / (24 * C_Re * e_k[i * C_M + j]) * (1 * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
				2 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
				1 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy));
	}

	//Для Г7.
#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr - i + C_qq;
		a = i * C_M + j;
		f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (8 * e_k[i * C_M + j]) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j] / C_hy - v_k1[i * C_M + j - 1] / C_hy) -
			-P(gamma, sigma_k[a], e_k[i * C_M + j]) / (16 * e_k[i * C_M + j]) * (-u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j - 1] / C_hy) +

			Mu(gamma, e_k[i * C_M + j]) / (24 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * (1 * -u_k1[i * C_M + j] * -u_k1[i * C_M + j] + 3 * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

			Mu(gamma, e_k[i * C_M + j]) / (24 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * (3 * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) + 1 * -v_k1[i * C_M + j] * -v_k1[i * C_M + j]) +

			Mu(gamma, e_k[i * C_M + j]) / (16 * C_Re * e_k[i * C_M + j]) * (1 * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
				2 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
				1 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hy)) +

			Mu(gamma, e_k[i * C_M + j]) / (24 * C_Re * e_k[i * C_M + j]) * (1 * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
				2 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
				1 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j] / C_hy));
	}

	//Для S_qq,C_N/2+C_q.
	i = C_qq + C_w - 1;
	j = C_cntr + i - C_qq;
	a = i * C_M + j;
	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (8 * e_k[i * C_M + j]) * (2 * u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + 2 * v_k1[i * C_M + j + 1] / C_hy - v_k1[i * C_M + j] / C_hy - v_k1[i * C_M + j - 1] / C_hy) -
		P(gamma, sigma_k[a], e_k[i * C_M + j]) / (16 * e_k[i * C_M + j]) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j] / C_hy) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * (4 * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) + 3 * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * (4 * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + 2 * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) + v_k1[i * C_M + j] * v_k1[i * C_M + j]) +

		Mu(gamma, e_k[i * C_M + j]) / (16 * C_Re * e_k[i * C_M + j]) * (2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
			2 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
			2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
			1 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy)) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_Re * e_k[i * C_M + j]) * (2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
			2 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
			2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
			1 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy));


	//Для S_qq,C_N/2-C_q.
	i = C_qq + C_w - 1;
	j = C_cntr - i + C_qq;
	a = i * C_M + j;
	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (8 * e_k[i * C_M + j]) * (2 * u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx + v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy - u_k1[(i - 1) * C_M + j] / C_hx - 2 * v_k1[i * C_M + j - 1] / C_hy) -
		P(gamma, sigma_k[a], e_k[i * C_M + j]) / (16 * e_k[i * C_M + j]) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * (4 * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) * (u_k1[(i + 1) * C_M + j] - u_k1[i * C_M + j]) + 3 * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j])) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * (4 * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) + 2 * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + -v_k1[i * C_M + j] * -v_k1[i * C_M + j]) +

		Mu(gamma, e_k[i * C_M + j]) / (16 * C_Re * e_k[i * C_M + j]) * (2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
			1 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hy) +
			2 * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[(i + 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
			2 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy)) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_Re * e_k[i * C_M + j]) * (2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
			1 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j] / C_hy) +
			2 * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[(i + 1) * C_M + j] / C_hx - u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
			2 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));


	//Для S_qq,C_N/2
	i = C_qq;
	j = C_cntr;
	a = i * C_M + j;
	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - P(gamma, sigma_k[a], e_k[i * C_M + j]) / (8 * e_k[i * C_M + j]) * (2 * u_k1[i * C_M + j] / C_hx - 2 * u_k1[(i - 1) * C_M + j] / C_hx + v_k1[i * C_M + j + 1] / C_hy - v_k1[i * C_M + j - 1] / C_hy) -
		P(gamma, sigma_k[a], e_k[i * C_M + j]) / (16 * e_k[i * C_M + j]) * (-2 * u_k1[i * C_M + j] / C_hx + v_k1[i * C_M + j + 1] / C_hy - v_k1[i * C_M + j - 1] / C_hy) +

		Mu(gamma, e_k[i * C_M + j]) / (12 * C_hx * C_hx * C_Re * e_k[i * C_M + j]) * (2 * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) * (u_k1[i * C_M + j] - u_k1[(i - 1) * C_M + j]) + -u_k1[i * C_M + j] * -u_k1[i * C_M + j]) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_hy * C_hy * C_Re * e_k[i * C_M + j]) * (3 * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) * (v_k1[i * C_M + j + 1] - v_k1[i * C_M + j]) + 3 * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1]) * (v_k1[i * C_M + j] - v_k1[i * C_M + j - 1])) +

		Mu(gamma, e_k[i * C_M + j]) / (16 * C_Re * e_k[i * C_M + j]) * (2 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) +
			2 * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (v_k1[i * C_M + j] / C_hx - v_k1[(i - 1) * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
			1 * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j + 1] / C_hy - u_k1[i * C_M + j] / C_hy) +
			1 * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy) * (-v_k1[i * C_M + j] / C_hx + u_k1[i * C_M + j] / C_hy - u_k1[i * C_M + j - 1] / C_hy)) +

		Mu(gamma, e_k[i * C_M + j]) / (24 * C_Re * e_k[i * C_M + j]) * (2 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) +
			2 * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (u_k1[i * C_M + j] / C_hx - u_k1[(i - 1) * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
			1 * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j + 1] / C_hy + v_k1[i * C_M + j] / C_hy) +
			1 * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy) * (-u_k1[i * C_M + j] / C_hx - v_k1[i * C_M + j] / C_hy + v_k1[i * C_M + j - 1] / C_hy));
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
inline void energy_d()
{
	int i = 0;
	int j = 0;
	int a;

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = 1; i < C_qq + 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		for (j = C_cntr + i - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		for (j = C_cntr - i + C_qq; j > 0; j--)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = C_qq + C_w - 1; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}
}

//Вектор B = A*Xk1
// m = C_M
inline void energy_b(double* e_k1, const int m)
{
	int i = 0;
	int j = 0;
	int a;

	//Для внутренних узлов

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < m - 1; j++)
		{
			a = i * m + j;
			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < m - 1; j++)
		{
			a = i * m + j;

			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = i * m + j;

			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
		}
	}

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < m - 1; j++)
		{
			a = i * m + j;

			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
				A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
		}
	}

	//Для Г5. l = C_q-1; m = 1,...,C_q-1;
	i = C_qq + C_w - 1;
#pragma omp parallel for private(j, a)
	for (j = C_cntr - C_q + 2; j < C_cntr + C_q - 1; j++)
	{
		a = i * m + j;
		B[a] = A[a][1] * e_k1[a - 1] + A[a][2] * e_k1[i * m + j] + A[a][3] * e_k1[a + 1] +
			A[a][4] * e_k1[(i + 1) * m + j];
	}

	//Для Г6. l = 1,...,C_q-1; m = C_q-1;

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr + i - C_qq;
		a = i * m + j;
		B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][2] * e_k1[i * m + j] + A[a][3] * e_k1[a + 1]
			+ A[a][5] * e_k1[(i - 1) * m + j - 1] + A[a][8] * e_k1[(i + 1) * m + j + 1];
	}

	//Для Г7.

#pragma omp parallel for private(i, j, a)
	for (i = C_qq + 1; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr - i + C_qq;
		a = i * m + j;
		B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + j - 1] + A[a][2] * e_k1[a]
			+ A[a][6] * e_k1[(i - 1) * m + j + 1] + A[a][7] * e_k1[(i + 1) * m + j - 1];
	}

	//Для S_qq,C_N/2+C_q.
	i = C_qq + C_w - 1;
	j = C_cntr + i - C_qq;
	a = i * m + j;
	B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
		A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j]
		+ A[a][5] * e_k1[(i - 1) * m + j - 1];

	//Для S_qq,C_N/2-C_q.
	i = C_qq + C_w - 1;
	j = C_cntr - i + C_qq;
	a = i * m + j;
	B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
		A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j]
		+ A[a][6] * e_k1[(i - 1) * m + j + 1];

	//Для S_qq,C_N/2
	i = C_qq;
	j = C_cntr;
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
inline double energy_jakobi(double* e_k1, double* e2, const int m,
                            const int qq_i,
                            const int w_i, const int m1)
{
	int i = 0;
	int j = 0;
	int a;

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = 1; i < qq_i + 1; i++)
	{
		for (j = 1; j < m - 1; j++)
		{
			a = i * m + j;
			e2[a] = e_k1[a] - D[a] * (B[a] - f[a]);
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
	{
		for (j = C_cntr + i - qq_i; j < m - 1; j++)
		{
			a = i * m + j;
			e2[a] = e_k1[a] - D[a] * (B[a] - f[a]);
		}
	}

#pragma omp parallel for private(i, j, a)
	for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
	{
		for (j = C_cntr - i + qq_i; j > 0; j--)
		{
			a = i * m + j;
			e2[a] = e_k1[a] - D[a] * (B[a] - f[a]);
		}
	}

#pragma omp parallel for collapse(2) private(i, j, a)
	for (i = qq_i + w_i - 1; i < m1 - 1; i++)
	{
		for (j = 1; j < m - 1; j++)
		{
			a = i * m + j;
			e2[a] = e_k1[a] - D[a] * (B[a] - f[a]);
		}
	}

	return 0;
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
                  const int w_i, const int m1,
                  const int q_i)
{
	int i = 0;
	int j = 0;
	int bl = 1;
	int a;
	int c;
	int s_e;

	energy_a(gamma, sigma_k1, e_k);
	energy_d();
	energy_f(gamma, sigma_k, sigma_k1, u_k, v_k, u_k1, v_k1, e_k);

	s_e = 0;
	while (bl)
	{
		energy_b(e_k1, m);
		energy_jakobi(e_k1, e2, m, qq_i, w_i, m1);

		c = 0;

#pragma omp parallel for collapse(2) private(i, j, a)
		for (i = 1; i < qq_i + 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				if (fabs(e_k1[a] - e2[a]) <= C_epsilon)
				{
#pragma omp atomic
					c++;
				}
			}
		}

#pragma omp parallel for private(i, j, a)
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			for (j = C_cntr + i - qq_i; j < m - 1; j++)
			{
				a = i * m + j;
				if (fabs(e_k1[a] - e2[a]) <= C_epsilon)
				{
#pragma omp atomic
					c++;
				}
			}
		}

#pragma omp parallel for private(i, j, a)
		for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			for (j = C_cntr - i + qq_i; j > 0; j--)
			{
				a = i * m + j;
				if (fabs(e_k1[a] - e2[a]) <= C_epsilon)
				{
#pragma omp atomic
					c++;
				}
			}
		}

#pragma omp parallel for collapse(2) private(i, j, a)
		for (i = qq_i + w_i - 1; i < m1 - 1; i++)
		{
			for (j = 1; j < m - 1; j++)
			{
				a = i * m + j;
				if (fabs(e_k1[a] - e2[a]) <= C_epsilon)
				{
#pragma omp atomic
					c++;
				}
			}
		}

		if (c == (C_N1 - 1) * (n - 1) - (2 + (q_i - 2 - 1) * 2) / 2 * (q_i - 2))
		{
			bl = 0;
		}
		else if (s_e > 20)
		{
			bl = 0;
		}

		else
		{
#pragma omp parallel for collapse(2) private(i, j, a)
			for (i = 1; i < qq_i + 1; i++)
			{
				for (j = 1; j < m - 1; j++)
				{
					a = i * m + j;
					e_k1[a] = e2[a];
				}
			}

#pragma omp parallel for private(i, j, a)
			for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
			{
				for (j = C_cntr + i - qq_i; j < m - 1; j++)
				{
					a = i * m + j;
					e_k1[a] = e2[a];
				}
			}

#pragma omp parallel for private(i, j, a)		
			for (i = qq_i + 1; i < qq_i + w_i - 1; i++)
			{
				for (j = C_cntr - i + qq_i; j > 0; j--)
				{
					a = i * m + j;
					e_k1[a] = e2[a];
				}
			}

#pragma omp parallel for collapse(2) private(i, j, a)
			for (i = qq_i + w_i - 1; i < m1 - 1; i++)
			{
				for (j = 1; j < m - 1; j++)
				{
					a = i * m + j;
					e_k1[a] = e2[a];
				}
			}
		}
		s_e++;
	}
	return s_e;
}
