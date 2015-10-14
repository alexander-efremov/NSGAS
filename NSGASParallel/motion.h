/*----- Функция заполняет элементы матрицы, составленной для двух уравнении движения.----*/
inline void motion_a(const double gamma, double* sigma_k1, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;

	const double c_coef = (C_hx * C_hy * C_Re);

	///////////////////////////////////////////////////Уравнение для u
	//Для внутренних узлов.
#pragma omp parallel 
	{
#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][4] = -2 * (gamma , Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr + i + 1 - C_qq;
			a = i * C_M + j;
			A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
			A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (4 * C_hy * C_hy * C_Re) - (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (8 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re)
				+ (Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[(i + 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (6 * C_hx * C_hx * C_Re)
				+ (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j - 1)])) / (4 * C_hy * C_hy * C_Re) + (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j - 1)])) / (8 * C_hy * C_hy * C_Re);
			A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
			A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (6 * C_hx * C_hx * C_Re);
			A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
			A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
			A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;
			A[a][9] = (1. / 4. - 1. / 8.) * Mu(gamma, e_k[i * C_M + (j - 1)]) / c_coef;
			A[a][11] = (1. / 12. - 1. / 6.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / c_coef;
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
			A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re)
				+ (Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[(i + 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (6 * C_hx * C_hx * C_Re)
				+ (Mu(gamma, e_k[i * C_M + (j + 1)]) + Mu(gamma, e_k[i * C_M + j])) / (4 * C_hy * C_hy * C_Re) + (Mu(gamma, e_k[i * C_M + (j + 1)]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (8 * C_hy * C_hy * C_Re);
			A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (4 * C_hy * C_hy * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (8 * C_hy * C_hy * C_Re);
			A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (6 * C_hx * C_hx * C_Re);
			A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
			A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
			A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
			A[a][10] = (1. / 8. - 1. / 4.) * Mu(gamma, e_k[i * C_M + (j + 1)]) / c_coef;
			A[a][11] = (1. / 6. - 1. / 12.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / c_coef;
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				a = i * C_M + j;
				A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
				A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;
			}
		}
	} // #pragma omp parallel

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;
	a = i * C_M + j;
	A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
	A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
	A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
	A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;
	a = i * C_M + j;
	A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / c_coef;
	A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / c_coef;
	A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / c_coef;
	A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / c_coef;

	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов. l,m = 1,...,n-1
#pragma omp parallel 
	{
#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = C_M2 + i * C_M + j;
				A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr + i + 1 - C_qq;
			a = C_M2 + i * C_M + j;
			A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
			A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re) - (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (6 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re)
				+ (Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (4 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[(i + 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (8 * C_hx * C_hx * C_Re)
				+ (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j - 1)])) / (3 * C_hy * C_hy * C_Re) + (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j - 1)])) / (6 * C_hy * C_hy * C_Re);
			A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
			A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (4 * C_hx * C_hx * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (8 * C_hx * C_hx * C_Re);
			A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
			A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
			A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;
			A[a][9] = (1. / 12. - 1. / 6.) * Mu(gamma, e_k[i * C_M + (j - 1)]) / c_coef;
			A[a][11] = (1. / 4. - 1. / 8.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / c_coef;
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr - i - 1 + C_qq;
			a = C_M2 + i * C_M + j;
			A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
			A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re)
				+ (Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (4 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[(i + 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (8 * C_hx * C_hx * C_Re)
				+ (Mu(gamma, e_k[i * C_M + (j + 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re) + (Mu(gamma, e_k[i * C_M + (j + 1)]) + 2 * Mu(gamma, e_k[i * C_M + j])) / (6 * C_hy * C_hy * C_Re);
			A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (6 * C_hy * C_hy * C_Re);
			A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (4 * C_hx * C_hx * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (8 * C_hx * C_hx * C_Re);
			A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
			A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
			A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
			A[a][10] = (1. / 6. - 1. / 12.) * Mu(gamma, e_k[i * C_M + (j + 1)]) / c_coef;
			A[a][11] = (1. / 8. - 1. / 4.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / c_coef;
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				a = C_M2 + i * C_M + j;
				A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				a = C_M2 + i * C_M + j;
				A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = C_M2 + i * C_M + j;
				A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
				A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
				A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
				A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
				A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
				A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
				A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;
			}
		}
	} // #pragma omp parallel

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;
	a = C_M2 + i * C_M + j;
	A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
	A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
	A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
	A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;
	a = C_M2 + i * C_M + j;
	A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / c_coef;
	A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / c_coef;
	A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / c_coef;
	A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / c_coef;
}


//Вектор правых частей системы уравнений
inline void motion_f(double gamma, double* sigma_k, double* sigma_k1, double* u_k, double* v_k, double* e_k)
{
	int i = 0;
	int j = 0;

	///////////////////////////////////////////////////Уравнение для u
	//Для внутренних узлов. l,m = 1,...,n-1
#pragma omp parallel
	{
#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				f[i * C_M + j] = uX_k[i * C_M + j] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				f[i * C_M + j] = uX_k[i * C_M + j] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				f[i * C_M + j] = uX_k[i * C_M + j] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				f[i * C_M + j] = uX_k[i * C_M + j] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
			}
		}
	} // #pragma omp parallel
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
// qq_i = C_qq
// m_i = C_M
// cntr_i = C_cnrt
// w_i = C_w
// m2_i = C_M2
inline void motion_d(const int qq_i, const int m_i, const int cntr_i, const int w_i, const int m2_i)
{
	int i = 0;
	int j = 0;

#pragma omp parallel
	{
#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < qq_i; i++)
		{
			for (j = 1; j < m_i - 1; j++)
			{
				D[i * m_i + j] = 1 / A[i * m_i + j][2];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
			{
				D[i * m_i + j] = 1 / A[i * m_i + j][2];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i - i - 1 + qq_i; j > 0; j--)
			{
				D[i * m_i + j] = 1 / A[i * m_i + j][2];
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = qq_i + w_i; i < C_M1 - 1; i++)
		{
			for (j = 1; j < m_i - 1; j++)
			{
				D[i * m_i + j] = 1 / A[i * m_i + j][2];
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < qq_i; i++)
		{
			for (j = 1; j < m_i - 1; j++)
			{
				D[m2_i + i * m_i + j] = 1 / A[m2_i + i * m_i + j][2];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
			{
				D[m2_i + i * m_i + j] = 1 / A[m2_i + i * m_i + j][2];
			}
		}

#pragma omp for private(i, j) nowait
		for (i = qq_i; i < qq_i + w_i; i++)
		{
			for (j = cntr_i - i - 1 + qq_i; j > 0; j--)
			{
				D[m2_i + i * m_i + j] = 1 / A[m2_i + i * m_i + j][2];
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = qq_i + w_i; i < C_M1 - 1; i++)
		{
			for (j = 1; j < m_i - 1; j++)
			{
				D[m2_i + i * m_i + j] = 1 / A[m2_i + i * m_i + j][2];
			}
		}
	} // #pragma omp parallel
}

//Вектор B = A*Xk1
inline void motion_b(double* u_k1, double* v_k1)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u
	//Для внутренних узлов
#pragma omp parallel
	{
#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr + i + 1 - C_qq;
			a = i * C_M + j;
			B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)] +
				A[a][9] * v_k1[i * C_M + (j - 1)] +
				A[a][11] * v_k1[(i + 1) * C_M + j];
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr - i - 1 + C_qq;
			a = i * C_M + j;
			B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][10] * v_k1[i * C_M + (j + 1)] +
				A[a][11] * v_k1[(i + 1) * C_M + j];
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				a = i * C_M + j;

				B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];
			}
		}
	} // #pragma omp parallel

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;
	a = i * C_M + j;
	B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
		A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
		A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
		A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
		A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
		A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;
	a = i * C_M + j;
	B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
		A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
		A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
		A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
		A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
		A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

	//Уравнение для v

	//Для внутренних узлов
#pragma omp parallel
	{
#pragma omp for private(i, j, a) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = C_M2 + i * C_M + j;
				B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * u_k1[(i + 1) * C_M + j + 1];
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr + i + 1 - C_qq;
			a = C_M2 + i * C_M + j;
			B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
				A[a][8] * u_k1[(i + 1) * C_M + j + 1] +
				A[a][9] * u_k1[i * C_M + j - 1] +
				A[a][11] * u_k1[(i + 1) * C_M + j];
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w - 1; i++)
		{
			j = C_cntr - i - 1 + C_qq;
			a = C_M2 + i * C_M + j;
			B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
				A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
				A[a][10] * u_k1[i * C_M + j + 1] +
				A[a][11] * u_k1[(i + 1) * C_M + j];
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 2 - C_qq; j < C_M - 1; j++)
			{
				a = C_M2 + i * C_M + j;
				B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * u_k1[(i + 1) * C_M + j + 1];
			}
		}

#pragma omp for private(i, j, a) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 2 + C_qq; j > 0; j--)
			{
				a = C_M2 + i * C_M + j;
				B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * u_k1[(i + 1) * C_M + j + 1];
			}
		}

#pragma omp for collapse(2) private(i, j, a) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = C_M2 + i * C_M + j;
				B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * u_k1[(i + 1) * C_M + j + 1];
			}
		}
	} // #pragma omp parallel

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;
	a = C_M2 + i * C_M + j;
	B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
		A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
		A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
		A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
		A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
		A[a][8] * u_k1[(i + 1) * C_M + j + 1];
	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;
	a = C_M2 + i * C_M + j;
	B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
		A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
		A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
		A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
		A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
		A[a][8] * u_k1[(i + 1) * C_M + j + 1];
}

//Метод Якоби
inline void motion_jakobi(double* u2, double* u_k1, double* v2, double* v_k1)
{
	int i = 0;
	int j = 0;
#pragma omp parallel
	{
#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{				
				u2[i * C_M + j] = u_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{				
				u2[i * C_M + j] = u_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{				
				u2[i * C_M + j] = u_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{				
				u2[i * C_M + j] = u_k1[i * C_M + j] - D[i * C_M + j] * (B[i * C_M + j] - f[i * C_M + j]);
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				v2[i * C_M + j] = v_k1[i * C_M + j] - D[C_M2 + i * C_M + j] * (B[C_M2 + i * C_M + j] - f[C_M2 + i * C_M + j]);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				v2[i * C_M + j] = v_k1[i * C_M + j] - D[C_M2 + i * C_M + j] * (B[C_M2 + i * C_M + j] - f[C_M2 + i * C_M + j]);
			}
		}

#pragma omp for private(i, j) nowait
		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				v2[i * C_M + j] = v_k1[i * C_M + j] - D[C_M2 + i * C_M + j] * (B[C_M2 + i * C_M + j] - f[C_M2 + i * C_M + j]);
			}
		}

#pragma omp for collapse(2) private(i, j) nowait
		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				v2[i * C_M + j] = v_k1[i * C_M + j] - D[C_M2 + i * C_M + j] * (B[C_M2 + i * C_M + j] - f[C_M2 + i * C_M + j]);
			}
		}
	} // #pragma omp parallel
}

// m_i = C_M
// m1_i = C_M1
// qq_i = C_qq
// w_i = C_w
// cntr_i = C_cntr
inline int motion(const double gamma, const int m_i, const int m1_i, const int m2_i, const int qq_i,
                  const int w_i, const int cntr_i, double* sigma_k1, double* sigma_k,
                  double* u_k, double* v_k, double* u_k1, double* v_k1, double* u2, double* v2, double* e_k)
{
	int i = 0;
	int j = 0;
	int c_u1;
	int c_u2;
	int c_u3;
	int c_u4;
	int c_v1;
	int c_v2;
	int c_v3;
	int c_v4;

	/*---------------------------------------------*/

	motion_a(gamma, sigma_k1, e_k);
	motion_d(qq_i,	m_i, cntr_i, w_i, m2_i);
	motion_f(gamma, sigma_k, sigma_k1, u_k, v_k, e_k);

	int s_m = 0;
	for (s_m = 0; s_m <= 20; ++s_m)
	{
		motion_b(u_k1, v_k1);
		motion_jakobi(u2, u_k1, v2, v_k1);

		c_u1 = 0;
		c_u2 = 0;
		c_u3 = 0;
		c_u4 = 0;
		c_v1 = 0;
		c_v2 = 0;
		c_v3 = 0;
		c_v4 = 0;
#pragma omp parallel
		{
#pragma omp for collapse(2) private(i, j) reduction(+:c_u1, c_v1) nowait
			for (i = 1; i < qq_i; i++)
			{
				for (j = 1; j < m_i - 1; j++)
				{
					if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= C_epsilon)
					{
						++c_u1;
					}
					if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= C_epsilon)
					{
						++c_v1;
					}
				}
			}

#pragma omp for private(i, j) reduction(+:c_u2, c_v2) nowait
			for (i = qq_i; i < qq_i + w_i; i++)
			{
				for (j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
				{
					if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= C_epsilon)
					{
						++c_u2;
					}
					if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= C_epsilon)
					{
						++c_v2;
					}
				}
			}

#pragma omp for private(i, j) reduction(+:c_u3, c_v3) nowait
			for (i = qq_i; i < qq_i + w_i; i++)
			{
				for (j = cntr_i - i - 1 + qq_i; j > 0; j--)
				{
					if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= C_epsilon)
					{
						++c_u3;
					}
					if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= C_epsilon)
					{
						++c_v3;
					}
				}
			}

#pragma omp for collapse(2) private(i, j) reduction(+:c_u4, c_v4) nowait
			for (i = qq_i + w_i; i < m1_i - 1; i++)
			{
				for (j = 1; j < m_i - 1; j++)
				{
					if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= C_epsilon)
					{
						++c_u4;
					}
					if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= C_epsilon)
					{
						++c_v4;
					}
				}
			}
		} // #pragma omp parallel

		if (c_u1 + c_u2 + c_u3 + c_u4 == (C_N1 - 1) * (C_N - 1) - (2 + (C_q - 1) * 2) / 2 * C_q && c_v1 + c_v2 + c_v3 + c_v4 >= (C_N1 - 1) * (C_N - 1) - (2 + (C_q - 1) * 2) / 2 * C_q)
		{
			break;
		}
#pragma omp parallel
		{
#pragma omp for collapse(2) private(i, j) nowait
			for (i = 1; i < qq_i; i++)
			{
				for (j = 1; j < m_i - 1; j++)
				{
					u_k1[i * m_i + j] = u2[i * m_i + j];
					v_k1[i * m_i + j] = v2[i * m_i + j];
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i; i < qq_i + w_i; i++)
			{
				for (j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
				{
					u_k1[i * m_i + j] = u2[i * m_i + j];
					v_k1[i * m_i + j] = v2[i * m_i + j];
				}
			}

#pragma omp for private(i, j) nowait
			for (i = qq_i; i < qq_i + w_i; i++)
			{
				for (j = cntr_i - i - 1 + qq_i; j > 0; j--)
				{
					u_k1[i * m_i + j] = u2[i * m_i + j];
					v_k1[i * m_i + j] = v2[i * m_i + j];
				}
			}

#pragma omp for collapse(2) private(i, j) nowait
			for (i = qq_i + w_i; i < m1_i - 1; i++)
			{
				for (j = 1; j < m_i - 1; j++)
				{
					u_k1[i * m_i + j] = u2[i * m_i + j];
					v_k1[i * m_i + j] = v2[i * m_i + j];
				}
			}
		}
	} // #pragma omp parallel
	return s_m;
}
