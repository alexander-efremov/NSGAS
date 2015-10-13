/*----- Функция заполняет элементы матрицы, составленной для двух уравнении движения.----*/

inline double motion_a(double gamma, double* sigma_k1, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u

//#pragma omp  parallel for private(i, j, a)

	//Для внутренних узлов.
	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
			A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
			A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
			A[a][4] = -2 * (gamma, Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
			A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);
		}
	}


//#pragma omp  parallel for private(i, j, a)

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
		A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
		A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
		A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);
		A[a][9] = (1. / 4. - 1. / 8.) * Mu(gamma, e_k[i * C_M + (j - 1)]) / (C_hx * C_hy * C_Re);
		A[a][11] = (1. / 12. - 1. / 6.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / (C_hx * C_hy * C_Re);
	}

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;

	a = i * C_M + j;

	A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
	A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
	A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
	A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr - i - 1 + C_qq;

		a = i * C_M + j;

		A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
		A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
		A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re)
			+ (Mu(gamma, e_k[(i + 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[(i + 1) * C_M + j]) + 2 * Mu(gamma,e_k[i * C_M + j])) / (6 * C_hx * C_hx * C_Re)
			+ (Mu(gamma, e_k[i * C_M + (j + 1)]) + Mu(gamma, e_k[i * C_M + j])) / (4 * C_hy * C_hy * C_Re) + (Mu(gamma, e_k[i * C_M + (j + 1)]) + 2 * Mu(gamma,e_k[i * C_M + j])) / (8 * C_hy * C_hy * C_Re);
		A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (4 * C_hy * C_hy * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma,e_k[i * C_M + (j + 1)])) / (8 * C_hy * C_hy * C_Re);
		A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) - (2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma,e_k[(i + 1) * C_M + j])) / (6 * C_hx * C_hx * C_Re);
		A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
		A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
		A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
		A[a][10] = (1. / 8. - 1. / 4.) * Mu(gamma, e_k[i * C_M + (j + 1)]) / (C_hx * C_hy * C_Re);
		A[a][11] = (1. / 6. - 1. / 12.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / (C_hx * C_hy * C_Re);
	}

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;
	a = i * C_M + j;

	A[a][0] = -2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][1] = -(Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau + 2 * (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re) + (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][3] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (2 * C_hy * C_hy * C_Re);
	A[a][4] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (3 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
	A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
	A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
	A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);


//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j - 1)]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[i * C_M + j + 1]) / 4 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[i * C_M + j - 1]) / 4 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 6 - Mu(gamma, e_k[i * C_M + (j + 1)]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов. l,m = 1,...,n-1

//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

//#pragma omp  parallel for private(i, j, a)

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
		A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
		A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		A[a][9] = (1. / 12. - 1. / 6.) * Mu(gamma, e_k[i * C_M + (j - 1)]) / (C_hx * C_hy * C_Re);
		A[a][11] = (1. / 4. - 1. / 8.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / (C_hx * C_hy * C_Re);
	}

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;

	a = C_M2 + i * C_M + j;

	A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
	A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
	A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
	A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);


//#pragma omp  parallel for private(i, j, a)

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
		A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
		A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
		A[a][10] = (1. / 6. - 1. / 12.) * Mu(gamma, e_k[i * C_M + (j + 1)]) / (C_hx * C_hy * C_Re);
		A[a][11] = (1. / 8. - 1. / 4.) * Mu(gamma, e_k[(i + 1) * C_M + j]) / (C_hx * C_hy * C_Re);
	}

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;

	a = C_M2 + i * C_M + j;

	A[a][0] = -(Mu(gamma, e_k[(i - 1) * C_M + j]) + Mu(gamma, e_k[i * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][1] = -2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + Mu(gamma, e_k[i * C_M + j])) / (3 * C_hy * C_hy * C_Re);
	A[a][2] = sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau + (Mu(gamma, e_k[(i - 1) * C_M + j]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re) + 2 * (Mu(gamma, e_k[i * C_M + (j - 1)]) + 2 * Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][3] = -2 * (Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[i * C_M + (j + 1)])) / (3 * C_hy * C_hy * C_Re);
	A[a][4] = -(Mu(gamma, e_k[i * C_M + j]) + Mu(gamma, e_k[(i + 1) * C_M + j])) / (2 * C_hx * C_hx * C_Re);
	A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
	A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
	A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
	A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);


//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

//#pragma omp  parallel for private(i, j, a)

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
			A[a][5] = (Mu(gamma, e_k[i * C_M + j - 1]) / 6 - Mu(gamma, e_k[(i - 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
			A[a][6] = (Mu(gamma, e_k[(i - 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j + 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][7] = (Mu(gamma, e_k[(i + 1) * C_M + j]) / 4 - Mu(gamma, e_k[i * C_M + j - 1]) / 6) / (C_hx * C_hy * C_Re);
			A[a][8] = (Mu(gamma, e_k[i * C_M + j + 1]) / 6 - Mu(gamma, e_k[(i + 1) * C_M + j]) / 4) / (C_hx * C_hy * C_Re);
		}
	}

	return 0;
}


//Вектор правых частей системы уравнений
inline double motion_f(double gamma, double* sigma_k, double* sigma_k1, double* u_k, double* v_k, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u

//#pragma omp  parallel for private(i, j, a)

	//Для внутренних узлов. l,m = 1,...,n-1
	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			f[a] = uX_k[a] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			f[a] = uX_k[a] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = i * C_M + j;

			f[a] = uX_k[a] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			f[a] = uX_k[a] * sigma_k1[i * C_M + j] * sigma_k1[i * C_M + j] / C_tau - (P(gamma, sigma_k[(i + 1) * C_M + j], e_k[(i + 1) * C_M + j]) - P(gamma, sigma_k[(i - 1) * C_M + j], e_k[(i - 1) * C_M + j])) / (2 * C_hx);
		}
	}

	return 0;
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
inline double motion_d()
{
	int i = 0;
	int j = 0;
	int a;

//#pragma omp  parallel for private(i, j, a)

	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = C_M2 + i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;
			D[a] = 1 / A[a][2];
		}
	}

	return 0;
}

//Вектор B = A*Xk1
inline double motion_b(double* u_k1, double* v_k1)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u

	//Для внутренних узлов

//#pragma omp  parallel for private(i, j, a)

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


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		j = C_cntr + i + 1 - C_qq;

		a = i * C_M + j;

		B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
			A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
			A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
			A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
			//A[a][7]*v_k1[(i+1)*M+(j-1)] +
			A[a][8] * v_k1[(i + 1) * C_M + (j + 1)] +
			A[a][9] * v_k1[i * C_M + (j - 1)] +
			A[a][11] * v_k1[(i + 1) * C_M + j];
	}

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;

	a = i * C_M + j;

	B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
		A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
		A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
		A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
		A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
		A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];


//#pragma omp  parallel for private(i, j, a)

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

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;

	a = i * C_M + j;

	B[a] = A[a][0] * u_k1[(i - 1) * C_M + j] + A[a][1] * u_k1[i * C_M + (j - 1)] + A[a][2] * u_k1[i * C_M + j] +
		A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
		A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
		A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
		A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
		A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];


//#pragma omp  parallel for private(i, j, a)

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

//#pragma omp  parallel for private(i, j, a)

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

//#pragma omp  parallel for private(i, j, a)

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

	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов

//#pragma omp  parallel for private(i, j, a)

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

//#pragma omp  parallel for private(i, j, a)

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

	i = C_qq + C_w - 1;
	j = C_cntr + i + 1 - C_qq;
	a = C_M2 + i * C_M + j;

	B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
		A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
		A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
		A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
		A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
		A[a][8] * u_k1[(i + 1) * C_M + j + 1];


//#pragma omp  parallel for private(i, j, a)

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

	i = C_qq + C_w - 1;
	j = C_cntr - i - 1 + C_qq;

	a = C_M2 + i * C_M + j;

	B[a] = A[a][0] * v_k1[(i - 1) * C_M + j] + A[a][1] * v_k1[i * C_M + j - 1] + A[a][2] * v_k1[i * C_M + j] +
		A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
		A[a][5] * u_k1[(i - 1) * C_M + j - 1] +
		A[a][6] * u_k1[(i - 1) * C_M + j + 1] +
		A[a][7] * u_k1[(i + 1) * C_M + j - 1] +
		A[a][8] * u_k1[(i + 1) * C_M + j + 1];


//#pragma omp  parallel for private(i, j, a)

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

//#pragma omp  parallel for private(i, j, a)

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

//#pragma omp  parallel for private(i, j, a)

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

	return 0;
}

//Метод Якоби
inline double motion_jakobi(double* u2, double* u_k1, double* v2, double* v_k1)
{
	int i = 0;
	int j = 0;
	int a;

//#pragma omp  parallel for private(i, j, a)

	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			u2[i * C_M + j] = u_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			u2[i * C_M + j] = u_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = i * C_M + j;

			u2[i * C_M + j] = u_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			u2[i * C_M + j] = u_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;

			v2[i * C_M + j] = v_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;

			v2[i * C_M + j] = v_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w; i++)
	{
		for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
		{
			a = C_M2 + i * C_M + j;

			v2[i * C_M + j] = v_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;

			v2[i * C_M + j] = v_k1[i * C_M + j] - D[a] * (B[a] - f[a]);
		}
	}

	return 0;
}

//Метод Гаусса-Зейделя (последовательных смещений)
inline double motion_Zeidel(double* u_k1, double* v_k1, double* u2, double* v2)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u

	//Для внутренних узлов. l,m = 1,...,n-1

//#pragma omp  parallel for private(i, j, a)

	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] /*+ A[a][2]*u_k1[i*C_M+j]*/ + A[a][3] * u_k1[i * C_M + j + 1] +
					A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + j] +
					A[a][6] * v_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * v_k1[(i + 1) * C_M + j] +
					A[a][8] * v_k1[(i + 1) * C_M + j + 1];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < C_N)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)] /*+ A[a][2]*u_k1[i*C_M+j]*/ +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_N)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + j - 1] /*+ A[a][2]*u_k1[i*C_M+j]*/ +
					A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * v_k1[(i - 1) * C_M + j] +
					A[a][7] * v_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * v_k1[(i + 1) * C_M + j];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		for (j = 1; j <= C_cntr - i - 1 + C_qq; j++)
		{
			a = i * C_M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][3] * u_k1[i * C_M + j + 1] +
					A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + j] +
					A[a][6] * v_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * v_k1[(i + 1) * C_M + j] +
					A[a][8] * v_k1[(i + 1) * C_M + j + 1];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < C_cntr - i - 1 + C_qq)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)]  +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_cntr - i - 1 + C_qq)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)]  +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][10] * v_k1[i * C_M + (j + 1)] +
					A[a][11] * v_k1[(i + 1) * C_M + j];
				//A[a][8]*v_k1[(i+1)*C_M+(j+1)];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = C_qq + C_w - 1;

//#pragma omp  parallel for private(j, a)

	for (j = 1; j <= C_cntr - i - 1 + C_qq; j++)
	{
		a = i * C_M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][3] * u_k1[i * C_M + j + 1] +
				A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + j] +
				A[a][6] * v_k1[(i - 1) * C_M + j + 1] +
				A[a][7] * v_k1[(i + 1) * C_M + j] +
				A[a][8] * v_k1[(i + 1) * C_M + j + 1];

			u2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > 0 && j < C_cntr - i - 1 + C_qq)
		{
			B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)]+
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

			u2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == C_cntr - i - 1 + C_qq)
		{
			B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)] +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

			u2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			if (j == C_cntr + i + 1 - C_qq)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)]+
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)] +
					A[a][9] * v_k1[i * C_M + (j - 1)] +
					A[a][11] * v_k1[(i + 1) * C_M + j];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > C_cntr + i + 1 - C_qq && j < C_N)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_N)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + j - 1]  +
					A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * v_k1[(i - 1) * C_M + j] +
					A[a][7] * v_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * v_k1[(i + 1) * C_M + j];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = C_qq + C_w - 1;

//#pragma omp  parallel for private(j, a)

	for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
	{
		a = i * C_M + j;

		if (j == C_cntr + i + 1 - C_qq)
		{
			B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)] /*+ A[a][2]*u_k1[i*C_M+j]*/ +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

			u2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > C_cntr + i + 1 - C_qq && j < C_N)
		{
			B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)] /*+ A[a][2]*u_k1[i*C_M+j]*/ +
				A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

			u2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == C_N)
		{
			B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + j - 1] /*+ A[a][2]*u_k1[i*C_M+j]*/ +
				A[a][4] * u_k1[(i + 1) * C_M + j] +
				A[a][5] * v_k1[(i - 1) * C_M + j - 1] +
				A[a][6] * v_k1[(i - 1) * C_M + j] +
				A[a][7] * v_k1[(i + 1) * C_M + j - 1] +
				A[a][8] * v_k1[(i + 1) * C_M + j];

			u2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = i * C_M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][3] * u_k1[i * C_M + j + 1] +
					A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + j] +
					A[a][6] * v_k1[(i - 1) * C_M + j + 1] +
					A[a][7] * v_k1[(i + 1) * C_M + j] +
					A[a][8] * v_k1[(i + 1) * C_M + j + 1];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < C_N)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + (j - 1)] +
					A[a][3] * u_k1[i * C_M + (j + 1)] + A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * C_M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * C_M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * C_M + (j + 1)];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_N)
			{
				B[a] = A[a][0] * u2[(i - 1) * C_M + j] + A[a][1] * u2[i * C_M + j - 1]  +
					A[a][4] * u_k1[(i + 1) * C_M + j] +
					A[a][5] * v_k1[(i - 1) * C_M + j - 1] +
					A[a][6] * v_k1[(i - 1) * C_M + j] +
					A[a][7] * v_k1[(i + 1) * C_M + j - 1] +
					A[a][8] * v_k1[(i + 1) * C_M + j];

				u2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов. l,m = 1,...,n-1

//#pragma omp  parallel for private(i, j, a)

	for (i = 1; i < C_qq; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j]+ A[a][3] * v_k1[i * C_M + j + 1] +
					A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < C_N)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}


			if (j == C_N)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] /*+ A[a][2]*v_k1[i*C_M+j]*/ +
					A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{		
		for (j = 1; j <= C_cntr - i - 1 + C_qq; j++)
		{
			a = C_M2 + i * C_M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j]+ A[a][3] * v_k1[i * C_M + j + 1] +
					A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < C_cntr - i - 1 + C_qq)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_cntr - i - 1 + C_qq)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][10] * u2[i * C_M + j + 1] +
					A[a][11] * u2[(i + 1) * C_M + j];
				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = C_qq + C_w - 1;

//#pragma omp  parallel for private(j, a)

	for (j = 1; j <= C_cntr - i - 1 + C_qq; j++)
	{
		a = C_M2 + i * C_M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * v2[(i - 1) * C_M + j]  + A[a][3] * v_k1[i * C_M + j + 1] +
				A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u2[(i - 1) * C_M + j] +
				A[a][6] * u2[(i - 1) * C_M + j + 1] +
				A[a][7] * u2[(i + 1) * C_M + j] +
				A[a][8] * u2[(i + 1) * C_M + j + 1];

			v2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > 0 && j < C_cntr - i - 1 + C_qq)
		{
			B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u2[(i - 1) * C_M + j - 1] +
				A[a][6] * u2[(i - 1) * C_M + j + 1] +
				A[a][7] * u2[(i + 1) * C_M + j - 1] +
				A[a][8] * u2[(i + 1) * C_M + j + 1];

			v2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == C_cntr - i - 1 + C_qq)
		{
			B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1]  +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u2[(i - 1) * C_M + j - 1] +
				A[a][6] * u2[(i - 1) * C_M + j + 1] +
				A[a][7] * u2[(i + 1) * C_M + j - 1] +
				A[a][8] * u2[(i + 1) * C_M + j + 1];

			v2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq; i < C_qq + C_w - 1; i++)
	{
		for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;

			if (j == C_cntr + i + 1 - C_qq)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][8] * u2[(i + 1) * C_M + j + 1] +
					A[a][9] * u2[i * C_M + j - 1] +
					A[a][11] * u2[(i + 1) * C_M + j];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > C_cntr + i + 1 - C_qq && j < C_N)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_N)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
					A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = C_qq + C_w - 1;

//#pragma omp  parallel for private(j, a)

	for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
	{
		a = C_M2 + i * C_M + j;

		if (j == C_cntr + i + 1 - C_qq)
		{
			B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u2[(i - 1) * C_M + j - 1] +
				A[a][6] * u2[(i - 1) * C_M + j + 1] +
				A[a][7] * u2[(i + 1) * C_M + j - 1] +
				A[a][8] * u2[(i + 1) * C_M + j + 1];

			v2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > C_cntr + i + 1 - C_qq && j < C_N)
		{
			B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
				A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u2[(i - 1) * C_M + j - 1] +
				A[a][6] * u2[(i - 1) * C_M + j + 1] +
				A[a][7] * u2[(i + 1) * C_M + j - 1] +
				A[a][8] * u2[(i + 1) * C_M + j + 1];

			v2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == C_N)
		{
			B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1]  +
				A[a][4] * v_k1[(i + 1) * C_M + j] +
				A[a][5] * u2[(i - 1) * C_M + j - 1] +
				A[a][6] * u2[(i - 1) * C_M + j] +
				A[a][7] * u2[(i + 1) * C_M + j - 1] +
				A[a][8] * u2[(i + 1) * C_M + j];

			v2[i * C_M + j] = D[a] * (f[a] - B[a]);
		}
	}


//#pragma omp  parallel for private(i, j, a)

	for (i = C_qq + C_w; i < C_M1 - 1; i++)
	{
		for (j = 1; j < C_M - 1; j++)
		{
			a = C_M2 + i * C_M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j]+ A[a][3] * v_k1[i * C_M + j + 1] +
					A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < C_N)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1]  +
					A[a][3] * v_k1[i * C_M + j + 1] + A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j + 1] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j + 1];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == C_N)
			{
				B[a] = A[a][0] * v2[(i - 1) * C_M + j] + A[a][1] * v2[i * C_M + j - 1] +
					A[a][4] * v_k1[(i + 1) * C_M + j] +
					A[a][5] * u2[(i - 1) * C_M + j - 1] +
					A[a][6] * u2[(i - 1) * C_M + j] +
					A[a][7] * u2[(i + 1) * C_M + j - 1] +
					A[a][8] * u2[(i + 1) * C_M + j];

				v2[i * C_M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	return 0;
}


inline int motion(const double gamma, double* sigma_k1, double* sigma_k, double* u_k, double* v_k, double* u_k1, double* v_k1, double* u2, double* v2, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;
	int bl = 1;
	int c_u;
	int c_v;
	int s_m = 0;

	/*---------------------------------------------*/

	motion_a(gamma, sigma_k1, e_k);
	motion_d();
	motion_f(gamma, sigma_k, sigma_k1, u_k, v_k, e_k);
		
	while (bl)
	{
		motion_b(u_k1, v_k1);
		motion_jakobi(u2, u_k1, v2, v_k1);

		c_u = 0;
		c_v = 0;

//#pragma omp  parallel for private(i, j, a)

		for (i = 1; i < C_qq; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				if (fabs(u_k1[a] - u2[a]) <= C_epsilon)
				{
					c_u += 1;
				}

				if (fabs(v_k1[a] - v2[a]) <= C_epsilon)
				{
					c_v += 1;
				}
			}
		}

//#pragma omp  parallel for private(i, j, a)

		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				if (fabs(u_k1[a] - u2[a]) <= C_epsilon)
				{
					c_u += 1;
				}

				if (fabs(v_k1[a] - v2[a]) <= C_epsilon)
				{
					c_v += 1;
				}
			}
		}

//#pragma omp  parallel for private(i, j, a)

		for (i = C_qq; i < C_qq + C_w; i++)
		{
			for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
			{
				a = i * C_M + j;
				if (fabs(u_k1[a] - u2[a]) <= C_epsilon)
				{
					c_u += 1;
				}

				if (fabs(v_k1[a] - v2[a]) <= C_epsilon)
				{
					c_v += 1;
				}
			}
		}

//#pragma omp  parallel for private(i, j, a)

		for (i = C_qq + C_w; i < C_M1 - 1; i++)
		{
			for (j = 1; j < C_M - 1; j++)
			{
				a = i * C_M + j;
				if (fabs(u_k1[a] - u2[a]) <= C_epsilon)
				{
					c_u += 1;
				}

				if (fabs(v_k1[a] - v2[a]) <= C_epsilon)
				{
					c_v += 1;
				}
			}
		}


		if (c_u == (C_N1 - 1) * (C_N - 1) - (2 + (C_q - 1) * 2) / 2 * C_q && c_v >= (C_N1 - 1) * (C_N - 1) - (2 + (C_q - 1) * 2) / 2 * C_q)
		{
			bl = 0;
		}
		else if (s_m > 20)
		{
		
			bl = 0;
		}

		else
		{

//#pragma omp  parallel for private(i, j, a)

			for (i = 1; i < C_qq; i++)
			{
				for (j = 1; j < C_M - 1; j++)
				{
					a = i * C_M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}

//#pragma omp  parallel for private(i, j, a)

			for (i = C_qq; i < C_qq + C_w; i++)
			{
				for (j = C_cntr + i + 1 - C_qq; j < C_M - 1; j++)
				{
					a = i * C_M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}

//#pragma omp  parallel for private(i, j, a)

			for (i = C_qq; i < C_qq + C_w; i++)
			{
				for (j = C_cntr - i - 1 + C_qq; j > 0; j--)
				{
					a = i * C_M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}

//#pragma omp  parallel for private(i, j, a)

			for (i = C_qq + C_w; i < C_M1 - 1; i++)
			{
				for (j = 1; j < C_M - 1; j++)
				{
					a = i * C_M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}
		}
		s_m++;
	}
	return s_m;
}
