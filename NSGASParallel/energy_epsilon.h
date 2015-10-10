/*----- Функция заполняет элементы матрицы для уравнения энергии----*/

inline double energy_a(double gamma, double* sigma_k1, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;
			A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
			A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				gamma / (2 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) +
				gamma / (2 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j - 1]));
			A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
			A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));
		}
	}

	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr + i + 1 - qq;

		a = i * M + j;

		A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
		A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1])) - gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])) - gamma / (8 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j - 1]) + 2 * Mu(gamma, e_k[i * M + j]));
		A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
			gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
			gamma / (8 * hx * hx * PrRe) * (8 * Mu(gamma, e_k[i * M + j]) + 3 * Mu(gamma, e_k[(i + 1) * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
			gamma / (8 * hy * hy * PrRe) * (8 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[i * M + j + 1]) + 3 * Mu(gamma, e_k[i * M + j - 1]));
		A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
		A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j])) - gamma / (4 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])) - gamma / (8 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]));
	}

	i = qq + w - 1;
	j = cntr + i + 1 - qq;

	a = i * M + j;

	A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
	A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
		gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
		gamma / (2 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) +
		gamma / (2 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j - 1]));
	A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
	A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr - i - 1 + qq;

		a = i * M + j;

		A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
		A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
		A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
			gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
			gamma / (8 * hx * hx * PrRe) * (8 * Mu(gamma, e_k[i * M + j]) + 3 * Mu(gamma, e_k[(i + 1) * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
			gamma / (8 * hy * hy * PrRe) * (8 * Mu(gamma, e_k[i * M + j]) + 3 * Mu(gamma, e_k[i * M + j + 1]) + 4 * Mu(gamma, e_k[i * M + j - 1]));
		A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j])) - gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])) - gamma / (8 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]));
		A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j])) - gamma / (4 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])) - gamma / (8 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]));
	}

	i = qq + w - 1;
	j = cntr - i - 1 + qq;

	a = i * M + j;

	A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
	A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
		gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
		gamma / (2 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) +
		gamma / (2 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j - 1]));
	A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
	A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 2 - qq; j < M - 1; j++)
		{
			a = i * M + j;

			A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
			A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				gamma / (2 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) +
				gamma / (2 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j - 1]));
			A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
			A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));
		}
	}


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 2 + qq; j > 0; j--)
		{
			a = i * M + j;

			A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
			A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				gamma / (2 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) +
				gamma / (2 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j - 1]));
			A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
			A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));
		}
	}


	for (i = qq + w; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
			A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau - gamma / (2 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				gamma / (2 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				gamma / (2 * hx * hx * PrRe) * (2 * Mu(gamma,e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) +
				gamma / (2 * hy * hy * PrRe) * (2 * Mu(gamma,e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j - 1]));
			A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
			A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));
		}
	}


	//Для Г5. l = q-1; m = 1,...,q-1;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < cntr + q - 1; j++)
	{
		a = i * M + j;

		A[a][1] = gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
		A[a][2] = sigma_k1[a] * sigma_k1[a] / (2 * tau) - gamma / (4 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - 2 * e_k[(i + 1) * M + j]) -
			gamma / (4 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
			gamma / (4 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + 2 * Mu(gamma, e_k[(i + 1) * M + j])) +
			gamma / (4 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j + 1]));
		A[a][3] = -gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
		A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[(i + 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
	}													  

	//Для Г6. l = 1,...,q-1; m = q-1;
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;

		a = i * M + j;

		A[a][0] = gamma / (4 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
			- gamma / (4 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) - gamma / (8 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) + 2 * Mu(gamma, e_k[(i - 1) * M + j]));
		A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[(i - 1) * M + j]) -
			gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[i * M + j + 1]) +
			gamma / (8 * hx * hx * PrRe) * (4 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
			gamma / (8 * hy * hy * PrRe) * (4 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[i * M + j + 1]));
		A[a][3] = -gamma / (4 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) - gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j])
			- gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j])) - gamma / (8 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j + 1]) + Mu(gamma, e_k[i * M + j]));
		A[a][5] = -gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
		A[a][8] = gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;

		a = i * M + j;

		A[a][0] = gamma / (4 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
			- gamma / (4 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i - 1) * M + j])) - gamma / (8 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) + 2 * Mu(gamma, e_k[(i - 1) * M + j]));
		A[a][1] = gamma / (4 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) + gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1])
			- gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j - 1])) - gamma / (8 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) + 2 * Mu(gamma, e_k[i * M + j - 1]));
		A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[(i - 1) * M + j]) -
			gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[i * M + j - 1]) +
			gamma / (8 * hx * hx * PrRe) * (4 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
			gamma / (8 * hy * hy * PrRe) * (4 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[i * M + j - 1]));
		A[a][6] = -gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
		A[a][7] = gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
	}

	//Для S_qq,N/2+q.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;

	A[a][0] = gamma / (4 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
		- gamma / (4 * hx * hx * PrRe) * (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])) - gamma / (8 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j]));
	A[a][1] = gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (7 * e_k[i * M + j] - 4 * e_k[(i + 1) * M + j] - 3 * e_k[(i - 1) * M + j]) -
		gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (7 * e_k[i * M + j] - 4 * e_k[i * M + j + 1] - 2 * e_k[i * M + j - 1]) +
		gamma / (8 * hx * hx * PrRe) * (7 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[(i + 1) * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
		gamma / (8 * hy * hy * PrRe) * (7 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[i * M + j + 1]) + 2 * Mu(gamma, e_k[i * M + j - 1])) + gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
	A[a][3] = -gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
	A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));
	A[a][5] = -gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);


	//Для S_qq,N/2-q.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;

	A[a][0] = gamma / (4 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
		- gamma / (4 * hx * hx * PrRe) * (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])) - gamma / (8 * hx * hx * PrRe) * (2 * Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j]));
	A[a][1] = gamma / (2 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (7 * e_k[i * M + j] - 4 * e_k[(i + 1) * M + j] - 3 * e_k[(i - 1) * M + j]) -
		gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (7 * e_k[i * M + j] - 2 * e_k[i * M + j + 1] - 4 * e_k[i * M + j - 1]) +
		gamma / (8 * hx * hx * PrRe) * (7 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[(i + 1) * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
		gamma / (8 * hy * hy * PrRe) * (7 * Mu(gamma, e_k[i * M + j]) + 2 * Mu(gamma, e_k[i * M + j + 1]) + 4 * Mu(gamma, e_k[i * M + j - 1])) + gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
	A[a][3] = -gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])));
	A[a][4] = -gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[(i + 1) * M + j])));
	A[a][6] = -gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);


	//Для S_qq,N/2
	i = qq;
	j = cntr;

	a = i * M + j;

	A[a][0] = gamma / (2 * hx * hx * PrRe) * (Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (Mu(gamma, e_k[(i - 1) * M + j]) + Mu(gamma, e_k[i * M + j])));
	A[a][1] = gamma / (4 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) + gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1])
		- gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j])) - gamma / (8 * hy * hy * PrRe) * (2 * Mu(gamma, e_k[i * M + j - 1]) + Mu(gamma, e_k[i * M + j]));
	A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * tau) + 1 / (2 * tau)) - gamma / (8 * hx * hx * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (6 * e_k[i * M + j] - 4 * e_k[(i - 1) * M + j]) -
		gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (6 * e_k[i * M + j] - 3 * e_k[i * M + j + 1] - 3 * e_k[i * M + j - 1]) +
		gamma / (8 * hx * hx * PrRe) * (6 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[(i - 1) * M + j])) +
		gamma / (8 * hy * hy * PrRe) * (6 * Mu(gamma, e_k[i * M + j]) + 4 * Mu(gamma, e_k[i * M + j + 1]) + 4 * Mu(gamma, e_k[i * M + j - 1])) - gamma / (tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
	A[a][3] = -gamma / (4 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) - gamma / (8 * hy * hy * PrRe) * Mu(gamma, e_k[i * M + j]) / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j])
		- gamma / (4 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) + Mu(gamma, e_k[i * M + j + 1])) - gamma / (8 * hy * hy * PrRe) * (Mu(gamma, e_k[i * M + j]) + 2 * Mu(gamma, e_k[i * M + j + 1]));
	A[a][7] = gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);
	A[a][8] = gamma / (2 * tg * hx * hy * PrRe) * Mu(gamma, e_k[i * M + j]);


	return 0;
}

//Вектор правых частей системы линейных уравнений
inline double energy_f(double gamma, double* sigma_k, double* sigma_k1, double* u_k, double* v_k, double* u_k1, double* v_k1, double* e_k)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / tau - P(gamma, sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				Mu(gamma, e_k[i * M + j]) / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				Mu(gamma, e_k[i * M + j]) / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				Mu(gamma, e_k[i * M + j]) / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

					Mu(gamma, e_k[i * M + j]) / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < M - 1; j++)
		{
			a = i * M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / tau - P(gamma, sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				Mu(gamma, e_k[i * M + j]) / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				Mu(gamma, e_k[i * M + j]) / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				Mu(gamma, e_k[i * M + j]) / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

					Mu(gamma, e_k[i * M + j]) / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / tau - P(gamma, sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				Mu(gamma, e_k[i * M + j]) / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				Mu(gamma, e_k[i * M + j]) / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				Mu(gamma, e_k[i * M + j]) / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

					Mu(gamma, e_k[i * M + j]) / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}

	for (i = qq + w; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / tau - P(gamma, sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				Mu(gamma, e_k[i * M + j]) / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				Mu(gamma, e_k[i * M + j]) / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				Mu(gamma, e_k[i * M + j]) / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

					Mu(gamma, e_k[i * M + j]) / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}

	//Для Г5. l = q-1; m = 1,...,q-1;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < cntr + q - 1; j++)
	{
		a = i * M + j;

		f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / (2 * tau) - P(gamma, sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * ((2 * u_k1[(i + 1) * M + j] - 2 * u_k1[i * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

			Mu(gamma, e_k[i * M + j]) / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j])) +

			Mu(gamma, e_k[i * M + j]) / (12 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

			Mu(gamma, e_k[i * M + j]) / (8 * Re * e_k[i * M + j]) * (((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j + 1] - u_k1[i * M + j]) / hy) * ((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j + 1] - u_k1[i * M + j]) / hy) +
				((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j] - u_k1[i * M + j - 1]) / hy) * ((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j] - u_k1[i * M + j - 1]) / hy)) +

				Mu(gamma, e_k[i * M + j]) / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
				(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
	}

	//Для Г6. l = 1,...,q-1; m = q-1;
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;

		a = i * M + j;

		f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - P(gamma, sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j + 1] / hy - v_k1[i * M + j] / hy) -
			-P(gamma, sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (-u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j + 1] / hy) +

			Mu(gamma, e_k[i * M + j]) / (24 * hx * hx * Re * e_k[i * M + j]) * (1 * -u_k1[i * M + j] * -u_k1[i * M + j] + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

			Mu(gamma, e_k[i * M + j]) / (24 * hy * hy * Re * e_k[i * M + j]) * (3 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + 1 * v_k1[i * M + j] * v_k1[i * M + j]) +

			Mu(gamma, e_k[i * M + j]) / (16 * Re * e_k[i * M + j]) * (1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
				2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
				1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy)) +

				Mu(gamma, e_k[i * M + j]) / (24 * Re * e_k[i * M + j]) * (1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
				2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
				1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy));
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;

		a = i * M + j;

		f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - P(gamma, sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy - v_k1[i * M + j - 1] / hy) -
			-P(gamma, sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (-u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j - 1] / hy) +

			Mu(gamma, e_k[i * M + j]) / (24 * hx * hx * Re * e_k[i * M + j]) * (1 * -u_k1[i * M + j] * -u_k1[i * M + j] + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

			Mu(gamma, e_k[i * M + j]) / (24 * hy * hy * Re * e_k[i * M + j]) * (3 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1]) + 1 * -v_k1[i * M + j] * -v_k1[i * M + j]) +

			Mu(gamma, e_k[i * M + j]) / (16 * Re * e_k[i * M + j]) * (1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
				2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
				1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy)) +

				Mu(gamma, e_k[i * M + j]) / (24 * Re * e_k[i * M + j]) * (1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
				2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
				1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy));
	}

	//Для S_qq,N/2+q.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;

	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - P(gamma, sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (2 * u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + 2 * v_k1[i * M + j + 1] / hy - v_k1[i * M + j] / hy - v_k1[i * M + j - 1] / hy) -
		P(gamma, sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) +

		Mu(gamma, e_k[i * M + j]) / (24 * hx * hx * Re * e_k[i * M + j]) * (4 * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

		Mu(gamma, e_k[i * M + j]) / (24 * hy * hy * Re * e_k[i * M + j]) * (4 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + 2 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1]) + v_k1[i * M + j] * v_k1[i * M + j]) +

		Mu(gamma, e_k[i * M + j]) / (16 * Re * e_k[i * M + j]) * (2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
			1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy)) +

			Mu(gamma, e_k[i * M + j]) / (24 * Re * e_k[i * M + j]) * (2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
			1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy));


	//Для S_qq,N/2-q.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;

	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - P(gamma, sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (2 * u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx + v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy - u_k1[(i - 1) * M + j] / hx - 2 * v_k1[i * M + j - 1] / hy) -
		P(gamma, sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy) +

		Mu(gamma, e_k[i * M + j]) / (24 * hx * hx * Re * e_k[i * M + j]) * (4 * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

		Mu(gamma, e_k[i * M + j]) / (24 * hy * hy * Re * e_k[i * M + j]) * (4 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1]) + 2 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + -v_k1[i * M + j] * -v_k1[i * M + j]) +

		Mu(gamma, e_k[i * M + j]) / (16 * Re * e_k[i * M + j]) * (2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy) +
			2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
			2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

			Mu(gamma, e_k[i * M + j]) / (24 * Re * e_k[i * M + j]) * (2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) +
			2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
			2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));


	//Для S_qq,N/2
	i = qq;
	j = cntr;

	a = i * M + j;

	f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * tau) + 1 / (2 * tau)) - P(gamma, sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (2 * u_k1[i * M + j] / hx - 2 * u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j + 1] / hy - v_k1[i * M + j - 1] / hy) -
		P(gamma, sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (-2 * u_k1[i * M + j] / hx + v_k1[i * M + j + 1] / hy - v_k1[i * M + j - 1] / hy) +

		Mu(gamma, e_k[i * M + j]) / (12 * hx * hx * Re * e_k[i * M + j]) * (2 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) + -u_k1[i * M + j] * -u_k1[i * M + j]) +

		Mu(gamma, e_k[i * M + j]) / (24 * hy * hy * Re * e_k[i * M + j]) * (3 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + 3 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

		Mu(gamma, e_k[i * M + j]) / (16 * Re * e_k[i * M + j]) * (2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
			2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

			Mu(gamma, e_k[i * M + j]) / (24 * Re * e_k[i * M + j]) * (2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
			2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));

	return 0;
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
inline double energy_d()
{
	int i = 0, j = 0, a;

	for (i = 1; i < qq + 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = cntr + i - qq; j < M - 1; j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = cntr - i + qq; j > 0; j--)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq + w - 1; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	return 0;
}

//Вектор B = A*Xk1
inline double energy_b(double* e_k1)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < M - 1; j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;

			B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];
		}
	}

	for (i = qq + w; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];
		}
	}

	//Для Г5. l = q-1; m = 1,...,q-1;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < cntr + q - 1; j++)
	{
		a = i * M + j;

		B[a] = A[a][1] * e_k1[i * M + j - 1] + A[a][2] * e_k1[i * M + j] + A[a][3] * e_k1[i * M + j + 1] +
			A[a][4] * e_k1[(i + 1) * M + j];
	}

	//Для Г6. l = 1,...,q-1; m = q-1;
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;

		a = i * M + j;

		B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][2] * e_k1[i * M + j] + A[a][3] * e_k1[i * M + j + 1]
			+ A[a][5] * e_k1[(i - 1) * M + j - 1] + A[a][8] * e_k1[(i + 1) * M + j + 1];
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;

		a = i * M + j;

		B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + j - 1] + A[a][2] * e_k1[i * M + j]
			+ A[a][6] * e_k1[(i - 1) * M + j + 1] + A[a][7] * e_k1[(i + 1) * M + j - 1];
	}

	//Для S_qq,N/2+q.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;

	B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
		A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j]
		+ A[a][5] * e_k1[(i - 1) * M + j - 1];

	//Для S_qq,N/2-q.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;

	B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
		A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j]
		+ A[a][6] * e_k1[(i - 1) * M + j + 1];


	//Для S_qq,N/2
	i = qq;
	j = cntr;

	a = i * M + j;

	B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
		A[a][3] * e_k1[i * M + (j + 1)]
		+ A[a][7] * e_k1[(i + 1) * M + j - 1] + A[a][8] * e_k1[(i + 1) * M + j + 1];

	return 0;
}

//Метод Якоби
inline double energy_jakobi(double* e_k1, double* e2)
{
	int i = 0, j = 0, a;

	for (i = 1; i < qq + 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			e2[i * M + j] = e_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = cntr + i - qq; j < M - 1; j++)
		{
			a = i * M + j;

			e2[i * M + j] = e_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = cntr - i + qq; j > 0; j--)
		{
			a = i * M + j;

			e2[i * M + j] = e_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq + w - 1; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			e2[i * M + j] = e_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	return 0;
}

//Метод Гаусса-Зейделя
inline double energy_zeidel(double* e_k1, double* e2)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < M - 1)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] /*+ A[a][2]*e_k1[i*M+j]*/ +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	i = qq;
	for (j = 1; j < M - 1; j++)
	{
		a = i * M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][3] * e_k1[i * M + j + 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > 0 && j < M - 1)
		{
			if (j == cntr)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
					A[a][3] * e_k1[i * M + (j + 1)]
					+ A[a][7] * e_k1[(i + 1) * M + j - 1] + A[a][8] * e_k1[(i + 1) * M + j + 1];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			else
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}

		if (j == N)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = 1; j <= cntr - i + qq; j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][3] * e_k1[i * M + j + 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < cntr - i + qq)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)]+
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == cntr - i + qq)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] + A[a][6] * e2[(i - 1) * M + j + 1] + A[a][7] * e_k1[(i + 1) * M + j - 1];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = cntr + i - qq; j < M - 1; j++)
		{
			a = i * M + j;

			if (j == cntr + i - qq)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][3] * e_k1[i * M + j + 1]
					+ A[a][5] * e2[(i - 1) * M + j - 1] + A[a][8] * e_k1[(i + 1) * M + j + 1];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > cntr + i - qq && j < M - 1)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = qq + w - 1;
	for (j = 1; j < M - 1; j++)
	{
		a = i * M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][3] * e_k1[i * M + j + 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > 0 && j < cntr - w + 1)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == cntr - w + 1)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j]
				+ A[a][6] * e2[(i - 1) * M + j + 1];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j >= cntr - q + 2 && j < cntr + q - 1)
		{
			B[a] = A[a][1] * e2[i * M + j - 1] + A[a][3] * e_k1[i * M + j + 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == cntr + w - 1)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j]
				+ A[a][5] * e2[(i - 1) * M + j - 1];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j > cntr + w - 1 && j < N)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == N)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq + w; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][3] * e_k1[i * M + j + 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j > 0 && j < M - 1)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] +
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	return 0;
}


inline int energy(const double gamma, double* sigma_k1, double* sigma_k, double* u_k, double* v_k, double* u_k1, double* v_k1, double* e2, double* e_k, double* e_k1)
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
		energy_b(e_k1);
		energy_jakobi(e_k1, e2);

		c = 0;
		for (i = 1; i < qq + 1; i++)
		{
			for (j = 1; j < M - 1; j++)
			{
				a = i * M + j;
				if (fabs(e_k1[a] - e2[a]) <= epsilon)
				{
					c += 1;
				}
			}
		}

		for (i = qq + 1; i < qq + w - 1; i++)
		{
			for (j = cntr + i - qq; j < M - 1; j++)
			{
				a = i * M + j;
				if (fabs(e_k1[a] - e2[a]) <= epsilon)
				{
					c += 1;
				}
			}
		}

		for (i = qq + 1; i < qq + w - 1; i++)
		{
			for (j = cntr - i + qq; j > 0; j--)
			{
				a = i * M + j;
				if (fabs(e_k1[a] - e2[a]) <= epsilon)
				{
					c += 1;
				}
			}
		}

		for (i = qq + w - 1; i < M1 - 1; i++)
		{
			for (j = 1; j < M - 1; j++)
			{
				a = i * M + j;
				if (fabs(e_k1[a] - e2[a]) <= epsilon)
				{
					c += 1;
				}
			}
		}

		if (c == (N1 - 1) * (N - 1) - (2 + (q - 2 - 1) * 2) / 2 * (q - 2))
		{
			bl = 0;
		}
		else if (s_e > 20)
		{
			bl = 0;
		}

		else
		{
			for (i = 1; i < qq + 1; i++)
			{
				for (j = 1; j < M - 1; j++)
				{
					a = i * M + j;
					e_k1[a] = e2[a];
				}
			}

			for (i = qq + 1; i < qq + w - 1; i++)
			{
				for (j = cntr + i - qq; j < M - 1; j++)
				{
					a = i * M + j;
					e_k1[a] = e2[a];
				}
			}

			for (i = qq + 1; i < qq + w - 1; i++)
			{
				for (j = cntr - i + qq; j > 0; j--)
				{
					a = i * M + j;
					e_k1[a] = e2[a];
				}
			}

			for (i = qq + w - 1; i < M1 - 1; i++)
			{
				for (j = 1; j < M - 1; j++)
				{
					a = i * M + j;
					e_k1[a] = e2[a];
				}
			}
		}
		s_e++;
	}

	return s_e;
}
