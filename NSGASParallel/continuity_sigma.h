inline double continuity(double* sigma_k1, double* u_k, double* v_k)
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
			sigma_k1[a] = sigmaX_k[a] / tau / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < M - 1; j++)
		{
			a = i * M + j;
			sigma_k1[a] = sigmaX_k[a] / tau / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;
			sigma_k1[a] = sigmaX_k[a] / tau / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}

	for (i = qq + w; i < M1 - 1; i++)
	{
		for (j = 1; j < M - 1; j++)
		{
			a = i * M + j;
			sigma_k1[a] = sigmaX_k[a] / tau / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}

	//Для Г5. l = w-1; m = 1,...,q-2;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < cntr + q - 1; j++)
	{
		a = i * M + j;
		sigma_k1[a] = sigmaX_k[a] / (2 * tau) / (1 / (2 * tau) + (u_k[(i + 1) * M + j] - u_k[i * M + j]) / (4 * hx)
			+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (8 * hy));
	}

	//Для Г6.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;
		a = i * M + j;
		sigma_k1[a] = sigmaX_k[a] * (1 / (4 * tau) + 1 / (4 * tau)) / (1 / (4 * tau) + 1 / (4 * tau) + (u_k[i * M + j] - u_k[(i - 1) * M + j]) / (8 * hx)
			- u_k[(i - 1) * M + j] / (16 * hx) + (v_k[i * M + j + 1] - v_k[i * M + j]) / (8 * hy) + v_k[i * M + j + 1] / (16 * hy));
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;
		a = i * M + j;
		sigma_k1[a] = sigmaX_k[a] * (1 / (4 * tau) + 1 / (4 * tau)) / (1 / (4 * tau) + 1 / (4 * tau) + (u_k[i * M + j] - u_k[(i - 1) * M + j]) / (8 * hx)
			- u_k[(i - 1) * M + j] / (16 * hx) + (v_k[i * M + j] - v_k[i * M + j - 1]) / (8 * hy) - v_k[i * M + j - 1] / (16 * hy));
	}

	//Для S_w-1q-1.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;
	sigma_k1[a] = sigmaX_k[a] * (3 / (4 * tau) + 1 / (8 * tau)) / (3 / (4 * tau) + 1 / (8 * tau) + (2 * u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j] - u_k[i * M + j]) / (8 * hx)
		+ (u_k[i * M + j] - u_k[(i - 1) * M + j]) / (16 * hx) + (2 * v_k[i * M + j + 1] - v_k[i * M + j - 1] - v_k[i * M + j]) / (8 * hy) + v_k[i * M + j] / (16 * hy));

	//Для S.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;
	sigma_k1[a] = sigmaX_k[a] * (3 / (4 * tau) + 1 / (8 * tau)) / (3 / (4 * tau) + 1 / (8 * tau) + (2 * u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j] - u_k[i * M + j]) / (8 * hx)
		+ (u_k[i * M + j] - u_k[(i - 1) * M + j]) / (16 * hx) + (v_k[i * M + j + 1] - 2 * v_k[i * M + j - 1] + v_k[i * M + j]) / (8 * hy) - v_k[i * M + j] / (16 * hy));

	//Для S_qq_0
	i = qq;
	j = cntr;

	a = i * M + j;
	sigma_k1[a] = sigmaX_k[a] * (1 / (4 * tau) + 1 / (2 * tau)) / (1 / (4 * tau) + 1 / (2 * tau) + (u_k[i * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
		+ (u_k[(i + 1) * M + j] - u_k[i * M + j]) / (8 * hx) + (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (8 * hy) + (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (16 * hy));

	return 0;
}
