/* ----- Функция заполняет элементы матрицы, составленной для двух уравнении движения.---- */
inline void mtn_calculate_common(const double gamma, const int m_i, const int m1_i, const int m2_i, const int qq_i, const int w_i, 
	const int cntr_i, const double tau_d, double* sigma_k1, double* e_k, const double* e_k_mu)
{
	const double c_coef = C_hx * C_hy * C_Re;
	const double c_coef2 = 3 * C_hx * C_hx * C_Re;
	const double c_coef3 = 2 * C_hy * C_hy * C_Re;
	const double c_coef4 = 2 * C_hx * C_hx * C_Re;
	const double c_coef5 = 3 * C_hy * C_hy * C_Re;
	const double c_coef6 = 4 * C_hy * C_hy * C_Re;
	const double c_coef7 = 6 * C_hx * C_hx * C_Re;
	const double c_coef8 = 8 * C_hy * C_hy * C_Re;	

	// Уравнение для u	
#pragma omp parallel 
	{
#pragma omp for collapse(2) nowait

		for (int i = 1; i < m1_i - 1; ++i)
		{
			for (int j = 1; j < m_i - 1; ++j)
			{
				if (i < qq_i || i >= qq_i + w_i)
				{
					// u
					int a = i * m_i + j;
					A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
					A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3;
					A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
					A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
					A[a][4] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2;
					A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
					A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
					A[a][7] = (e_k_mu[i * m_i + j - 1] / 4 - e_k_mu[(i + 1) * m_i + j] / 6) / c_coef;
					A[a][8] = (e_k_mu[(i + 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j + 1)] / 4) / c_coef;
					D[a] = 1 / A[a][2];
					f[a] = uX_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - (P(gamma, sigma_k[(i + 1) * m_i + j], e_k[(i + 1) * m_i + j]) - P(gamma, sigma_k[(i - 1) * m_i + j], e_k[(i - 1) * m_i + j])) / (2 * C_hx);

					// v
					a += m2_i;
					A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
					A[a][1] = -2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5;
					A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
					A[a][3] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
					A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4;
					A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
					A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
					A[a][7] = (e_k_mu[(i + 1) * m_i + j] / 4 - e_k_mu[i * m_i + j - 1] / 6) / c_coef;
					A[a][8] = (e_k_mu[i * m_i + j + 1] / 6 - e_k_mu[(i + 1) * m_i + j] / 4) / c_coef;
					D[a] = 1 / A[a][2];
				}
			}
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i - 1; i++)
		{
			// u
			int j = cntr_i + i + 1 - qq_i;
			int a = i * m_i + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef6 - (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j]) / c_coef8;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3
				+ (e_k_mu[(i + 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2 + (e_k_mu[(i + 1) * m_i + j] + 2 * e_k_mu[i * m_i + j]) / c_coef7
				+ (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j - 1)]) / c_coef6 + (2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j - 1)]) / c_coef8;
			A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
			A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 - (2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef7;
			A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
			A[a][8] = (e_k_mu[(i + 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j + 1)] / 4) / c_coef;
			A[a][9] = (1. / 4. - 1. / 8.) * e_k_mu[i * m_i + (j - 1)] / c_coef;
			A[a][11] = (1. / 12. - 1. / 6.) * e_k_mu[(i + 1) * m_i + j] / c_coef;

			// v
			a += m2_i;
			A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
			A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5 - (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j]) / (6 * C_hy * C_hy * C_Re);
			A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5
				+ (e_k_mu[(i + 1) * m_i + j] + e_k_mu[i * m_i + j]) / (4 * C_hx * C_hx * C_Re) + (e_k_mu[(i + 1) * m_i + j] + 2 * e_k_mu[i * m_i + j]) / (8 * C_hx * C_hx * C_Re)
				+ (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j - 1)]) / c_coef5 + (2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j - 1)]) / (6 * C_hy * C_hy * C_Re);
			A[a][3] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
			A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / (4 * C_hx * C_hx * C_Re) - (2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / (8 * C_hx * C_hx * C_Re);
			A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
			A[a][8] = (e_k_mu[i * m_i + j + 1] / 6 - e_k_mu[(i + 1) * m_i + j] / 4) / c_coef;
			A[a][9] = (1. / 12. - 1. / 6.) * e_k_mu[i * m_i + (j - 1)] / c_coef;
			A[a][11] = (1. / 4. - 1. / 8.) * e_k_mu[(i + 1) * m_i + j] / c_coef;

			// u
			j = cntr_i - i - 1 + qq_i;
			a = i * m_i + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3
				+ (e_k_mu[(i + 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2 + (e_k_mu[(i + 1) * m_i + j] + 2 * e_k_mu[i * m_i + j]) / c_coef7
				+ (e_k_mu[i * m_i + (j + 1)] + e_k_mu[i * m_i + j]) / c_coef6 + (e_k_mu[i * m_i + (j + 1)] + 2 * e_k_mu[i * m_i + j]) / c_coef8;
			A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef6 - (2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef8;
			A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 - (2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef7;
			A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
			A[a][7] = (e_k_mu[i * m_i + j - 1] / 4 - e_k_mu[(i + 1) * m_i + j] / 6) / c_coef;
			A[a][10] = (1. / 8. - 1. / 4.) * e_k_mu[i * m_i + (j + 1)] / c_coef;
			A[a][11] = (1. / 6. - 1. / 12.) * e_k_mu[(i + 1) * m_i + j] / c_coef;

			// v
			a += m2_i;
			A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
			A[a][1] = -2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5;
			A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5
				+ (e_k_mu[(i + 1) * m_i + j] + e_k_mu[i * m_i + j]) / (4 * C_hx * C_hx * C_Re) + (e_k_mu[(i + 1) * m_i + j] + 2 * e_k_mu[i * m_i + j]) / (8 * C_hx * C_hx * C_Re)
				+ (e_k_mu[i * m_i + (j + 1)] + e_k_mu[i * m_i + j]) / c_coef5 + (e_k_mu[i * m_i + (j + 1)] + 2 * e_k_mu[i * m_i + j]) / (6 * C_hy * C_hy * C_Re);
			A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5 - (2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / (6 * C_hy * C_hy * C_Re);
			A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / (4 * C_hx * C_hx * C_Re) - (2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / (8 * C_hx * C_hx * C_Re);
			A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
			A[a][7] = (e_k_mu[(i + 1) * m_i + j] / 4 - e_k_mu[i * m_i + j - 1] / 6) / c_coef;
			A[a][10] = (1. / 6. - 1. / 12.) * e_k_mu[i * m_i + (j + 1)] / c_coef;
			A[a][11] = (1. / 8. - 1. / 4.) * e_k_mu[(i + 1) * m_i + j] / c_coef;
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i; i++)
		{
			for (int j = cntr_i + i + 2 - qq_i; j < m_i - 1; j++)
			{
				// u
				int a = i * m_i + j;
				A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
				A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3;
				A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
				A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
				A[a][4] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2;
				A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
				A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
				A[a][7] = (e_k_mu[i * m_i + j - 1] / 4 - e_k_mu[(i + 1) * m_i + j] / 6) / c_coef;
				A[a][8] = (e_k_mu[(i + 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j + 1)] / 4) / c_coef;

				// v
				a += m2_i;
				A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
				A[a][1] = -2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5;
				A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
				A[a][3] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
				A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4;
				A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
				A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
				A[a][7] = (e_k_mu[(i + 1) * m_i + j] / 4 - e_k_mu[i * m_i + j - 1] / 6) / c_coef;
				A[a][8] = (e_k_mu[i * m_i + j + 1] / 6 - e_k_mu[(i + 1) * m_i + j] / 4) / c_coef;
			}
			for (int j = cntr_i - i - 2 + qq_i; j > 0; j--)
			{
				// u
				int a = i * m_i + j;
				A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
				A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3;
				A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
				A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
				A[a][4] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2;
				A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
				A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
				A[a][7] = (e_k_mu[i * m_i + j - 1] / 4 - e_k_mu[(i + 1) * m_i + j] / 6) / c_coef;
				A[a][8] = (e_k_mu[(i + 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j + 1)] / 4) / c_coef;

				// v
				a += m2_i;
				A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
				A[a][1] = -2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5;
				A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
				A[a][3] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
				A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4;
				A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
				A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
				A[a][7] = (e_k_mu[(i + 1) * m_i + j] / 4 - e_k_mu[i * m_i + j - 1] / 6) / c_coef;
				A[a][8] = (e_k_mu[i * m_i + j + 1] / 6 - e_k_mu[(i + 1) * m_i + j] / 4) / c_coef;
			}
		}
#pragma omp single nowait
		{
			// u
			int i = qq_i + w_i - 1;
			int j = cntr_i + i + 1 - qq_i;
			int a = i * m_i + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
			A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
			A[a][4] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2;
			A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
			A[a][7] = (e_k_mu[i * m_i + j - 1] / 4 - e_k_mu[(i + 1) * m_i + j] / 6) / c_coef;
			A[a][8] = (e_k_mu[(i + 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j + 1)] / 4) / c_coef;
			// v
			a += m2_i;
			A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
			A[a][1] = -2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5;
			A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
			A[a][3] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
			A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4;
			A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
			A[a][7] = (e_k_mu[(i + 1) * m_i + j] / 4 - e_k_mu[i * m_i + j - 1] / 6) / c_coef;
			A[a][8] = (e_k_mu[i * m_i + j + 1] / 6 - e_k_mu[(i + 1) * m_i + j] / 4) / c_coef;

			// u
			j = cntr_i - i - 1 + qq_i;
			a = i * m_i + j;
			A[a][0] = -2 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef2;
			A[a][1] = -(e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef3;
			A[a][2] = sigma_k1[a] * sigma_k1[a] / tau_d + 2 * (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2 + (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
			A[a][3] = -(e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef3;
			A[a][4] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef2;
			A[a][5] = (e_k_mu[(i - 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j - 1)] / 4) / c_coef;
			A[a][6] = (e_k_mu[i * m_i + j + 1] / 4 - e_k_mu[(i - 1) * m_i + j] / 6) / c_coef;
			A[a][7] = (e_k_mu[i * m_i + j - 1] / 4 - e_k_mu[(i + 1) * m_i + j] / 6) / c_coef;
			A[a][8] = (e_k_mu[(i + 1) * m_i + j] / 6 - e_k_mu[i * m_i + (j + 1)] / 4) / c_coef;
			// v
			a += m2_i;
			A[a][0] = -(e_k_mu[(i - 1) * m_i + j] + e_k_mu[i * m_i + j]) / c_coef4;
			A[a][1] = -2 * (e_k_mu[i * m_i + (j - 1)] + e_k_mu[i * m_i + j]) / c_coef5;
			A[a][2] = sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / tau_d + (e_k_mu[(i - 1) * m_i + j] + 2 * e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4 + 2 * (e_k_mu[i * m_i + (j - 1)] + 2 * e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
			A[a][3] = -2 * (e_k_mu[i * m_i + j] + e_k_mu[i * m_i + (j + 1)]) / c_coef5;
			A[a][4] = -(e_k_mu[i * m_i + j] + e_k_mu[(i + 1) * m_i + j]) / c_coef4;
			A[a][5] = (e_k_mu[i * m_i + j - 1] / 6 - e_k_mu[(i - 1) * m_i + j] / 4) / c_coef;
			A[a][6] = (e_k_mu[(i - 1) * m_i + j] / 4 - e_k_mu[i * m_i + j + 1] / 6) / c_coef;
			A[a][7] = (e_k_mu[(i + 1) * m_i + j] / 4 - e_k_mu[i * m_i + j - 1] / 6) / c_coef;
			A[a][8] = (e_k_mu[i * m_i + j + 1] / 6 - e_k_mu[(i + 1) * m_i + j] / 4) / c_coef;
		}
	} // #pragma omp parallel

	// Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
#pragma omp parallel for 
	for (int i = qq_i; i < qq_i + w_i; i++)
	{
		for (int j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
		{
			D[i * m_i + j] = 1 / A[i * m_i + j][2];
			D[m2_i + i * m_i + j] = 1 / A[m2_i + i * m_i + j][2];
			//Вектор правых частей системы уравнений
			f[i * m_i + j] = uX_k[i * m_i + j] * sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / C_tau - (P(gamma, sigma_k[(i + 1) * m_i + j], e_k[(i + 1) * m_i + j]) - P(gamma, sigma_k[(i - 1) * m_i + j], e_k[(i - 1) * m_i + j])) / (2 * C_hx);
		}
		for (int j = cntr_i - i - 1 + qq_i; j > 0; j--)
		{
			D[i * m_i + j] = 1 / A[i * m_i + j][2];
			D[m2_i + i * m_i + j] = 1 / A[m2_i + i * m_i + j][2];
			//Вектор правых частей системы уравнений
			f[i * m_i + j] = uX_k[i * m_i + j] * sigma_k1[i * m_i + j] * sigma_k1[i * m_i + j] / C_tau - (P(gamma, sigma_k[(i + 1) * m_i + j], e_k[(i + 1) * m_i + j]) - P(gamma, sigma_k[(i - 1) * m_i + j], e_k[(i - 1) * m_i + j])) / (2 * C_hx);
		}
	} // #pragma omp parallel
}

//Вектор B = A*Xk1
inline void mtn_calculate_jakobi(const int qq_i, const int m_i, const int w_i, const int cntr_i, const int m1_i, const int m2_i, double* u_k1, double* v_k1)
{
#pragma omp parallel
	{
#pragma omp for collapse(2) nowait
		for (int i = 1; i < m1_i - 1; i++)
		{
			for (int j = 1; j < m_i - 1; j++)
			{
				if (i < qq_i || i >= qq_i + w_i)
				{
					// u
					int a = i * m_i + j;
					B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
						A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
						A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
						A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
						A[a][7] * v_k1[(i + 1) * m_i + (j - 1)] +
						A[a][8] * v_k1[(i + 1) * m_i + (j + 1)];

					u2[a] = u_k1[a] - D[a] * (B[a] - f[a]);

					// v
					a += m2_i;
					B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
						A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
						A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
						A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
						A[a][7] * u_k1[(i + 1) * m_i + j - 1] +
						A[a][8] * u_k1[(i + 1) * m_i + j + 1];

					v2[a - m2_i] = v_k1[a - m2_i] - D[a] * (B[a] - f[a]);
				}
			}
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i - 1; i++)
		{
			// u
			int j = cntr_i + i + 1 - qq_i;
			int a = i * m_i + j;
			B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
				A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
				A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
				A[a][8] * v_k1[(i + 1) * m_i + (j + 1)] +
				A[a][9] * v_k1[i * m_i + (j - 1)] +
				A[a][11] * v_k1[(i + 1) * m_i + j];
			// v
			a += m2_i;
			B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
				A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
				A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
				A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
				A[a][8] * u_k1[(i + 1) * m_i + j + 1] +
				A[a][9] * u_k1[i * m_i + j - 1] +
				A[a][11] * u_k1[(i + 1) * m_i + j];

			// u
			j = cntr_i - i - 1 + qq_i;
			a = i * m_i + j;
			B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
				A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
				A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * m_i + (j - 1)] +
				A[a][10] * v_k1[i * m_i + (j + 1)] +
				A[a][11] * v_k1[(i + 1) * m_i + j];
			// v
			a += m2_i;
			B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
				A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
				A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
				A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
				A[a][7] * u_k1[(i + 1) * m_i + j - 1] +
				A[a][10] * u_k1[i * m_i + j + 1] +
				A[a][11] * u_k1[(i + 1) * m_i + j];
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i; i++)
		{
			for (int j = cntr_i + i + 2 - qq_i; j < m_i - 1; j++)
			{
				// u
				int a = i * m_i + j;
				B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
					A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
					A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * m_i + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * m_i + (j + 1)];
				// v
				a += m2_i;
				B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
					A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
					A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
					A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
					A[a][7] * u_k1[(i + 1) * m_i + j - 1] +
					A[a][8] * u_k1[(i + 1) * m_i + j + 1];
			}
			for (int j = cntr_i - i - 2 + qq_i; j > 0; j--)
			{
				// u
				int a = i * m_i + j;
				B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
					A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
					A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * m_i + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * m_i + (j + 1)];
				// v
				a += m2_i;
				B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
					A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
					A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
					A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
					A[a][7] * u_k1[(i + 1) * m_i + j - 1] +
					A[a][8] * u_k1[(i + 1) * m_i + j + 1];
			}
		}
#pragma omp single nowait
		{
			// u
			int i = qq_i + w_i - 1;
			int j = cntr_i + i + 1 - qq_i;
			int a = i * m_i + j;
			B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
				A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
				A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * m_i + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * m_i + (j + 1)];

			// v
			a += m2_i;
			B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
				A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
				A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
				A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
				A[a][7] * u_k1[(i + 1) * m_i + j - 1] +
				A[a][8] * u_k1[(i + 1) * m_i + j + 1];

			// u
			j = cntr_i - i - 1 + qq_i;
			a = i * m_i + j;
			B[a] = A[a][0] * u_k1[(i - 1) * m_i + j] + A[a][1] * u_k1[i * m_i + (j - 1)] + A[a][2] * u_k1[i * m_i + j] +
				A[a][3] * u_k1[i * m_i + (j + 1)] + A[a][4] * u_k1[(i + 1) * m_i + j] +
				A[a][5] * v_k1[(i - 1) * m_i + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * m_i + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * m_i + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * m_i + (j + 1)];
			// v
			a += m2_i;
			B[a] = A[a][0] * v_k1[(i - 1) * m_i + j] + A[a][1] * v_k1[i * m_i + j - 1] + A[a][2] * v_k1[i * m_i + j] +
				A[a][3] * v_k1[i * m_i + j + 1] + A[a][4] * v_k1[(i + 1) * m_i + j] +
				A[a][5] * u_k1[(i - 1) * m_i + j - 1] +
				A[a][6] * u_k1[(i - 1) * m_i + j + 1] +
				A[a][7] * u_k1[(i + 1) * m_i + j - 1] +
				A[a][8] * u_k1[(i + 1) * m_i + j + 1];
		}
	} // #pragma omp parallel
	// Метод Якоби остаток
#pragma omp parallel
	{
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i; i++)
		{
			for (int j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
			{
				u2[i * m_i + j] = u_k1[i * m_i + j] - D[i * m_i + j] * (B[i * m_i + j] - f[i * m_i + j]);
				v2[i * m_i + j] = v_k1[i * m_i + j] - D[m2_i + i * m_i + j] * (B[m2_i + i * m_i + j] - f[m2_i + i * m_i + j]);
			}
			for (int j = cntr_i - i - 1 + qq_i; j > 0; j--)
			{
				u2[i * m_i + j] = u_k1[i * m_i + j] - D[i * m_i + j] * (B[i * m_i + j] - f[i * m_i + j]);
				v2[i * m_i + j] = v_k1[i * m_i + j] - D[m2_i + i * m_i + j] * (B[m2_i + i * m_i + j] - f[m2_i + i * m_i + j]);
			}
		}
	} // #pragma omp parallel
}

// m_i = C_M
// m1_i = C_M1
// qq_i = C_qq
// w_i = C_w
// cntr_i = C_cntr
inline int motion(const double gamma, const double epsilon_d, const int m_i, const int m1_i, const int m2_i, const int qq_i,
                  const int w_i, const int cntr_i,
				  const int q_i, 
				  const int n_i, 
				  const int n1_i, 
				  const double tau_d, double* sigma_k1, 				  
				  double* u_k1, double* v_k1, double* u2, double* v2, double* e_k, const double* e_k_mu)
{
	int c_u1;
	int c_u2;
	int c_v1;
	int c_v2;							
	const int break_value = (n1_i - 1) * (n_i - 1) - (2 + (q_i - 1) * 2) / 2 * q_i;

	mtn_calculate_common(gamma, m_i, m1_i, m2_i, qq_i, w_i, cntr_i, tau_d, sigma_k1, e_k, e_k_mu);

	int s_m = 0;
	for (s_m = 0; s_m <= 20; ++s_m)
	{
		mtn_calculate_jakobi(qq_i, m_i, w_i, cntr_i, m1_i, m2_i, u_k1, v_k1);

		c_u1 = 0;
		c_u2 = 0;
		c_v1 = 0;
		c_v2 = 0;
#pragma omp parallel
		{
#pragma omp for collapse(2) reduction(+:c_u1, c_v1) nowait
			for (int i = 1; i < m1_i - 1; i++)
			{
				for (int j = 1; j < m_i - 1; j++)
				{
					if (i < qq_i || i >= qq_i + w_i)
					{
						if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= epsilon_d)
						{
							++c_u1;
						}
						if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= epsilon_d)
						{
							++c_v1;
						}
					}
				}
			}

#pragma omp for reduction(+:c_u2, c_v2) nowait
			for (int i = qq_i; i < qq_i + w_i; i++)
			{
				for (int j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
				{
					if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= epsilon_d)
					{
						++c_u2;
					}
					if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= epsilon_d)
					{
						++c_v2;
					}
				}
				for (int j = cntr_i - i - 1 + qq_i; j > 0; j--)
				{
					if (fabs(u_k1[i * m_i + j] - u2[i * m_i + j]) <= epsilon_d)
					{
						++c_u2;
					}
					if (fabs(v_k1[i * m_i + j] - v2[i * m_i + j]) <= epsilon_d)
					{
						++c_v2;
					}
				}
			}
		} // #pragma omp parallel

		if (c_u1 + c_u2 == break_value && c_v1 + c_v2 >= break_value)
		{
			break;
		}

#pragma omp parallel
		{
#pragma omp for collapse(2) nowait
			for (int i = 1; i < m1_i - 1; i++)
			{
				for (int j = 1; j < m_i - 1; j++)
				{
					if (i < qq_i || i >= qq_i + w_i)
					{
						u_k1[i * m_i + j] = u2[i * m_i + j];
						v_k1[i * m_i + j] = v2[i * m_i + j];
					}
				}
			}

#pragma omp for nowait
			for (int i = qq_i; i < qq_i + w_i; i++)
			{
				for (int j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
				{
					u_k1[i * m_i + j] = u2[i * m_i + j];
					v_k1[i * m_i + j] = v2[i * m_i + j];
				}
				for (int j = cntr_i - i - 1 + qq_i; j > 0; j--)
				{
					u_k1[i * m_i + j] = u2[i * m_i + j];
					v_k1[i * m_i + j] = v2[i * m_i + j];
				}
			}
		}
	} // #pragma omp parallel
	return s_m;
}
