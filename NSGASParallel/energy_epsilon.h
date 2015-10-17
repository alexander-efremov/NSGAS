/*----- Функция заполняет элементы матрицы для уравнения энергии----*/
// This method requires refactoring: change i * m + j to a
// m = C_M
// qq = C_qq
// w_i = C_w
// cntr_i = C_cntr
inline void nrg_calc_matrix_a(const double gamma, const int m_i, const int m1_i, const int qq_i, const int w_i, const int cntr_i,
                              const int q_i, double* sigma_k1, double* e_k, const double* e_k_mu)
{
	const double c_coef1 = 2 * C_hx * C_hx * C_PrRe;
	const double c_coef2 = 2 * C_hy * C_hy * C_PrRe;
	const double c_coef3 = 4 * C_hy * C_hy * C_PrRe;
	const double c_coef4 = 8 * C_hy * C_hy * C_PrRe;
	const double c_coef5 = 2 * C_tg * C_hx * C_hy * C_PrRe;
	const double c_coef6 = 4 * C_hx * C_hx * C_PrRe;
	const double c_coef7 = 8 * C_hx * C_hx * C_PrRe;

#pragma omp parallel 
	{
#pragma omp for collapse(2) nowait
		for (int i = 1; i < m1_i - 1; i++)
		{
			for (int j = 1; j < m_i - 1; j++)
			{
				if (i < qq_i || i >= qq_i + w_i)
				{
					int a = i * m_i + j;
					A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
					A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
					A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
						gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
						gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j] + e_k_mu[(i - 1) * m_i + j]) +
						gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
					A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
					A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));
					D[a] = 1 / A[a][2];
					f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[(i - 1) * m_i + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +
						e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[a]) * (u_k1[(i + 1) * m_i + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +
						e_k_mu[a] / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +
						e_k_mu[a] / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
							(v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
							(v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
							(v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
						e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
							(u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
							(u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
							(u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
				}
			}
		}
#pragma omp for nowait 
		for (int i = qq_i; i < qq_i + w_i - 1; i++)
		{
			int j = cntr_i + i + 1 - qq_i;
			int a = i * m_i + j;
			A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
			A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1])) - gamma / c_coef3 * (e_k_mu[a - 1] + e_k_mu[a]) - gamma / c_coef4 * (e_k_mu[a - 1] + 2 * e_k_mu[a]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
				gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / c_coef7 * (8 * e_k_mu[a] + 3 * e_k_mu[(i + 1) * m_i + j] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (8 * e_k_mu[a] + 4 * e_k_mu[a + 1] + 3 * e_k_mu[a - 1]);
			A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a])) - gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]) - gamma / c_coef7 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]);

			j = cntr_i - i - 1 + qq_i;
			a = i * m_i + j;
			A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
			A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
				gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / c_coef7 * (8 * e_k_mu[a] + 3 * e_k_mu[(i + 1) * m_i + j] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (8 * e_k_mu[a] + 3 * e_k_mu[a + 1] + 4 * e_k_mu[a - 1]);
			A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a])) - gamma / c_coef3 * (e_k_mu[a] + e_k_mu[a + 1]) - gamma / c_coef4 * (2 * e_k_mu[a] + e_k_mu[a + 1]);
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a])) - gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]) - gamma / c_coef7 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]);
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i; i++)
		{
			for (int j = cntr_i + i + 2 - qq_i; j < m_i - 1; j++)
			{
				int a = i * m_i + j;
				A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
				A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
					gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j] + e_k_mu[(i - 1) * m_i + j]) +
					gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
				A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
				A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));
			}
			for (int j = cntr_i - i - 2 + qq_i; j > 0; j--)
			{
				int a = i * m_i + j;
				A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
				A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
				A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
					gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
					gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j] + e_k_mu[(i - 1) * m_i + j]) +
					gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
				A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
				A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));
			}
		}
#pragma omp for nowait
		for (int j = cntr_i - q_i + 2; j < cntr_i + q_i - 1; j++)
		{
			//Для Г5. l = C_q-1; m = 1,...,C_q-1;	
			int i = qq_i + w_i - 1;
			int a = i * m_i + j;
			A[a][1] = gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - gamma / c_coef6 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - 2 * e_k[(i + 1) * m_i + j]) -
				gamma / c_coef3 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / c_coef6 * (2 * e_k_mu[a] + 2 * e_k_mu[(i + 1) * m_i + j]) +
				gamma / c_coef3 * (2 * e_k_mu[a] + e_k_mu[a - 1] + e_k_mu[a + 1]);
			A[a][3] = -gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[(i + 1) * m_i + j] + e_k_mu[a]));
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / (2 * C_tau) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * ((2 * u_k1[(i + 1) * m_i + j] - 2 * u_k1[a]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

				e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[a]) * (u_k1[(i + 1) * m_i + j] - u_k1[a])) +

				e_k_mu[a] / (12 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

				e_k_mu[a] / (8 * C_Re * e_k[a]) * (((v_k1[(i + 1) * m_i + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) * ((v_k1[(i + 1) * m_i + j] - v_k1[a]) / C_hx + (u_k1[a + 1] - u_k1[a]) / C_hy) +
					((v_k1[(i + 1) * m_i + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy) * ((v_k1[(i + 1) * m_i + j] - v_k1[a]) / C_hx + (u_k1[a] - u_k1[a - 1]) / C_hy)) +

				e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					(u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		}
#pragma omp for nowait
		for (int i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			//Для Г6. l = 1,...,C_q-1; m = C_q-1;
			int j = cntr_i + i - qq_i;
			int a = i * m_i + j;
			A[a][0] = gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) + gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j])
				- gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i - 1) * m_i + j]) - gamma / c_coef7 * (e_k_mu[a] + 2 * e_k_mu[(i - 1) * m_i + j]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - gamma / c_coef7 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[(i - 1) * m_i + j]) -
				gamma / c_coef4 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[a + 1]) +
				gamma / c_coef7 * (4 * e_k_mu[a] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (4 * e_k_mu[a] + 4 * e_k_mu[a + 1]);
			A[a][3] = -gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) - gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a])
				- gamma / c_coef3 * (e_k_mu[a + 1] + e_k_mu[a]) - gamma / c_coef4 * (2 * e_k_mu[a + 1] + e_k_mu[a]);
			A[a][5] = -gamma / c_coef5 * e_k_mu[a];
			A[a][8] = gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a] / C_hy) -
				-P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a + 1] / C_hy) +

				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (1 * -u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +

				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 1 * v_k1[a] * v_k1[a]) +

				e_k_mu[a] / (16 * C_Re * e_k[a]) * (1 * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy)) +

				e_k_mu[a] / (24 * C_Re * e_k[a]) * (1 * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy));

			// Для Г7.
			j = cntr_i - i + qq_i;
			a = i * m_i + j;
			A[a][0] = gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) + gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j])
				- gamma / c_coef6 * (e_k_mu[a] + e_k_mu[(i - 1) * m_i + j]) - gamma / c_coef7 * (e_k_mu[a] + 2 * e_k_mu[(i - 1) * m_i + j]);
			A[a][1] = gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) + gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1])
				- gamma / c_coef3 * (e_k_mu[a] + e_k_mu[a - 1]) - gamma / c_coef4 * (e_k_mu[a] + 2 * e_k_mu[a - 1]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - gamma / c_coef7 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[(i - 1) * m_i + j]) -
				gamma / c_coef4 * e_k_mu[a] / e_k[a] * (4 * e_k[a] - 3 * e_k[a - 1]) +
				gamma / c_coef7 * (4 * e_k_mu[a] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (4 * e_k_mu[a] + 4 * e_k_mu[a - 1]);
			A[a][6] = -gamma / c_coef5 * e_k_mu[a];
			A[a][7] = gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (4 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				-P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a - 1] / C_hy) +

				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (1 * -u_k1[a] * -u_k1[a] + 3 * (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +

				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + 1 * -v_k1[a] * -v_k1[a]) +

				e_k_mu[a] / (16 * C_Re * e_k[a]) * (1 * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx - u_k1[a] / C_hy)) +

				e_k_mu[a] / (24 * C_Re * e_k[a]) * (1 * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a] / C_hy));
		}
#pragma omp single nowait
		{
			int i = qq_i + w_i - 1;
			int j = cntr_i + i + 1 - qq_i;
			int a = i * m_i + j;
			A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
			A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
				gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j] + e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
			A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));

			j = cntr_i - i - 1 + qq_i;
			a = i * m_i + j;
			A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
			A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] / C_tau - gamma / c_coef1 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[(i + 1) * m_i + j] - e_k[(i - 1) * m_i + j]) -
				gamma / c_coef2 * e_k_mu[a] / e_k[a] * (2 * e_k[a] - e_k[a + 1] - e_k[a - 1]) +
				gamma / c_coef1 * (2 * e_k_mu[a] + e_k_mu[(i + 1) * m_i + j] + e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef2 * (2 * e_k_mu[a] + e_k_mu[a + 1] + e_k_mu[a - 1]);
			A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));

			//Для S_qq, C_N / 2 + C_q.
			j = cntr_i + i - qq_i;
			a = i * m_i + j;
			A[a][0] = gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) + gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j])
				- gamma / c_coef6 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]) - gamma / c_coef7 * (2 * e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]);
			A[a][1] = gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - gamma / c_coef7 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 4 * e_k[(i + 1) * m_i + j] - 3 * e_k[(i - 1) * m_i + j]) -
				gamma / c_coef4 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 4 * e_k[a + 1] - 2 * e_k[a - 1]) +
				gamma / c_coef7 * (7 * e_k_mu[a] + 4 * e_k_mu[(i + 1) * m_i + j] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (7 * e_k_mu[a] + 4 * e_k_mu[a + 1] + 2 * e_k_mu[a - 1]) + gamma / c_coef5 * e_k_mu[a];
			A[a][3] = -gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));
			A[a][5] = -gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + 2 * v_k1[a + 1] / C_hy - v_k1[a] / C_hy - v_k1[a - 1] / C_hy) -
				P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a] / C_hy) +
				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (4 * (u_k1[(i + 1) * m_i + j] - u_k1[a]) * (u_k1[(i + 1) * m_i + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +
				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (4 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 2 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + v_k1[a] * v_k1[a]) +
				e_k_mu[a] / (16 * C_Re * e_k[a]) * (2 * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					2 * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy)) +
				e_k_mu[a] / (24 * C_Re * e_k[a]) * (2 * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					2 * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy));

			//Для S_qq, C_N/2 - C_q.
			j = cntr_i - i + qq_i;
			a = i * m_i + j;
			A[a][0] = gamma / c_coef6 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) + gamma / c_coef7 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j])
				- gamma / c_coef6 * (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]) - gamma / c_coef7 * (2 * e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]);
			A[a][1] = gamma / c_coef2 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) - (e_k_mu[a - 1] + e_k_mu[a]));
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - gamma / c_coef7 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 4 * e_k[(i + 1) * m_i + j] - 3 * e_k[(i - 1) * m_i + j]) -
				gamma / c_coef4 * e_k_mu[a] / e_k[a] * (7 * e_k[a] - 2 * e_k[a + 1] - 4 * e_k[a - 1]) +
				gamma / c_coef7 * (7 * e_k_mu[a] + 4 * e_k_mu[(i + 1) * m_i + j] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (7 * e_k_mu[a] + 2 * e_k_mu[a + 1] + 4 * e_k_mu[a - 1]) + gamma / c_coef5 * e_k_mu[a];
			A[a][3] = -gamma / c_coef3 * (e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) + (e_k_mu[a] + e_k_mu[a + 1]));
			A[a][4] = -gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[(i + 1) * m_i + j] - e_k[a]) + (e_k_mu[a] + e_k_mu[(i + 1) * m_i + j]));
			A[a][6] = -gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (3 / (4 * C_tau) + 1 / (8 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx + v_k1[a + 1] / C_hy + v_k1[a] / C_hy - u_k1[(i - 1) * m_i + j] / C_hx - 2 * v_k1[a - 1] / C_hy) -
				P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy) +
				e_k_mu[a] / (24 * C_hx * C_hx * C_Re * e_k[a]) * (4 * (u_k1[(i + 1) * m_i + j] - u_k1[a]) * (u_k1[(i + 1) * m_i + j] - u_k1[a]) + 3 * (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +
				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (4 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1]) + 2 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + -v_k1[a] * -v_k1[a]) +
				e_k_mu[a] / (16 * C_Re * e_k[a]) * (2 * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					1 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx - u_k1[a] / C_hy) +
					2 * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
				e_k_mu[a] / (24 * C_Re * e_k[a]) * (2 * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					1 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a] / C_hy) +
					2 * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));

			//Для S_qq,C_N/2
			i = qq_i;
			j = cntr_i;
			a = i * m_i + j;
			A[a][0] = gamma / c_coef1 * (e_k_mu[a] / e_k[a] * (e_k[a] - e_k[(i - 1) * m_i + j]) - (e_k_mu[(i - 1) * m_i + j] + e_k_mu[a]));
			A[a][1] = gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1]) + gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a] - e_k[a - 1])
				- gamma / c_coef3 * (e_k_mu[a - 1] + e_k_mu[a]) - gamma / c_coef4 * (2 * e_k_mu[a - 1] + e_k_mu[a]);
			A[a][2] = sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - gamma / c_coef7 * e_k_mu[a] / e_k[a] * (6 * e_k[a] - 4 * e_k[(i - 1) * m_i + j]) -
				gamma / c_coef4 * e_k_mu[a] / e_k[a] * (6 * e_k[a] - 3 * e_k[a + 1] - 3 * e_k[a - 1]) +
				gamma / c_coef7 * (6 * e_k_mu[a] + 4 * e_k_mu[(i - 1) * m_i + j]) +
				gamma / c_coef4 * (6 * e_k_mu[a] + 4 * e_k_mu[a + 1] + 4 * e_k_mu[a - 1]) - gamma / (C_tg * C_hx * C_hy * C_PrRe) * e_k_mu[a];
			A[a][3] = -gamma / c_coef3 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a]) - gamma / c_coef4 * e_k_mu[a] / e_k[a] * (e_k[a + 1] - e_k[a])
				- gamma / c_coef3 * (e_k_mu[a] + e_k_mu[a + 1]) - gamma / c_coef4 * (e_k_mu[a] + 2 * e_k_mu[a + 1]);
			A[a][7] = gamma / c_coef5 * e_k_mu[a];
			A[a][8] = gamma / c_coef5 * e_k_mu[a];
			f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] * (1 / (4 * C_tau) + 1 / (2 * C_tau)) - P(gamma, sigma_k[a], e_k[a]) / (8 * e_k[a]) * (2 * u_k1[a] / C_hx - 2 * u_k1[(i - 1) * m_i + j] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) -
				P(gamma, sigma_k[a], e_k[a]) / (16 * e_k[a]) * (-2 * u_k1[a] / C_hx + v_k1[a + 1] / C_hy - v_k1[a - 1] / C_hy) +
				e_k_mu[a] / (12 * C_hx * C_hx * C_Re * e_k[a]) * (2 * (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j]) + -u_k1[a] * -u_k1[a]) +
				e_k_mu[a] / (24 * C_hy * C_hy * C_Re * e_k[a]) * (3 * (v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + 3 * (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +
				e_k_mu[a] / (16 * C_Re * e_k[a]) * (2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
					2 * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					1 * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (-v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
					1 * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (-v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +
				e_k_mu[a] / (24 * C_Re * e_k[a]) * (2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
					2 * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					1 * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (-u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
					1 * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (-u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
		} // #pragma omp single		
	} // #pragma omp parallel

#pragma omp parallel 
	{
		//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
#pragma omp for nowait
		for (int i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			for (int j = cntr_i + i - qq_i; j < m_i - 1; j++)
			{
				D[i * m_i + j] = 1 / A[i * m_i + j][2];
			}
			for (int j = cntr_i - i + qq_i; j > 0; j--)
			{
				D[i * m_i + j] = 1 / A[i * m_i + j][2];
			}
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i; i++)
		{
			for (int j = cntr_i + i + 1 - qq_i; j < m_i - 1; j++)
			{
				int a = i * m_i + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[(i - 1) * m_i + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[a]) * (u_k1[(i + 1) * m_i + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +

					e_k_mu[a] / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					e_k_mu[a] / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
			for (int j = cntr_i - i - 1 + qq_i; j > 0; j--)
			{
				int a = i * m_i + j;
				f[a] = eR_k[a] * sigma_k1[a] * sigma_k1[a] / C_tau - P(gamma, sigma_k[a], e_k[a]) / (4 * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[(i - 1) * m_i + j]) / C_hx + (v_k1[a + 1] - v_k1[a - 1]) / C_hy) +

					e_k_mu[a] / (6 * C_hx * C_hx * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] - u_k1[a]) * (u_k1[(i + 1) * m_i + j] - u_k1[a]) + (u_k1[a] - u_k1[(i - 1) * m_i + j]) * (u_k1[a] - u_k1[(i - 1) * m_i + j])) +

					e_k_mu[a] / (6 * C_hy * C_hy * C_Re * e_k[a]) * ((v_k1[a + 1] - v_k1[a]) * (v_k1[a + 1] - v_k1[a]) + (v_k1[a] - v_k1[a - 1]) * (v_k1[a] - v_k1[a - 1])) +

					e_k_mu[a] / (8 * C_Re * e_k[a]) * ((v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a + 1] / C_hy - u_k1[a] / C_hy) +
						(v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[(i + 1) * m_i + j] / C_hx - v_k1[a] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) +
						(v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy) * (v_k1[a] / C_hx - v_k1[(i - 1) * m_i + j] / C_hx + u_k1[a] / C_hy - u_k1[a - 1] / C_hy)) +

					e_k_mu[a] / (12 * C_Re * e_k[a]) * ((u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a + 1] / C_hy + v_k1[a] / C_hy) +
						(u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[(i + 1) * m_i + j] / C_hx - u_k1[a] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) +
						(u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy) * (u_k1[a] / C_hx - u_k1[(i - 1) * m_i + j] / C_hx - v_k1[a] / C_hy + v_k1[a - 1] / C_hy));
			}
		}
	}
}

//Вектор B = A*Xk1
// m = C_M
// qq_i = C_qq
// w_i = C_w
// cntr_i = C_cntr
// m1_i = C_M1
inline void nrg_calc_energy(double* e_k1, const int m, const int qq_i, const int q_i, const int w_i, const int cntr_i, const int m1_i)
{
#pragma omp parallel 
	{
		//Для внутренних узлов
#pragma omp for collapse(2) nowait
		for (int i = 1; i < m1_i - 1; i++)
		{
			for (int j = 1; j < m - 1; j++)
			{
				if (i < qq_i || i >= qq_i + w_i)
				{
					int a = i * m + j;
					B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
						A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
					e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
				}
			}
		}
#pragma omp for nowait
		for (int i = qq_i; i < qq_i + w_i; i++)
		{
			for (int j = cntr_i + i + 1 - qq_i; j < m - 1; j++)
			{
				int a = i * m + j;
				B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
			}
			for (int j = cntr_i - i - 1 + qq_i; j > 0; j--)
			{
				int a = i * m + j;
				B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + (j - 1)] + A[a][2] * e_k1[a] +
					A[a][3] * e_k1[i * m + (j + 1)] + A[a][4] * e_k1[(i + 1) * m + j];
			}
		}
#pragma omp for nowait
		for (int j = cntr_i - q_i + 2; j < cntr_i + q_i - 1; j++)
		{
			//Для Г5. l = C_q-1; m = 1,...,C_q-1;
			int i = qq_i + w_i - 1;
			int a = i * m + j;
			B[a] = A[a][1] * e_k1[a - 1] + A[a][2] * e_k1[i * m + j] + A[a][3] * e_k1[a + 1] +
				A[a][4] * e_k1[(i + 1) * m + j];
		}
#pragma omp for nowait
		for (int i = qq_i + 1; i < qq_i + w_i - 1; i++)
		{
			//Для Г6. l = 1,...,C_q-1; m = C_q-1;
			int j = cntr_i + i - qq_i;
			int a = i * m + j;
			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][2] * e_k1[i * m + j] + A[a][3] * e_k1[a + 1]
				+ A[a][5] * e_k1[(i - 1) * m + j - 1] + A[a][8] * e_k1[(i + 1) * m + j + 1];
			//Для Г7.
			j = cntr_i - i + qq_i;
			a = i * m + j;
			B[a] = A[a][0] * e_k1[(i - 1) * m + j] + A[a][1] * e_k1[i * m + j - 1] + A[a][2] * e_k1[a]
				+ A[a][6] * e_k1[(i - 1) * m + j + 1] + A[a][7] * e_k1[(i + 1) * m + j - 1];
		}
#pragma omp single
		{
			//Для S_qq,C_N/2+C_q.
			int i = qq_i + w_i - 1;
			int j = cntr_i + i - qq_i;
			int a = i * m + j;
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
	} // #pragma omp parallel 
#pragma omp parallel for
	//Метод Якоби
	for (int i = qq_i + 1; i < qq_i + w_i - 1; i++)
	{
		for (int j = cntr_i + i - qq_i; j < m - 1; j++)
		{
			e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
		}
		for (int j = cntr_i - i + qq_i; j > 0; j--)
		{
			e2[i * m + j] = e_k1[i * m + j] - D[i * m + j] * (B[i * m + j] - f[i * m + j]);
		}
	}
}

// m = C_M
// n = C_N
// qq_i = C_qq
// w_i = C_w
// m1 = C_M1
// q_i = C_q
inline int energy(const double gamma,
                  double* sigma_k1, double* e2, double* e_k, double* e_k1,
                  const int m_i,
                  const int n, const int qq_i,
                  const int w_i, const int m1_i,
                  const int n1_i,
                  const double epsilon_d,
                  const int q_i, const int cntr_i, const double* e_k_mu)
{
	int c;
	const int c_br = (n1_i - 1) * (n - 1);
	nrg_calc_matrix_a(gamma, m_i, m1_i, qq_i, w_i, cntr_i, q_i, sigma_k1, e_k, e_k_mu);
	int s_e = 0;
	for (s_e = 0; s_e <= 20; ++s_e)
	{
		nrg_calc_energy(e_k1, m_i, qq_i, q_i, w_i, cntr_i, m1_i);
		c = 0;

#pragma omp parallel for collapse(2) reduction(+:c)
		for (int i = 1; i < m1_i - 1; i++)
		{
			for (int j = 1; j < m_i - 1; j++)
			{
				if (fabs(e_k1[i * m_i + j] - e2[i * m_i + j]) <= epsilon_d)
				{
					++c;
				}
			}
		}

		if (c == c_br)
		{
			break;
		}

#pragma omp parallel for collapse(2)
		for (int i = 1; i < m1_i - 1; i++)
		{
			for (int j = 1; j < m_i - 1; j++)
			{
				e_k1[i * m_i + j] = e2[i * m_i + j];
			}
		}
	}

	return s_e;
}
