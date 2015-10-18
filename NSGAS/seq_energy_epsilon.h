/*----- Функция заполняет элементы матрицы для уравнения энергии----*/

inline double energy_A(double* Sigma_k1, double* e_k)
{
	int i = 0, j = 0, a = 0;


	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;
			A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
			A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				Gamma / (2 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j] + e_k_mu[(i - 1) * M + j]) +
				Gamma / (2 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1] + e_k_mu[i * M + j - 1]);
			A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
			A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));
		}
	}

	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr + i + 1 - qq;

		a = i * M + j;

		A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
		A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1])) - Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]) - Gamma / (8 * hy * hy * Pr * Re) * (e_k_mu[i * M + j - 1] + 2 * e_k_mu[i * M + j]);
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
			Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
			Gamma / (8 * hx * hx * Pr * Re) * (8 * e_k_mu[i * M + j] + 3 * e_k_mu[(i + 1) * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
			Gamma / (8 * hy * hy * Pr * Re) * (8 * e_k_mu[i * M + j] + 4 * e_k_mu[i * M + j + 1] + 3 * e_k_mu[i * M + j - 1]);
		A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
		A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j])) - Gamma / (4 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]) - Gamma / (8 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]);
	}

	i = qq + w - 1;
	j = cntr + i + 1 - qq;

	a = i * M + j;

	A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
	A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
		Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
		Gamma / (2 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j] + e_k_mu[(i - 1) * M + j]) +
		Gamma / (2 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1] + e_k_mu[i * M + j - 1]);
	A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
	A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr - i - 1 + qq;

		a = i * M + j;

		A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
		A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
			Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
			Gamma / (8 * hx * hx * Pr * Re) * (8 * e_k_mu[i * M + j] + 3 * e_k_mu[(i + 1) * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
			Gamma / (8 * hy * hy * Pr * Re) * (8 * e_k_mu[i * M + j] + 3 * e_k_mu[i * M + j + 1] + 4 * e_k_mu[i * M + j - 1]);
		A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j])) - Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]) - Gamma / (8 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]);
		A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j])) - Gamma / (4 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]) - Gamma / (8 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]);
	}


	i = qq + w - 1;
	j = cntr - i - 1 + qq;

	a = i * M + j;

	A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
	A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
		Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
		Gamma / (2 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j] + e_k_mu[(i - 1) * M + j]) +
		Gamma / (2 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1] + e_k_mu[i * M + j - 1]);
	A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
	A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 2 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
			A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				Gamma / (2 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j] + e_k_mu[(i - 1) * M + j]) +
				Gamma / (2 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1] + e_k_mu[i * M + j - 1]);
			A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
			A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));
		}
	}


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 2 + qq; j > 0; j--)
		{
			a = i * M + j;

			A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
			A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				Gamma / (2 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j] + e_k_mu[(i - 1) * M + j]) +
				Gamma / (2 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1] + e_k_mu[i * M + j - 1]);
			A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
			A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));
		}
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
			A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau - Gamma / (2 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[(i + 1) * M + j] - e_k[(i - 1) * M + j]) -
				Gamma / (2 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
				Gamma / (2 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j] + e_k_mu[(i - 1) * M + j]) +
				Gamma / (2 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j + 1] + e_k_mu[i * M + j - 1]);
			A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
			A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));
		}
	}


	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			/*A[a][0] = Gamma/(2*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) + (Mu(e_k[i*M+j])-Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(4*e_k[i*M+j]-2*e_k[(i-1)*M+j]-e_k[i*M+j-1]-e_k[i*M+j+1]) +
				Gamma/(4*h*h*Pr*Re)*(2*Mu(e_k[(i-1)*M+j])+Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[i*M+j+1])) );*/

	/*	A[a][0] = Gamma/(2*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(4*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-2*e_k[(i-1)*M+j]) -
                                                        Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[i*M+j-1]-e_k[i*M+j+1])+
                                                        Gamma/(4*hx*hx*Pr*Re)*(2*Mu(e_k[i*M+j])+2*Mu(e_k[(i-1)*M+j])) +
                                                        Gamma/(4*hy*hy*Pr*Re)*(2*Mu(e_k[i*M+j])+Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(4*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[i*M+j+1])) );

		}*/

	//Для Г2. l = q,...,n-1; m = 0;
	/*j = 0;
        for(i = 1; i < M1-1; i++)
        {
			a = i*M + j;

			/*A[a][0] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(4*e_k[i*M+j]-e_k[(i-1)*M+j]-e_k[(i+1)*M+j]-2*e_k[i*M+j+1]) +
				Gamma/(4*h*h*Pr*Re)*(2*Mu(e_k[i*M+j+1])+Mu(e_k[(i-1)*M+j])+Mu(e_k[(i+1)*M+j]));
			A[a][3] = - Gamma/(2*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])-Mu(e_k[i*M+j])) );
			A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[(i+1)*M+j])) );*/

	/*	A[a][0] = Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i-1)*M+j]-e_k[(i+1)*M+j]) -
                                                        Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-2*e_k[i*M+j+1])+
                                                        Gamma/(4*hx*hx*Pr*Re)*(2*Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])+Mu(e_k[(i+1)*M+j])) +
                                                        Gamma/(4*hy*hy*Pr*Re)*(2*Mu(e_k[i*M+j])+2*Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(2*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])+Mu(e_k[i*M+j])) );
			A[a][4] = - Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[(i+1)*M+j])) );

		}*/

	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			A[a][1] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(4*e_k[i*M+j]-2*e_k[(i+1)*M+j]-e_k[i*M+j+1]-e_k[i*M+j-1]) +
				Gamma/(4*h*h*Pr*Re)*(2*Mu(e_k[(i+1)*M+j])+Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[i*M+j+1])) );
			A[a][4] = - Gamma/(2*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[(i+1)*M+j])-Mu(e_k[i*M+j])) );
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			/*A[a][0] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(2*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) + (Mu(e_k[i*M+j])-Mu(e_k[i*M+j-1])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(4*e_k[i*M+j]-2*e_k[i*M+j-1]-e_k[(i+1)*M+j]-e_k[(i-1)*M+j]) +
				Gamma/(4*h*h*Pr*Re)*(2*Mu(e_k[i*M+j-1])+Mu(e_k[(i-1)*M+j])+Mu(e_k[(i+1)*M+j]));
			A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[(i+1)*M+j])) );*/

	/*	A[a][0] = Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(2*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j])+Mu(e_k[i*M+j-1])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i+1)*M+j]-e_k[(i-1)*M+j]) -
                                                        Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-2*e_k[i*M+j-1])+
                                                        Gamma/(4*hx*hx*Pr*Re)*(2*Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])+Mu(e_k[(i+1)*M+j])) +
                                                        Gamma/(4*hy*hy*Pr*Re)*(2*Mu(e_k[i*M+j])+2*Mu(e_k[i*M+j-1]));
			A[a][4] = - Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[(i+1)*M+j])) );
		}*/

	//Для Г5. l = q-1; m = 1,...,q-1;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < (cntr + q - 1); j++)
	{
		a = i * M + j;

		/*A[a][1] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(4*e_k[i*M+j]-2*e_k[(i+1)*M+j]-e_k[i*M+j+1]-e_k[i*M+j-1]) +
				Gamma/(4*h*h*Pr*Re)*(2*Mu(e_k[(i+1)*M+j])+Mu(e_k[i*M+j-1])+Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[i*M+j+1])) );
			A[a][4] = - Gamma/(2*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[(i+1)*M+j])-Mu(e_k[i*M+j])) );*/

		A[a][1] = Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] / (2 * tau) - Gamma / (4 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - 2 * e_k[(i + 1) * M + j]) -
			Gamma / (4 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (2 * e_k[i * M + j] - e_k[i * M + j + 1] - e_k[i * M + j - 1]) +
			Gamma / (4 * hx * hx * Pr * Re) * (2 * e_k_mu[i * M + j] + 2 * e_k_mu[(i + 1) * M + j]) +
			Gamma / (4 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j] + e_k_mu[i * M + j - 1] + e_k_mu[i * M + j + 1]);
		A[a][3] = - Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
		A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[(i + 1) * M + j] + e_k_mu[i * M + j]));
	}

	//Для Г6. l = 1,...,q-1; m = q-1;
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;

		a = i * M + j;

		A[a][0] = Gamma / (4 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
			- Gamma / (4 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[(i - 1) * M + j]) - Gamma / (8 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] + 2 * e_k_mu[(i - 1) * M + j]);
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[(i - 1) * M + j]) -
			Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[i * M + j + 1]) +
			Gamma / (8 * hx * hx * Pr * Re) * (4 * e_k_mu[i * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
			Gamma / (8 * hy * hy * Pr * Re) * (4 * e_k_mu[i * M + j] + 4 * e_k_mu[i * M + j + 1]);
		A[a][3] = - Gamma / (4 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) - Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j])
			- Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j + 1] + e_k_mu[i * M + j]) - Gamma / (8 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j + 1] + e_k_mu[i * M + j]);
		//A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[(i+1)*M+j])) );
		A[a][5] = -Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];
		A[a][8] = Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;

		a = i * M + j;

		A[a][0] = Gamma / (4 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
			- Gamma / (4 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[(i - 1) * M + j]) - Gamma / (8 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] + 2 * e_k_mu[(i - 1) * M + j]);
		A[a][1] = Gamma / (4 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i) * M + j - 1]) + Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i) * M + j - 1])
			- Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[(i) * M + j - 1]) - Gamma / (8 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] + 2 * e_k_mu[(i) * M + j - 1]);
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[(i - 1) * M + j]) -
			Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (4 * e_k[i * M + j] - 3 * e_k[i * M + j - 1]) +
			Gamma / (8 * hx * hx * Pr * Re) * (4 * e_k_mu[i * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
			Gamma / (8 * hy * hy * Pr * Re) * (4 * e_k_mu[i * M + j] + 4 * e_k_mu[i * M + j - 1]);
		//A[a][3] = - 3*Gamma/(8*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])+Mu(e_k[i*M+j])) );
		//A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[i*M+j])+Mu(e_k[(i+1)*M+j])) );
		A[a][6] = -Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];
		A[a][7] = Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];
	}

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i+1)*M+j]-e_k[i*M+j+1]) +
				Gamma/(4*h*h*Pr*Re)*(Mu(e_k[(i+1)*M+j])+Mu(e_k[i*M+j+1])-2*Mu(e_k[i*M+j]));
			A[a][3] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])-Mu(e_k[i*M+j])) );
			A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[(i+1)*M+j])-Mu(e_k[i*M+j])) );

*/
	//Для S_0q-1.
	/*		i = 0;
		j = q-1;

			a = i*M + j;

			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i+1)*M+j]-e_k[i*M+j+1]) +
				Gamma/(4*h*h*Pr*Re)*(Mu(e_k[(i+1)*M+j])+Mu(e_k[i*M+j+1])-2*Mu(e_k[i*M+j]));
			A[a][3] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])-Mu(e_k[i*M+j])) );
			A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[(i+1)*M+j])-Mu(e_k[i*M+j])) );

*/
	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			A[a][1] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) + (Mu(e_k[i*M+j])-Mu(e_k[i*M+j-1])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i+1)*M+j]-e_k[i*M+j-1]) +
				Gamma/(4*h*h*Pr*Re)*(Mu(e_k[(i+1)*M+j])+Mu(e_k[i*M+j-1])-2*Mu(e_k[i*M+j]));
			A[a][4] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[(i+1)*M+j]-e_k[i*M+j]) + (Mu(e_k[(i+1)*M+j])-Mu(e_k[i*M+j])) );

*/
	//Для S_nn.
	/*i = N1;
		j = N;

			a = i*M + j;

			/*A[a][0] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) + (Mu(e_k[i*M+j])-Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) + (Mu(e_k[i*M+j])-Mu(e_k[i*M+j-1])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i-1)*M+j]-e_k[i*M+j-1]) +
				Gamma/(4*h*h*Pr*Re)*(Mu(e_k[(i-1)*M+j])+Mu(e_k[i*M+j-1])-2*Mu(e_k[i*M+j]));*/

	/*          A[a][0] = Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(4*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j])+Mu(e_k[i*M+j-1])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) -
                                                        Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) +
                                                        Gamma/(4*hx*hx*Pr*Re)*(Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) +
                                                        Gamma/(4*hy*hy*Pr*Re)*(Mu(e_k[i*M+j])+Mu(e_k[i*M+j-1]));
*/
	/*A[a][0] = Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) + (Mu(e_k[i*M+j])-Mu(e_k[(i-1)*M+j])) );
			A[a][1] = Gamma/(4*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) - (Mu(e_k[i*M+j])+Mu(e_k[i*M+j-1])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) -
                                                        Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j-1]) +
                                                        Gamma/(4*hx*hx*Pr*Re)*(Mu(e_k[(i-1)*M+j])) +
                                                        Gamma/(4*hy*hy*Pr*Re)*(Mu(e_k[i*M+j-1]));*/


	//Для S_n0.
	/*i = N1;
		j = 0;

			a = i*M + j;

			/*A[a][0] = Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) + (Mu(e_k[i*M+j])-Mu(e_k[(i-1)*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*h*h*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(2*e_k[i*M+j]-e_k[(i-1)*M+j]-e_k[i*M+j+1]) +
				Gamma/(4*h*h*Pr*Re)*(Mu(e_k[(i-1)*M+j])+Mu(e_k[i*M+j+1])-2*Mu(e_k[i*M+j]));
			A[a][3] = - Gamma/(4*h*h*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])-Mu(e_k[i*M+j])) );*/

	/*			A[a][0] = Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) - (Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) -
                                                    Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j+1]) +
				Gamma/(4*hx*hx*Pr*Re)*(Mu(e_k[i*M+j])+Mu(e_k[(i-1)*M+j])) +
				Gamma/(4*hy*hy*Pr*Re)*(Mu(e_k[i*M+j])+Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(4*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])+Mu(e_k[i*M+j])) );
*/
	/*A[a][0] = Gamma/(4*hx*hx*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) + (Mu(e_k[i*M+j])-Mu(e_k[(i-1)*M+j])) );
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau)  - Gamma/(4*hx*hx*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[(i-1)*M+j]) -
                                                    Gamma/(4*hy*hy*Pr*Re)*Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j]-e_k[i*M+j+1]) +
				Gamma/(4*hx*hx*Pr*Re)*(Mu(e_k[(i-1)*M+j])) +
				Gamma/(4*hy*hy*Pr*Re)*(Mu(e_k[i*M+j+1]));
			A[a][3] = - Gamma/(4*hy*hy*Pr*Re)*( Mu(e_k[i*M+j])/e_k[i*M+j]*(e_k[i*M+j+1]-e_k[i*M+j]) + (Mu(e_k[i*M+j+1])+Mu(e_k[i*M+j])) );*/


	//Для S_qq,N/2+q.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;

	A[a][0] = Gamma / (4 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
		- Gamma / (4 * hx * hx * Pr * Re) * (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]) - Gamma / (8 * hx * hx * Pr * Re) * (2 * e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]);
	A[a][1] = Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (7 * e_k[i * M + j] - 4 * e_k[(i + 1) * M + j] - 3 * e_k[(i - 1) * M + j]) -
		Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (7 * e_k[i * M + j] - 4 * e_k[i * M + j + 1] - 2 * e_k[i * M + j - 1]) +
		Gamma / (8 * hx * hx * Pr * Re) * (7 * e_k_mu[i * M + j] + 4 * e_k_mu[(i + 1) * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
		Gamma / (8 * hy * hy * Pr * Re) * (7 * e_k_mu[i * M + j] + 4 * e_k_mu[i * M + j + 1] + 2 * e_k_mu[i * M + j - 1]) + Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[(i) * M + j];
	A[a][3] = - Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1])) /*+ 3*Gamma/(2*sqrt(10)*hx*hy*Pr*Re)*Mu(e_k[(i)*M+j])*/;
	A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));
	A[a][5] = -Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];


	//Для S_qq,N/2-q.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;

	A[a][0] = Gamma / (4 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) + Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j])
		- Gamma / (4 * hx * hx * Pr * Re) * (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]) - Gamma / (8 * hx * hx * Pr * Re) * (2 * e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]);
	A[a][1] = Gamma / (2 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) - (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]));
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (7 * e_k[i * M + j] - 4 * e_k[(i + 1) * M + j] - 3 * e_k[(i - 1) * M + j]) -
		Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (7 * e_k[i * M + j] - 2 * e_k[i * M + j + 1] - 4 * e_k[i * M + j - 1]) +
		Gamma / (8 * hx * hx * Pr * Re) * (7 * e_k_mu[i * M + j] + 4 * e_k_mu[(i + 1) * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
		Gamma / (8 * hy * hy * Pr * Re) * (7 * e_k_mu[i * M + j] + 2 * e_k_mu[i * M + j + 1] + 4 * e_k_mu[i * M + j - 1]) + Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[(i) * M + j];
	A[a][3] = - Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]));
	A[a][4] = - Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[(i + 1) * M + j] - e_k[i * M + j]) + (e_k_mu[i * M + j] + e_k_mu[(i + 1) * M + j]));
	A[a][6] = -Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];


	//Для S_qq,N/2
	i = qq;
	j = cntr;

	a = i * M + j;

	A[a][0] = Gamma / (2 * hx * hx * Pr * Re) * (e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[(i - 1) * M + j]) - (e_k_mu[(i - 1) * M + j] + e_k_mu[i * M + j]));
	A[a][1] = Gamma / (4 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1]) + Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j] - e_k[i * M + j - 1])
		- Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]) - Gamma / (8 * hy * hy * Pr * Re) * (2 * e_k_mu[i * M + j - 1] + e_k_mu[i * M + j]);
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] * (1 / (4 * tau) + 1 / (2 * tau)) - Gamma / (8 * hx * hx * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (6 * e_k[i * M + j] - 4 * e_k[(i - 1) * M + j]) -
		Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (6 * e_k[i * M + j] - 3 * e_k[i * M + j + 1] - 3 * e_k[i * M + j - 1]) +
		Gamma / (8 * hx * hx * Pr * Re) * (6 * e_k_mu[i * M + j] + 4 * e_k_mu[(i - 1) * M + j]) +
		Gamma / (8 * hy * hy * Pr * Re) * (6 * e_k_mu[i * M + j] + 4 * e_k_mu[i * M + j + 1] + 4 * e_k_mu[i * M + j - 1]) - Gamma / (tg * hx * hy * Pr * Re) * e_k_mu[(i) * M + j];
	A[a][3] = - Gamma / (4 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j]) - Gamma / (8 * hy * hy * Pr * Re) * e_k_mu[i * M + j] / e_k[i * M + j] * (e_k[i * M + j + 1] - e_k[i * M + j])
		- Gamma / (4 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] + e_k_mu[i * M + j + 1]) - Gamma / (8 * hy * hy * Pr * Re) * (e_k_mu[i * M + j] + 2 * e_k_mu[i * M + j + 1]);
	A[a][7] = Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];
	A[a][8] = Gamma / (2 * tg * hx * hy * Pr * Re) * e_k_mu[i * M + j];


	return 0;
}


//Вектор правых частей системы линейных уравнений
inline double energy_f(double* Sigma_k, double* Sigma_k1, double* u_k, double* v_k, double* u_k1, double* v_k1, double* e_k)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] / tau - P(Sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				e_k_mu[i * M + j] / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				e_k_mu[i * M + j] / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				e_k_mu[i * M + j] / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

				e_k_mu[i * M + j] / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] / tau - P(Sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				e_k_mu[i * M + j] / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				e_k_mu[i * M + j] / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				e_k_mu[i * M + j] / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

				e_k_mu[i * M + j] / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
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

			//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] / tau - P(Sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				e_k_mu[i * M + j] / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				e_k_mu[i * M + j] / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				e_k_mu[i * M + j] / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

				e_k_mu[i * M + j] / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] / tau - P(Sigma_k[a], e_k[i * M + j]) / (4 * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[(i - 1) * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

				e_k_mu[i * M + j] / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

				e_k_mu[i * M + j] / (6 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

				e_k_mu[i * M + j] / (8 * Re * e_k[i * M + j]) * ((v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
					(v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
					(v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

				e_k_mu[i * M + j] / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
					(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
					(u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
		}
	}


	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;



		/*	f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(2*tau) - P(Sigma_k[0*M+100],e_k[0*M+100])/(8*e_k[i*M+j])*((2*u_k1[i*M+j]-2*u_k1[(i-1)*M+j])/hx+(v_k1[i*M+j+1]-v_k1[i*M+j-1])/hy) +

				Mu(e_k[i*M+j])/(6*hx*hx*Re*e_k[i*M+j])*( (u_k1[i*M+j]-u_k1[(i-1)*M+j])*(u_k1[i*M+j]-u_k1[(i-1)*M+j]) ) +

				Mu(e_k[i*M+j])/(12*hy*hy*Re*e_k[i*M+j])*( (v_k1[i*M+j+1]-v_k1[i*M+j])*(v_k1[i*M+j+1]-v_k1[i*M+j]) + (v_k1[i*M+j]-v_k1[i*M+j-1])*(v_k1[i*M+j]-v_k1[i*M+j-1]) ) +

				Mu(e_k[i*M+j])/(8*Re*e_k[i*M+j])*( ( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j+1]/hy-u_k1[i*M+j]/hy )*( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j+1]/hy-u_k1[i*M+j]/hy ) +
													( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j]/hy-u_k1[i*M+j-1]/hy )*( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j]/hy-u_k1[i*M+j-1]/hy ) ) +

				Mu(e_k[i*M+j])/(12*Re*e_k[i*M+j])*( ( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy )*( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy ) +
													( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy )*( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy ) );

		}*/

	//Для Г2. l = 1,...,n-1; m = 0;
	/*j = 0;
		for(i = 1; i < M1-1; i++)
		{
			a = i*M + j;



		/*	f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(2*tau) - P(Sigma_k[a],e_k[i*M+j])/(8*e_k[i*M+j])*((u_k1[(i+1)*M+j]-u_k1[(i-1)*M+j])/hx+(2*v_k1[i*M+j+1]-2*v_k1[i*M+j])/hy) +

				Mu(e_k[i*M+j])/(12*hx*hx*Re*e_k[i*M+j])*( (u_k1[(i+1)*M+j]-u_k1[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]) + (u_k1[i*M+j]-u_k1[(i-1)*M+j])*(u_k1[i*M+j]-u_k1[(i-1)*M+j]) ) +

				Mu(e_k[i*M+j])/(6*hy*hy*Re*e_k[i*M+j])*( (v_k1[i*M+j+1]-v_k1[i*M+j])*(v_k1[i*M+j+1]-v_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(8*Re*e_k[i*M+j])*( ( (v_k1[(i+1)*M+j]-v_k1[i*M+j])/hx+(u_k1[i*M+j+1]-u_k1[i*M+j])/hy )*( (v_k1[(i+1)*M+j]-v_k1[i*M+j])/hx+(u_k1[i*M+j+1]-u_k1[i*M+j])/hy ) +
													( (v_k1[i*M+j]-v_k1[(i-1)*M+j])/hx+(u_k1[i*M+j+1]-u_k1[i*M+j])/hy )*( (v_k1[i*M+j]-v_k1[(i-1)*M+j])/hx+(u_k1[i*M+j+1]-u_k1[i*M+j])/hy )  ) +

				Mu(e_k[i*M+j])/(12*Re*e_k[i*M+j])*( ( u_k1[(i+1)*M+j]/hx-u_k1[i*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy )*( u_k1[(i+1)*M+j]/hx-u_k1[i*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy ) +
													( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy )*( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy )  );

		}*/


	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;



			eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(2*tau) - P(Sigma_k[a],e_k[i*M+j])/(8*h*e_k[i*M+j])*(2*u_k1[(i+1)*M+j]-2*u_k1[i*M+j]+v_k1[i*M+j+1]-v_k1[i*M+j-1]) +

				Mu(e_k[i*M+j])/(6*h*h*Re*e_k[i*M+j])*( (u_k1[(i+1)*M+j]-u_k1[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (v_k1[i*M+j+1]-v_k1[i*M+j])*(v_k1[i*M+j+1]-v_k1[i*M+j]) + (v_k1[i*M+j]-v_k1[i*M+j-1])*(v_k1[i*M+j]-v_k1[i*M+j-1]) ) +

				Mu(e_k[i*M+j])/(8*h*h*Re*e_k[i*M+j])*( ( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j+1]-u_k1[i*M+j] )*( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j+1]-u_k1[i*M+j] ) +
													( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j]-u_k1[i*M+j-1] )*( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j]-u_k1[i*M+j-1] ) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( ( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j+1]+v_k1[i*M+j] )*( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j+1]+v_k1[i*M+j] ) +
													( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j]+v_k1[i*M+j-1] )*( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j]+v_k1[i*M+j-1] ) );

		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;



		/*	f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(2*tau) - P(Sigma_k[a],e_k[i*M+j])/(8*e_k[i*M+j])*((u_k1[(i+1)*M+j]-u_k1[(i-1)*M+j])/hx+(2*v_k1[i*M+j]-2*v_k1[i*M+j-1])/hy) +

				Mu(e_k[i*M+j])/(12*hx*hx*Re*e_k[i*M+j])*( (u_k1[(i+1)*M+j]-u_k1[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]) + (u_k1[i*M+j]-u_k1[(i-1)*M+j])*(u_k1[i*M+j]-u_k1[(i-1)*M+j]) ) +

				Mu(e_k[i*M+j])/(6*hy*hy*Re*e_k[i*M+j])*( (v_k1[i*M+j]-v_k1[i*M+j-1])*(v_k1[i*M+j]-v_k1[i*M+j-1]) ) +

				Mu(e_k[i*M+j])/(8*Re*e_k[i*M+j])*( ( (v_k1[(i+1)*M+j]-v_k1[i*M+j])/hx+(u_k1[i*M+j]-u_k1[i*M+j-1])/hy )*( (v_k1[(i+1)*M+j]-v_k1[i*M+j])/hx+(u_k1[i*M+j]-u_k1[i*M+j-1])/hy ) +
													( (v_k1[i*M+j]-v_k1[(i-1)*M+j])/hx+(u_k1[i*M+j]-u_k1[i*M+j-1])/hy )*( (v_k1[i*M+j]-v_k1[(i-1)*M+j])/hx+(u_k1[i*M+j]-u_k1[i*M+j-1])/hy ) ) +

				Mu(e_k[i*M+j])/(12*Re*e_k[i*M+j])*( ( u_k1[(i+1)*M+j]/hx-u_k1[i*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy )*( u_k1[(i+1)*M+j]/hx-u_k1[i*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy ) +
													( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy )*( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy ) );

		}*/

	//Для Г5. l = q-1; m = 1,...,q-1;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < (cntr + q - 1); j++)
	{
		a = i * M + j;

		//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

		f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] / (2 * tau) - P(Sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * ((2 * u_k1[(i + 1) * M + j] - 2 * u_k1[i * M + j]) / hx + (v_k1[i * M + j + 1] - v_k1[i * M + j - 1]) / hy) +

			e_k_mu[i * M + j] / (6 * hx * hx * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j])) +

			e_k_mu[i * M + j] / (12 * hy * hy * Re * e_k[i * M + j]) * ((v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

			e_k_mu[i * M + j] / (8 * Re * e_k[i * M + j]) * (((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j + 1] - u_k1[i * M + j]) / hy) * ((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j + 1] - u_k1[i * M + j]) / hy) +
				((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j] - u_k1[i * M + j - 1]) / hy) * ((v_k1[(i + 1) * M + j] - v_k1[i * M + j]) / hx + (u_k1[i * M + j] - u_k1[i * M + j - 1]) / hy)) +

			e_k_mu[i * M + j] / (12 * Re * e_k[i * M + j]) * ((u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
				(u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));
	}

	//Для Г6. l = 1,...,q-1; m = q-1;
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;

		a = i * M + j;

		//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

		f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - P(Sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (u_k1[(i) * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j + 1] / hy - v_k1[i * M + j] / hy) -
			- P(Sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (-u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j + 1] / hy) +

			e_k_mu[i * M + j] / (24 * hx * hx * Re * e_k[i * M + j]) * (1 * (-u_k1[i * M + j]) * (-u_k1[i * M + j]) + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

			e_k_mu[i * M + j] / (24 * hy * hy * Re * e_k[i * M + j]) * (3 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + 1 * (v_k1[i * M + j]) * (v_k1[i * M + j])) +

			e_k_mu[i * M + j] / (16 * Re * e_k[i * M + j]) * (1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
				2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
				1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy)) +

			e_k_mu[i * M + j] / (24 * Re * e_k[i * M + j]) * (1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
				2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
				1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy));
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;

		a = i * M + j;

		//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

		f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] * (1 / (4 * tau) + 1 / (4 * tau)) - P(Sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (u_k1[(i) * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy - v_k1[i * M + j - 1] / hy) -
			- P(Sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (-u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j - 1] / hy) +

			e_k_mu[i * M + j] / (24 * hx * hx * Re * e_k[i * M + j]) * (1 * (-u_k1[i * M + j]) * (-u_k1[i * M + j]) + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

			e_k_mu[i * M + j] / (24 * hy * hy * Re * e_k[i * M + j]) * (3 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1]) + 1 * (-v_k1[i * M + j]) * (-v_k1[i * M + j])) +

			e_k_mu[i * M + j] / (16 * Re * e_k[i * M + j]) * (1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
				2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
				1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy)) +

			e_k_mu[i * M + j] / (24 * Re * e_k[i * M + j]) * (1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
				2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
				1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy));
	}

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(4*tau) - P(Sigma_k[a],e_k[i*M+j])/(8*h*e_k[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]+v_k1[i*M+j+1]-v_k1[i*M+j]) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (u_k1[(i+1)*M+j]-u_k1[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (v_k1[i*M+j+1]-v_k1[i*M+j])*(v_k1[i*M+j+1]-v_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(8*h*h*Re*e_k[i*M+j])*( ( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j+1]-u_k1[i*M+j] )*( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j+1]-u_k1[i*M+j] ) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( ( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j+1]+v_k1[i*M+j] )*( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j+1]+v_k1[i*M+j] ) );
*/
	//Для S_0q-1.
	/*		i = 0;
		j = q-1;

			a = i*M + j;

			eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(4*tau) - P(Sigma_k[a],e_k[i*M+j])/(8*h*e_k[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]+v_k1[i*M+j+1]-v_k1[i*M+j]) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (u_k1[(i+1)*M+j]-u_k1[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (v_k1[i*M+j+1]-v_k1[i*M+j])*(v_k1[i*M+j+1]-v_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(8*h*h*Re*e_k[i*M+j])*( ( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j+1]-u_k1[i*M+j] )*( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j+1]-u_k1[i*M+j] ) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( ( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j+1]+v_k1[i*M+j] )*( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j+1]+v_k1[i*M+j] ) );
*/
	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

			f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(4*tau) - P(Sigma_k[a],e_k[i*M+j])/(8*h*e_k[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]+v_k1[i*M+j]-v_k1[i*M+j-1]) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (u_k1[(i+1)*M+j]-u_k1[i*M+j])*(u_k1[(i+1)*M+j]-u_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( (v_k1[i*M+j]-v_k1[i*M+j-1])*(v_k1[i*M+j]-v_k1[i*M+j-1]) ) +

				Mu(e_k[i*M+j])/(8*h*h*Re*e_k[i*M+j])*( ( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j]-u_k1[i*M+j-1] )*( v_k1[(i+1)*M+j]-v_k1[i*M+j]+u_k1[i*M+j]-u_k1[i*M+j-1] ) ) +

				Mu(e_k[i*M+j])/(12*h*h*Re*e_k[i*M+j])*( ( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j]+v_k1[i*M+j-1] )*( u_k1[(i+1)*M+j]-u_k1[i*M+j]-v_k1[i*M+j]+v_k1[i*M+j-1] ) );

*/
	//Для S_nn.
	/*i = N1;
		j = N;

			a = i*M + j;


/*			f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(4*tau) - P(Sigma_k[0*M+100],e_k[0*M+100])/(8*e_k[i*M+j])*(u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx+v_k1[i*M+j]/hy-v_k1[i*M+j-1]/hy) +

				Mu(e_k[i*M+j])/(12*hx*hx*Re*e_k[i*M+j])*( (u_k1[i*M+j]-u_k1[(i-1)*M+j])*(u_k1[i*M+j]-u_k1[(i-1)*M+j]) ) +

				Mu(e_k[i*M+j])/(12*hy*hy*Re*e_k[i*M+j])*( (v_k1[i*M+j]-v_k1[i*M+j-1])*(v_k1[i*M+j]-v_k1[i*M+j-1]) ) +

				Mu(e_k[i*M+j])/(8*Re*e_k[i*M+j])*( ( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j]/hy-u_k1[i*M+j-1]/hy )*( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j]/hy-u_k1[i*M+j-1]/hy ) ) +

				Mu(e_k[i*M+j])/(12*Re*e_k[i*M+j])*( ( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy )*( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j]/hy+v_k1[i*M+j-1]/hy ) );
*/

	//Для S_n0.
	/*i = N1;
		j = 0;

			a = i*M + j;



/*			f[a] = eR_k[a] * Sigma_k1[a]*Sigma_k1[a]/(4*tau) - P(Sigma_k[0*M+100],e_k[0*M+100])/(8*e_k[i*M+j])*(u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx+v_k1[i*M+j+1]/hy-v_k1[i*M+j]/hy) +

				Mu(e_k[i*M+j])/(12*hx*hx*Re*e_k[i*M+j])*( (u_k1[i*M+j]-u_k1[(i-1)*M+j])*(u_k1[i*M+j]-u_k1[(i-1)*M+j]) ) +

				Mu(e_k[i*M+j])/(12*hy*hy*Re*e_k[i*M+j])*( (v_k1[i*M+j+1]-v_k1[i*M+j])*(v_k1[i*M+j+1]-v_k1[i*M+j]) ) +

				Mu(e_k[i*M+j])/(8*Re*e_k[i*M+j])*( ( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j+1]/hy-u_k1[i*M+j]/hy )*( v_k1[i*M+j]/hx-v_k1[(i-1)*M+j]/hx+u_k1[i*M+j+1]/hy-u_k1[i*M+j]/hy ) ) +

				Mu(e_k[i*M+j])/(12*Re*e_k[i*M+j])*( ( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy )*( u_k1[i*M+j]/hx-u_k1[(i-1)*M+j]/hx-v_k1[i*M+j+1]/hy+v_k1[i*M+j]/hy ) );
*/

	//Для S_qq,N/2+q.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;

	//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

	f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - P(Sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (2 * u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + 2 * v_k1[i * M + j + 1] / hy - v_k1[i * M + j] / hy - v_k1[i * M + j - 1] / hy) -
		P(Sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) +

		e_k_mu[i * M + j] / (24 * hx * hx * Re * e_k[i * M + j]) * (4 * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

		e_k_mu[i * M + j] / (24 * hy * hy * Re * e_k[i * M + j]) * (4 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + 2 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1]) + (v_k1[i * M + j]) * (v_k1[i * M + j])) +

		e_k_mu[i * M + j] / (16 * Re * e_k[i * M + j]) * (2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
			1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy)) +

		e_k_mu[i * M + j] / (24 * Re * e_k[i * M + j]) * (2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
			1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy));


	//Для S_qq,N/2-q.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;

	//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

	f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] * (3 / (4 * tau) + 1 / (8 * tau)) - P(Sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (2 * u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx + v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy - u_k1[(i - 1) * M + j] / hx - 2 * v_k1[i * M + j - 1] / hy) -
		P(Sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy) +

		e_k_mu[i * M + j] / (24 * hx * hx * Re * e_k[i * M + j]) * (4 * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) * (u_k1[(i + 1) * M + j] - u_k1[i * M + j]) + 3 * (u_k1[i * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[i * M + j] - u_k1[(i - 1) * M + j])) +

		e_k_mu[i * M + j] / (24 * hy * hy * Re * e_k[i * M + j]) * (4 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1]) + 2 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + (-v_k1[i * M + j]) * (-v_k1[i * M + j])) +

		e_k_mu[i * M + j] / (16 * Re * e_k[i * M + j]) * (2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			1 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx - u_k1[i * M + j] / hy) +
			2 * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i + 1) * M + j] / hx - v_k1[i * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
			2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

		e_k_mu[i * M + j] / (24 * Re * e_k[i * M + j]) * (2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			1 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j] / hy) +
			2 * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i + 1) * M + j] / hx - u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
			2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));


	//Для S_qq,N/2
	i = qq;
	j = cntr;

	a = i * M + j;

	//eR_k[a] = e_k[a] - trajectory(i, j, e_k, u_k[a], v_k[a]);

	f[a] = eR_k[a] * Sigma_k1[a] * Sigma_k1[a] * (1 / (4 * tau) + 1 / (2 * tau)) - P(Sigma_k[a], e_k[i * M + j]) / (8 * e_k[i * M + j]) * (2 * u_k1[(i) * M + j] / hx - 2 * u_k1[(i - 1) * M + j] / hx + v_k1[i * M + j + 1] / hy - v_k1[i * M + j - 1] / hy) -
		P(Sigma_k[a], e_k[i * M + j]) / (16 * e_k[i * M + j]) * (-2 * u_k1[i * M + j] / hx + v_k1[(i) * M + j + 1] / hy - v_k1[i * M + j - 1] / hy) +

		e_k_mu[i * M + j] / (12 * hx * hx * Re * e_k[i * M + j]) * (2 * (u_k1[(i) * M + j] - u_k1[(i - 1) * M + j]) * (u_k1[(i) * M + j] - u_k1[(i - 1) * M + j]) + (-u_k1[i * M + j]) * (-u_k1[i * M + j])) +

		e_k_mu[i * M + j] / (24 * hy * hy * Re * e_k[i * M + j]) * (3 * (v_k1[i * M + j + 1] - v_k1[i * M + j]) * (v_k1[i * M + j + 1] - v_k1[i * M + j]) + 3 * (v_k1[i * M + j] - v_k1[i * M + j - 1]) * (v_k1[i * M + j] - v_k1[i * M + j - 1])) +

		e_k_mu[i * M + j] / (16 * Re * e_k[i * M + j]) * (2 * (v_k1[(i) * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) * (v_k1[(i) * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j] / hy - u_k1[i * M + j - 1] / hy) +
			2 * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (v_k1[i * M + j] / hx - v_k1[(i - 1) * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			1 * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) * (-v_k1[i * M + j] / hx + u_k1[i * M + j + 1] / hy - u_k1[i * M + j] / hy) +
			1 * (-v_k1[i * M + j] / hx + u_k1[(i) * M + j] / hy - u_k1[i * M + j - 1] / hy) * (-v_k1[i * M + j] / hx + u_k1[(i) * M + j] / hy - u_k1[i * M + j - 1] / hy)) +

		e_k_mu[i * M + j] / (24 * Re * e_k[i * M + j]) * (2 * (u_k1[(i) * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (u_k1[(i) * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) +
			2 * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (u_k1[i * M + j] / hx - u_k1[(i - 1) * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j + 1] / hy + v_k1[i * M + j] / hy) +
			1 * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy) * (-u_k1[i * M + j] / hx - v_k1[i * M + j] / hy + v_k1[i * M + j - 1] / hy));


	return 0;
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
inline double energy_D()
{
	int i = 0, j = 0, a;

	for (i = 1; i < qq + 1; i++)
	{
		for (j = 1; j < (M - 1); j++)
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
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	return 0;
}

//Вектор B = A*Xk1
inline double energy_B(double* e_k1)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
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

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * e_k1[(i - 1) * M + j] + A[a][1] * e_k1[i * M + (j - 1)] + A[a][2] * e_k1[i * M + j] +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];
		}
	}

	//Для Г1. l = n; m = 1,...,n-1;
	/*		i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			B[a] = A[a][0]*e_k1[(i-1)*M+j] + A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j] +
					A[a][3]*e_k1[i*M+j+1];

		}
*/
	//Для Г2. l = 1,...,n-1; m = 0;
	/*		j = 0;
		for(i = 1; i < M1-1; i++)
		{
			a = i*M + j;

			B[a] = A[a][0]*e_k1[(i-1)*M+j] + A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] +
					A[a][4]*e_k1[(i+1)*M+j];
		}

*/
	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			B[a] = A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] +
					A[a][4]*e_k1[(i+1)*M+j];
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*		j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			B[a] = A[a][0]*e_k1[(i-1)*M+j] + A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j] +
					A[a][4]*e_k1[(i+1)*M+j];
		}
*/
	//Для Г5. l = q-1; m = 1,...,q-1;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < (cntr + q - 1); j++)
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

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			B[a] = A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] + A[a][4]*e_k1[(i+1)*M+j];
*/
	//Для S_0q-1.
	/*		i = 0;
		j = q-1;

			a = i*M + j;

			B[a] = A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] + A[a][4]*e_k1[(i+1)*M+j];

*/
	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			B[a] = A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j] + A[a][4]*e_k1[(i+1)*M+j];

*/
	//Для S_nn.
	/*
		i = N1;
		j = N;

			a = i*M + j;

			B[a] = A[a][0]*e_k1[(i-1)*M+j] + A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j];

*/
	//Для S_n0.
	/*
		i = N1;
		j = 0;

			a = i*M + j;

			B[a] = A[a][0]*e_k1[(i-1)*M+j] + A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1];
*/

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
inline double energy_Jakobi(double* e_k1, double* e2)
{
	int i = 0, j = 0, a;

	for (i = 1; i < qq + 1; i++)
	{
		for (j = 1; j < (M - 1); j++)
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
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			e2[i * M + j] = e_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	return 0;
}

//Метод Гаусса-Зейделя
inline double energy_Zeidel(double* e_k1, double* e2)
{
	int i = 0, j = 0, a;

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < (M - 1)))
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
	for (j = 1; j < (M - 1); j++)
	{
		a = i * M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > 0) && (j < (M - 1)))
		{
			if (j == cntr)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
					A[a][3] * e_k1[i * M + (j + 1)]
					+ A[a][7] * e_k1[(i + 1) * M + j - 1] + A[a][8] * e_k1[(i + 1) * M + j + 1];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			else
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}

		if (j == N)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] /*+ A[a][2]*e_k1[i*M+j]*/ +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq + 1; i < qq + w - 1; i++)
	{
		//for(j = cntr-i-1+qq ; j > 0; j--)
		for (j = 1; j <= cntr - i + qq; j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < (cntr - i + qq)))
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
					A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == cntr - i + qq)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + j - 1] /*+ A[a][2]*e_k1[i*M+j]*/
					+ A[a][6] * e2[(i - 1) * M + j + 1] + A[a][7] * e_k1[(i + 1) * M + j - 1];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	for (i = qq + 1; i < qq + w - 1; i++)
	{
		for (j = cntr + i - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			if (j == cntr + i - qq)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1]
					+ A[a][5] * e2[(i - 1) * M + j - 1] + A[a][8] * e_k1[(i + 1) * M + j + 1];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > cntr + i - qq) && (j < (M - 1)))
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


	i = qq + w - 1;
	for (j = 1; j < (M - 1); j++)
	{
		a = i * M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > 0) && (j < (cntr - w + 1)))
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == cntr - w + 1)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j]
				+ A[a][6] * e2[(i - 1) * M + j + 1];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j >= (cntr - q + 2)) && (j < (cntr + q - 1)))
		{
			B[a] = A[a][1] * e2[i * M + j - 1] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
				A[a][4] * e_k1[(i + 1) * M + j];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == cntr + w - 1)
		{
			B[a] = A[a][0] * e2[(i - 1) * M + j] + A[a][1] * e2[i * M + (j - 1)] /*+ A[a][2]*e_k1[i*M+j]*/ +
				A[a][3] * e_k1[i * M + (j + 1)] + A[a][4] * e_k1[(i + 1) * M + j]
				+ A[a][5] * e2[(i - 1) * M + j - 1];

			e2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > cntr + w - 1) && (j < N))
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

		/*else
            {
                B[a] = A[a][0]*e2[(i-1)*M+j] + A[a][1]*e2[i*M+(j-1)] /*+ A[a][2]*e_k1[i*M+j]*/ /*+
			/*		A[a][3]*e_k1[i*M+(j+1)] +  A[a][4]*e_k1[(i+1)*M+j];

                e2[i*M + j] = D[a]*(f[a] - B[a]);
            }*/
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * e2[(i - 1) * M + j] /*+ A[a][2]*e_k1[i*M+j]*/ + A[a][3] * e_k1[i * M + j + 1] +
					A[a][4] * e_k1[(i + 1) * M + j];

				e2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < (M - 1)))
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

	/*i = N1;
        for(j = 0; j < (M); j++)
		{
			a = i*M + j;

			if(j == 0)
			{
			    B[a] = A[a][0]*e2[(i-1)*M+j] //+ A[a][2]*e_k1[i*M+j]
			    + A[a][3]*e_k1[i*M+j+1];

                e2[i*M + j] = D[a]*(f[a] - B[a]);
			}

            if( (j > 0) && (j < (M-1)) )
            {
                B[a] = A[a][0]*e2[(i-1)*M+j] + A[a][1]*e2[i*M+j-1] //+ A[a][2]*e_k1[i*M+j]
                +	A[a][3]*e_k1[i*M+j+1];

                e2[i*M + j] = D[a]*(f[a] - B[a]);
            }

            if(j == N)
            {
                B[a] = A[a][0]*e2[(i-1)*M+j] + A[a][1]*e2[i*M+j-1] ;//+ A[a][2]*e_k1[i*M+j];

                e2[i*M + j] = D[a]*(f[a] - B[a]);
            }

		}*/


	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			B[a] = A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] +
					A[a][4]*e_k1[(i+1)*M+j];
		}
*/

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			B[a] = A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] + A[a][4]*e_k1[(i+1)*M+j];
*/
	//Для S_0q-1.
	/*		i = 0;
		j = q-1;

			a = i*M + j;

			B[a] = A[a][2]*e_k1[i*M+j] + A[a][3]*e_k1[i*M+j+1] + A[a][4]*e_k1[(i+1)*M+j];

*/
	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			B[a] = A[a][1]*e_k1[i*M+j-1] + A[a][2]*e_k1[i*M+j] + A[a][4]*e_k1[(i+1)*M+j];

*/


	return 0;
}


inline double energy(double* Sigma_k1, double* Sigma_k, double* u_k, double* v_k, double* u_k1, double* v_k1, double* e2, double* e_k, double* e_k1)
{
	int i = 0, j = 0, a, bl = 1, c;

	/*---------------------------------------------*/


	energy_A(Sigma_k1, e_k);
	energy_D();


	energy_f(Sigma_k, Sigma_k1, u_k, v_k, u_k1, v_k1, e_k);

	s_e = 0;
	//c = 0;
	while (bl)
	{
		energy_B(e_k1);
		energy_Jakobi(e_k1, e2);

		//energy_Zeidel(e_k1, e2);

		c = 0;

		//	if(s >= 1)
		//	{

		for (i = 1; i < qq + 1; i++)
		{
			for (j = 1; j < (M - 1); j++)
			{
				a = i * M + j;
				if ((fabs(e_k1[a] - e2[a]) <= epsilon))
				{
					c += 1;
					//bl = 0;
				}
			}
		}

		for (i = qq + 1; i < qq + w - 1; i++)
		{
			for (j = cntr + i - qq; j < M - 1; j++)
			{
				a = i * M + j;
				if ((fabs(e_k1[a] - e2[a]) <= epsilon))
				{
					c += 1;
					//bl = 0;
				}
			}
		}

		for (i = qq + 1; i < qq + w - 1; i++)
		{
			for (j = cntr - i + qq; j > 0; j--)
			{
				a = i * M + j;
				if ((fabs(e_k1[a] - e2[a]) <= epsilon))
				{
					c += 1;
					//bl = 0;
				}
			}
		}

		for (i = qq + w - 1; i < M1 - 1; i++)
		{
			for (j = 1; j < (M - 1); j++)
			{
				a = i * M + j;
				if ((fabs(e_k1[a] - e2[a]) <= epsilon))
				{
					c += 1;
					//bl = 0;
				}
			}
		}


		//if(c >= (M2 - (q)*(q)-(N-q)))
		if (c == (N1 - 1) * (N - 1) - (2 + (q - 2 - 1) * 2) / 2 * (q - 2))
		{
			bl = 0;
		}

		else if (s_e > 20)
		{
			//fprintf(out,"\nk = %i   t = %f   h = %f   s = %i\n\n",k,tau,h,s1);
			bl = 0;
		}

		else
		{
			for (i = 1; i < qq + 1; i++)
			{
				for (j = 1; j < (M - 1); j++)
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
				for (j = 1; j < (M - 1); j++)
				{
					a = i * M + j;
					e_k1[a] = e2[a];
				}
			}
		}

		s_e += 1;
		//if (s_e >= 1)
		//    bl = 0;
	}

	return 0;
}
