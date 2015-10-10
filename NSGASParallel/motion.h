/*----- Функция заполняет элементы матрицы, составленной для двух уравнении движения.----*/

inline double motion_A(double* Sigma_k1, double* e_k)
{
	int i = 0, j = 0, a;

	///////////////////////////////////////////////////Уравнение для u

	//Для внутренних узлов.
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
			A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][4] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re);
			A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);
		}
	}


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr + i + 1 - qq;

		a = i * M + j;

		A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
		A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (4 * hy * hy * Re) - (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j])) / (8 * hy * hy * Re);
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re)
			+ (Mu(e_k[(i + 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[(i + 1) * M + j]) + 2 * Mu(e_k[i * M + j])) / (6 * hx * hx * Re)
			+ (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j - 1)])) / (4 * hy * hy * Re) + (2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j - 1)])) / (8 * hy * hy * Re);
		A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
		A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) - (2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (6 * hx * hx * Re);
		A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
		A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
		//A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[(i+1)*M+j])/6 )/(h*h*Re);
		A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);
		A[a][9] = ((1. / 4. - 1. / 8.) * Mu(e_k[i * M + (j - 1)])) / (hx * hy * Re);
		A[a][11] = ((1. / 12. - 1. / 6.) * Mu(e_k[(i + 1) * M + j])) / (hx * hy * Re);
	}

	i = qq + w - 1;
	j = cntr + i + 1 - qq;

	a = i * M + j;

	A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
	A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
	A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
	A[a][4] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re);
	A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
	A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
	A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
	A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr - i - 1 + qq;

		a = i * M + j;

		A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
		A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
		A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re)
			+ (Mu(e_k[(i + 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[(i + 1) * M + j]) + 2 * Mu(e_k[i * M + j])) / (6 * hx * hx * Re)
			+ (Mu(e_k[i * M + (j + 1)]) + Mu(e_k[i * M + j])) / (4 * hy * hy * Re) + (Mu(e_k[i * M + (j + 1)]) + 2 * Mu(e_k[i * M + j])) / (8 * hy * hy * Re);
		A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (4 * hy * hy * Re) - (2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (8 * hy * hy * Re);
		A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) - (2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (6 * hx * hx * Re);
		A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
		A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
		A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
		//A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(h*h*Re);
		A[a][10] = ((1. / 8. - 1. / 4.) * Mu(e_k[i * M + (j + 1)])) / (hx * hy * Re);
		A[a][11] = ((1. / 6. - 1. / 12.) * Mu(e_k[(i + 1) * M + j])) / (hx * hy * Re);
	}

	i = qq + w - 1;
	j = cntr - i - 1 + qq;

	a = i * M + j;

	A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
	A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
	A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
	A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
	A[a][4] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re);
	A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
	A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
	A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
	A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 2 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
			A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][4] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re);
			A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 2 + qq; j > 0; j--)
		{
			a = i * M + j;

			A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
			A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][4] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re);
			A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			A[a][0] = - 2 * (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (3 * hx * hx * Re);
			A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (2 * hy * hy * Re);
			A[a][2] = Sigma_k1[a] * Sigma_k1[a] / tau + 2 * (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re) + (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (2 * hy * hy * Re);
			A[a][4] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (3 * hx * hx * Re);
			A[a][5] = (Mu(e_k[(i - 1) * M + j]) / 6 - Mu(e_k[i * M + (j - 1)]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[i * M + j + 1]) / 4 - Mu(e_k[(i - 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[i * M + j - 1]) / 4 - Mu(e_k[(i + 1) * M + j]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[(i + 1) * M + j]) / 6 - Mu(e_k[i * M + (j + 1)]) / 4) / (hx * hy * Re);
		}
	}


	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			/*A[a][0] = - 2*(Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + 2*(Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(4*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(4*h*h*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(h*h*Re);*/

	/*A[a][0] = - 2*(Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + 2*(Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(4*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(4*hy*hy*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(hx*hy*Re);*/


	/*A[a][0] = - 2*(Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + 2*(Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(4*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(4*hy*hy*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 + Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][8] = (- Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(hx*hy*Re);
		}*/

	//Для Г2. l = 1,...,n-1; m = 0;
	/*j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Ro(i,j,d)/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+j+1])-Mu(e_k[i*M+j]))/(2*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(2*h*h*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i+1)*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+j+1])+Mu(e_k[i*M+j]))/(2*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) + Mu(e_k[i*M+j]))/(2*hy*hy*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*hx*hx*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i+1)*M+j])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(hx*hy*Re);
		}*/

	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			A[a][1] = - (Mu(e_k[i*M+j-1]) + Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + 2*(Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+j-1]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[i*M+j+1]))/(4*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+j+1]))/(4*h*h*Re);
			A[a][4] = - 2*(Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+j-1])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[(i+1)*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(h*h*Re);
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(2*h*h*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(2*h*h*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i-1)*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[(i+1)*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(2*hy*hy*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(2*hy*hy*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(3*hx*hx*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[(i+1)*M+j])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+j])/4 )/(hx*hy*Re);
		}*/

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			A[a][2] = Ro(i,j,d)/(4*tau) + (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][4] = - (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][5] = - Mu(e_k[i*M+j])/(12*h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i+1)*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(h*h*Re);
*/

	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			A[a][1] = - (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau) + (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][4] = - (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(h*h*Re);
			A[a][6] = Mu(e_k[i*M+j])/(12*h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[(i+1)*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[(i+1)*M+j])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);
*/

	//Для S_nn.
	/*i = N1;
		j = N;

			a = i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau) + (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i-1)*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][8] = - Mu(e_k[i*M+j])/(12*h*h*Re);*/

	/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 - Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][8] = - Mu(e_k[i*M+j])/(12*hx*hy*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+(j-1)])/4 )/(hx*hy*Re);
			A[a][6] = ( - Mu(e_k[i*M+j])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j-1])/4 + Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][8] = Mu(e_k[i*M+j])/(12*hx*hy*Re);
*/

	//Для S_n0.
	/*i = N1;
		j = 0;

			a = i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Ro(i,j,d)/(4*tau) + (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(3*h*h*Re) + (Mu(e_k[i*M+(j+1)]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(h*h*Re);
			A[a][7] = Mu(e_k[i*M+j])/(12*h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re);
			A[a][2] = Sigma_k1[a]*Sigma_k1[a]/(4*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(3*hx*hx*Re) + (Mu(e_k[i*M+(j+1)]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) + Mu(e_k[i*M+j]))/(4*hy*hy*Re);
			A[a][5] = ( Mu(e_k[(i-1)*M+j])/6 + Mu(e_k[i*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[i*M+j+1])/4 - Mu(e_k[(i-1)*M+j])/6 )/(hx*hy*Re);
			A[a][7] = - Mu(e_k[i*M+j])/(12*hx*hy*Re);
			A[a][8] = ( - Mu(e_k[i*M+j])/6 - Mu(e_k[i*M+(j+1)])/4 )/(hx*hy*Re);
*/

	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов. l,m = 1,...,n-1
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
			A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
			A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re);
			A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);
		}
	}

	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr + i + 1 - qq;

		a = M2 + i * M + j;

		A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
		A[a][1] = - (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re) - (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j])) / (6 * hy * hy * Re);
		A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re)
			+ (Mu(e_k[(i + 1) * M + j]) + Mu(e_k[i * M + j])) / (4 * hx * hx * Re) + (Mu(e_k[(i + 1) * M + j]) + 2 * Mu(e_k[i * M + j])) / (8 * hx * hx * Re)
			+ (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j - 1)])) / (3 * hy * hy * Re) + (2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j - 1)])) / (6 * hy * hy * Re);
		A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
		A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (4 * hx * hx * Re) - (2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (8 * hx * hx * Re);
		A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
		A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
		//A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(h*h*Re);
		A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);
		A[a][9] = ((1. / 12. - 1. / 6.) * Mu(e_k[i * M + (j - 1)])) / (hx * hy * Re);
		A[a][11] = ((1. / 4. - 1. / 8.) * Mu(e_k[(i + 1) * M + j])) / (hx * hy * Re);
	}

	i = qq + w - 1;
	j = cntr + i + 1 - qq;

	a = M2 + i * M + j;

	A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
	A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
	A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
	A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
	A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re);
	A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
	A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
	A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
	A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr - i - 1 + qq;

		a = M2 + i * M + j;

		A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
		A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
		A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re)
			+ (Mu(e_k[(i + 1) * M + j]) + Mu(e_k[i * M + j])) / (4 * hx * hx * Re) + (Mu(e_k[(i + 1) * M + j]) + 2 * Mu(e_k[i * M + j])) / (8 * hx * hx * Re)
			+ (Mu(e_k[i * M + (j + 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re) + (Mu(e_k[i * M + (j + 1)]) + 2 * Mu(e_k[i * M + j])) / (6 * hy * hy * Re);
		A[a][3] = - (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re) - (2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (6 * hy * hy * Re);
		A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (4 * hx * hx * Re) - (2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (8 * hx * hx * Re);
		A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
		A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
		A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
		//A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[(i+1)*M+j])/4 )/(h*h*Re);
		A[a][10] = ((1. / 6. - 1. / 12.) * Mu(e_k[i * M + (j + 1)])) / (hx * hy * Re);
		A[a][11] = ((1. / 8. - 1. / 4.) * Mu(e_k[(i + 1) * M + j])) / (hx * hy * Re);
	}

	i = qq + w - 1;
	j = cntr - i - 1 + qq;

	a = M2 + i * M + j;

	A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
	A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
	A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
	A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
	A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re);
	A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
	A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
	A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
	A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 2 - qq; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
			A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
			A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re);
			A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 2 + qq; j > 0; j--)
		{
			a = M2 + i * M + j;

			A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
			A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
			A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re);
			A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);
		}
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			A[a][0] = - (Mu(e_k[(i - 1) * M + j]) + Mu(e_k[i * M + j])) / (2 * hx * hx * Re);
			A[a][1] = - 2 * (Mu(e_k[i * M + (j - 1)]) + Mu(e_k[i * M + j])) / (3 * hy * hy * Re);
			A[a][2] = Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau + (Mu(e_k[(i - 1) * M + j]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re) + 2 * (Mu(e_k[i * M + (j - 1)]) + 2 * Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][3] = - 2 * (Mu(e_k[i * M + j]) + Mu(e_k[i * M + (j + 1)])) / (3 * hy * hy * Re);
			A[a][4] = - (Mu(e_k[i * M + j]) + Mu(e_k[(i + 1) * M + j])) / (2 * hx * hx * Re);
			A[a][5] = (Mu(e_k[i * M + j - 1]) / 6 - Mu(e_k[(i - 1) * M + j]) / 4) / (hx * hy * Re);
			A[a][6] = (Mu(e_k[(i - 1) * M + j]) / 4 - Mu(e_k[i * M + j + 1]) / 6) / (hx * hy * Re);
			A[a][7] = (Mu(e_k[(i + 1) * M + j]) / 4 - Mu(e_k[i * M + j - 1]) / 6) / (hx * hy * Re);
			A[a][8] = (Mu(e_k[i * M + j + 1]) / 6 - Mu(e_k[(i + 1) * M + j]) / 4) / (hx * hy * Re);
		}
	}


	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(2*h*h*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(2*h*h*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(3*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);*/

	/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(2*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(2*hx*hx*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(3*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(3*hy*hy*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[i*M+j])/4 )/(hx*hy*Re);*/


	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(2*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(2*hx*hx*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(3*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(3*hy*hy*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(hx*hy*Re);
			A[a][7] = (- Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 + Mu(e_k[i*M+j])/4 )/(hx*hy*Re);
		}*/

	//Для Г2. l = 1,...,n-1; m = 0;
	/*j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][2] = Ro(i,j,d)/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*h*h*Re) + 2*(Mu(e_k[i*M+(j+1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][3] = - 2*(Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[(i-1)*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[(i+1)*M+j])/4 )/(h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*hx*hx*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*hx*hx*Re) + 2*(Mu(e_k[i*M+(j+1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][3] = - 2*(Mu(e_k[i*M+j+1]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*hx*hx*Re);
			A[a][5] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[(i-1)*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[(i+1)*M+j])/4 )/(hx*hy*Re);
		}*/

	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(2*h*h*Re) + (Mu(e_k[i*M+(j-1)])+2*Mu(e_k[i*M+j])+Mu(e_k[i*M+(j+1)]))/(3*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j]) + Mu(e_k[i*M+(j+1)]))/(3*h*h*Re);
			A[a][4] = - (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(2*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[(i+1)*M+j])/4 )/(h*h*Re);
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][1] = - 2*(Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*h*h*Re) + 2*(Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[(i+1)*M+j])/4 )/(h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*hx*hx*Re);
			A[a][1] = - 2*(Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) + (Mu(e_k[(i-1)*M+j]) + 2*Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*hx*hx*Re) + 2*(Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][4] = - (Mu(e_k[i*M+j]) + Mu(e_k[(i+1)*M+j]))/(4*hx*hx*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(hx*hy*Re);
			A[a][8] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[(i+1)*M+j])/4 )/(hx*hy*Re);
		}*/

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = M2 + i*M + j;

			A[a][2] = Ro(i,j,d)/(4*tau) + (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re) + (Mu(e_k[i*M+(j+1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][4] = - (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][5] = - Mu(e_k[i*M+j])/(12*h*h*Re);
			A[a][6] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[(i+1)*M+j])/4 )/(h*h*Re);

*/
	//Для S_0n.
	/*		i = 0;
		j = N;

			a = M2 + i*M + j;

			A[a][1] = - (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) + (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re) + (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][4] = - (Mu(e_k[(i+1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);
			A[a][6] = Mu(e_k[i*M+j])/(12*h*h*Re);
			A[a][7] = ( Mu(e_k[(i+1)*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[(i+1)*M+j])/4 )/(h*h*Re);

*/
	//Для S_nn.
	/*i = N1;
		j = N;

			a = M2 + i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) + (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re) + (Mu(e_k[i*M+(j-1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j])/6 )/(h*h*Re);
			A[a][7] = ( Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(h*h*Re);
			A[a][8] = - Mu(e_k[i*M+j])/(12*h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*hx*hx*Re);
			A[a][1] = - (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*hx*hx*Re) + (Mu(e_k[i*M+(j-1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][5] = ( Mu(e_k[i*M+j-1])/6 - Mu(e_k[(i-1)*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 + Mu(e_k[i*M+j])/6 )/(hx*hy*Re);
			A[a][7] = ( - Mu(e_k[i*M+j])/4 - Mu(e_k[i*M+j-1])/6 )/(hx*hy*Re);
			A[a][8] =  Mu(e_k[i*M+j])/(12*hx*hy*Re);
*/

	//Для S_n0.
	/*i = N1;
		j = 0;

			a = M2 + i*M + j;

			/*A[a][0] = - (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re);
			A[a][2] = Ro(i,j,d)/(4*tau) + (Mu(e_k[(i-1)*M+j]) - Mu(e_k[i*M+j]))/(4*h*h*Re) + (Mu(e_k[i*M+(j+1)]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) - Mu(e_k[i*M+j]))/(3*h*h*Re);
			A[a][5] = ( Mu(e_k[i*M+j])/6 - Mu(e_k[(i-1)*M+j])/4 )/(h*h*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(h*h*Re);
			A[a][7] = Mu(e_k[i*M+j])/(12*h*h*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 - Mu(e_k[i*M+j])/4 )/(h*h*Re);*/

	/*	A[a][0] = - (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*hx*hx*Re);
			A[a][2] = Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) + (Mu(e_k[(i-1)*M+j]) + Mu(e_k[i*M+j]))/(4*hx*hx*Re) + (Mu(e_k[i*M+(j+1)]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][3] = - (Mu(e_k[i*M+j+1]) + Mu(e_k[i*M+j]))/(3*hy*hy*Re);
			A[a][5] = ( - Mu(e_k[i*M+j])/6 - Mu(e_k[(i-1)*M+j])/4 )/(hx*hy*Re);
			A[a][6] = ( Mu(e_k[(i-1)*M+j])/4 - Mu(e_k[i*M+j+1])/6 )/(hx*hy*Re);
			A[a][7] = - Mu(e_k[i*M+j])/(12*hx*hy*Re);
			A[a][8] = ( Mu(e_k[i*M+j+1])/6 + Mu(e_k[i*M+j])/4 )/(hx*hy*Re);
*/

	return 0;
}


//Вектор правых частей системы уравнений
inline double motion_f(double* Sigma_k, double* Sigma_k1, double* u_k, double* v_k, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u

	//Для внутренних узлов. l,m = 1,...,n-1
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			//uX_k[a] = u_k[a] - trajectory(i, j, u_k, u_k[a], v_k[a]);


			f[a] = uX_k[a] * Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau - (P(Sigma_k[(i + 1) * M + j], e_k[(i + 1) * M + j]) - P(Sigma_k[(i - 1) * M + j], e_k[(i - 1) * M + j])) / (2 * hx);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			//uX_k[a] = u_k[a] - trajectory(i, j, u_k, u_k[a], v_k[a]);

			f[a] = uX_k[a] * Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau - (P(Sigma_k[(i + 1) * M + j], e_k[(i + 1) * M + j]) - P(Sigma_k[(i - 1) * M + j], e_k[(i - 1) * M + j])) / (2 * hx);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;

			//uX_k[a] = u_k[a] - trajectory(i, j, u_k, u_k[a], v_k[a]);

			f[a] = uX_k[a] * Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau - (P(Sigma_k[(i + 1) * M + j], e_k[(i + 1) * M + j]) - P(Sigma_k[(i - 1) * M + j], e_k[(i - 1) * M + j])) / (2 * hx);
		}
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			//uX_k[a] = u_k[a] - trajectory(i, j, u_k, u_k[a], v_k[a]);

			f[a] = uX_k[a] * Sigma_k1[i * M + j] * Sigma_k1[i * M + j] / tau - (P(Sigma_k[(i + 1) * M + j], e_k[(i + 1) * M + j]) - P(Sigma_k[(i - 1) * M + j], e_k[(i - 1) * M + j])) / (2 * hx);
		}
	}


	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;



		/*	f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[0*M+100],e_k[0*M+100]) - P(Sigma_k[(i-1)*M+j],e_k[(i-1)*M+j]))/(2*hx);

		}*/

	//Для Г2. l = 1,...,n-1; m = 0;

	/*j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;



		/*	f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[(i+1)*M+j],e_k[(i+1)*M+j]) - P(Sigma_k[(i-1)*M+j],e_k[(i-1)*M+j]))/(4*hx);

		}*/

	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			uX_k[a] = u_k[a] - trajectory(i, j, u_k, u_k[a], v_k[a]);

			f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[(i+1)*M+j],e_k[(i+1)*M+j]) - P(Sigma_k[(i)*M+j],e_k[(i)*M+j]))/(2*h);
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;

	/*j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

        /*    f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[(i+1)*M+j],e_k[(i+1)*M+j]) - P(Sigma_k[(i-1)*M+j],e_k[(i-1)*M+j]))/(4*hx);
		}*/

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

//			f[a] = u_k[a] * Ro(i,j,d)/(4*tau) - (P(i+1,j,d-1,e_k[(i+1)*M+j]) - P(i,j,d-1,e_k[i*M+j]))/(4*h);

//Для S_0n.
/*		i = 0;
		j = N;

			a = i*M + j;

			uX_k[a] = u_k[a] - trajectory(i, j, u_k, u_k[a], v_k[a]);

			f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) - (P(Sigma_k[(i+1)*M+j],e_k[(i+1)*M+j]) - P(Sigma_k[i*M+j],e_k[i*M+j]))/(4*h);
*/
	//Для S_nn.
	/*i = N1;
		j = N;

			a = i*M + j;

		//	f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) - (P(Sigma_k[0*M+100],e_k[0*M+100]) - P(Sigma_k[(i-1)*M+j],e_k[(i-1)*M+j]))/(4*hx);

//Для S_n0.
		/*i = N1;
		j = 0;

			a = i*M + j;

		//	f[a] = uX_k[a] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) - (P(Sigma_k[0*M+100],e_k[0*M+100]) - P(Sigma_k[(i-1)*M+j],e_k[(i-1)*M+j]))/(4*hx);

///////////////////////////////////////////////////Уравнение для v

//Для внутренних узлов. l,m = 1,...,n-1
    for(i = 1; i < qq; i++)
	{
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			//vY_k[i*M + j] = v_k[i*M + j] - trajectory(i, j, v_k, u_k[i*M + j], v_k[i*M + j]);

			f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/tau - (P(Sigma_k[i*M+j+1],e_k[i*M+j+1]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(2*hy);
		}
	}

    for(i = qq; i < qq+w; i++)
	{
		for(j = cntr+i+1-qq ; j < (M-1); j++)
		{

			a = M2 + i*M + j;

			//vY_k[i*M + j] = v_k[i*M + j] - trajectory(i, j, v_k, u_k[i*M + j], v_k[i*M + j]);

			f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/tau - (P(Sigma_k[i*M+j+1],e_k[i*M+j+1]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(2*hy);
		}
	}

	for(i = qq; i < qq+w; i++)
	{
		for(j = cntr-i-1+qq ; j > 0; j--)
		{

            a = M2 + i*M + j;

			//vY_k[i*M + j] = v_k[i*M + j] - trajectory(i, j, v_k, u_k[i*M + j], v_k[i*M + j]);

			f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/tau - (P(Sigma_k[i*M+j+1],e_k[i*M+j+1]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(2*hy);

		}
	}


	for(i = qq+w; i < (M1-1); i++)
	{
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			//vY_k[i*M + j] = v_k[i*M + j] - trajectory(i, j, v_k, u_k[i*M + j], v_k[i*M + j]);

			f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/tau - (P(Sigma_k[i*M+j+1],e_k[i*M+j+1]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(2*hy);
		}
	}


//Для Г1. l = n; m = 1,...,n-1;
		/*i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;



		/*	f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[0*M+100],e_k[0*M+100]) - P(Sigma_k[0*M+100],e_k[0*M+100]))/(4*hy);

		}*/

	//Для Г2. l = 1,...,n-1; m = 0;
	/*j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;



		/*	f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[i*M+j+1],e_k[i*M+j+1]) - P(Sigma_k[i*M+j],e_k[i*M+j]))/(2*hy);
		}*/

	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			vY_k[i*M + j] = v_k[i*M + j] - trajectory(i, j, v_k, u_k[i*M + j], v_k[i*M + j]);


			f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[i*M+j+1],e_k[i*M+j+1]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(4*h);
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;



		/*	f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(2*tau) - (P(Sigma_k[i*M+j],e_k[i*M+j]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(2*hy);
		}*/

	//Для S_00.
	/*		i = 0;
		j = 0;

			a = M2 + i*M + j;



//			f[a] = v_k[i*M + j] * Ro(i,j,d)/(4*tau) - (P(i,j+1,d-1,e_k[i*M+j+1]) - P(i,j,d-1,e_k[i*M+j]))/(4*h);

//Для S_0n.
/*		i = 0;
		j = N;

			a = M2 + i*M + j;


			vY_k[i*M + j] = v_k[i*M + j] - trajectory(i, j, v_k, u_k[i*M + j], v_k[i*M + j]);

			f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) - (P(Sigma_k[i*M+j],e_k[i*M+j]) - P(Sigma_k[i*M+j-1],e_k[i*M+j-1]))/(4*h);
*/
	//Для S_nn.
	/*i = N1;
		j = N;

			a = M2 + i*M + j;



		//	f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) - (P(Sigma_k[0*M+100],e_k[0*M+100]) - P(Sigma_k[0*M+100],e_k[0*M+100]))/(4*hy);

//Для S_n0.
		/*i = N1;
		j = 0;

			a = M2 + i*M + j;

			/*if(s_itr < itr_tr)
			{

			if((u_k[i*M + j] < 0) && (v_k[i*M + j] > 0))
            {
                vY_k[i*M + j] = v_kk[(i)*M + j]; //- ( u_k[i*M + j]*tau*( (V(i+1,j,d-1)-v_k[i*M+j])/((i+1)*h - i*h) )
                                //+ v_k[i*M + j]*tau*( (v_k[i*M+j-1]-v_k[i*M+j])/((j-1)*h - j*h) ) );
            }

            else if(v_k[i*M + j] > 0)
            {
                vY_k[i*M + j] = v_kk[i*M + j] - trajectory(i, j, v_kk, u_k[i*M + j], 0);
            }
            else if(u_k[i*M + j] < 0)
            {
                vY_k[i*M + j] = v_kk[i*M + j] - trajectory(i, j, v_kk, 0, v_k[i*M + j]);
            }

            else
                vY_k[i*M + j] = v_kk[i*M + j] - trajectory(i, j, v_kk, u_k[i*M + j], v_k[i*M + j]);
			}*/

	//	f[a] = vY_k[i*M + j] * Sigma_k1[i*M+j]*Sigma_k1[i*M+j]/(4*tau) - (P(Sigma_k[0*M+100],e_k[0*M+100]) - P(Sigma_k[0*M+100],e_k[0*M+100]))/(4*hy);


	return 0;
}


//Обратная диагональная матрица для матрицы А. Представлена в виде вектора из элементов обратных элементам главной диагонали матрицы А
inline double motion_D()
{
	int i = 0;
	int j = 0;
	int a;

	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;
			D[a] = 1 / A[a][2];
		}
	}


	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = M2 + i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = M2 + i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;
			D[a] = 1 / A[a][2];
		}
	}

	return 0;
}

//Вектор B = A*Xk1
inline double motion_B(double* u_k1, double* v_k1)
{
	int i = 0;
	int j = 0;
	int a;

	///////////////////////////////////////////////////Уравнение для u

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];
		}
	}


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr + i + 1 - qq;

		a = i * M + j;

		B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
			A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
			A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
			A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
			//A[a][7]*v_k1[(i+1)*M+(j-1)] +
			A[a][8] * v_k1[(i + 1) * M + (j + 1)] +
			A[a][9] * v_k1[(i) * M + (j - 1)] +
			A[a][11] * v_k1[(i + 1) * M + (j)];
	}

	i = qq + w - 1;
	j = cntr + i + 1 - qq;

	a = i * M + j;

	B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
		A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
		A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
		A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
		A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
		A[a][8] * v_k1[(i + 1) * M + (j + 1)];


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr - i - 1 + qq;

		a = i * M + j;

		B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
			A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
			A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
			A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
			A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
			A[a][10] * v_k1[(i) * M + (j + 1)] +
			A[a][11] * v_k1[(i + 1) * M + (j)];
		//A[a][8]*v_k1[(i+1)*M+(j+1)];
	}

	i = qq + w - 1;
	j = cntr - i - 1 + qq;

	a = i * M + j;

	B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
		A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
		A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
		A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
		A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
		A[a][8] * v_k1[(i + 1) * M + (j + 1)];


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 2 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 2 + qq; j > 0; j--)
		{
			a = i * M + j;

			B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			B[a] = A[a][0] * u_k1[(i - 1) * M + j] + A[a][1] * u_k1[i * M + (j - 1)] + A[a][2] * u_k1[i * M + j] +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];
		}
	}

	//Для Г1. l = n; m = 1,...,n-1;
	/*		i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] +
					A[a][3]*u_k1[i*M+j+1] +
					A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[i*M+j-1] +
					A[a][8]*v_k1[i*M+j+1];
		}
*/
	//Для Г2. l = 1,...,n-1; m = 0;
	/*		j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] +
					A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[(i-1)*M+j] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[(i+1)*M+j] +
					A[a][8]*v_k1[(i+1)*M+j+1];
		}
*/
	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			B[a] = A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] +
					A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[i*M+j-1] +
					A[a][6]*v_k1[i*M+j+1] +
					A[a][7]*v_k1[(i+1)*M+j-1] +
					A[a][8]*v_k1[(i+1)*M+j+1];
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*		j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] +
					A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j] +
					A[a][7]*v_k1[(i+1)*M+j-1] +
					A[a][8]*v_k1[(i+1)*M+j];
		}
*/
	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			B[a] = A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] + A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[i*M+j] +
					A[a][6]*v_k1[i*M+j+1] +
					A[a][7]*v_k1[(i+1)*M+j] +
					A[a][8]*v_k1[(i+1)*M+j+1];
*/

	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			B[a] = A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] + A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[i*M+j-1] +
					A[a][6]*v_k1[i*M+j] +
					A[a][7]*v_k1[(i+1)*M+j-1] +
					A[a][8]*v_k1[(i+1)*M+j];

*/
	//Для S_nn.
	/*		i = N1;
		j = N;

			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] +
					A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j] +
					A[a][7]*v_k1[i*M+j-1] +
					A[a][8]*v_k1[i*M+j];
*/

	//Для S_n0.
	/*		i = N1;
		j = 0;

			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] +
					A[a][5]*v_k1[(i-1)*M+j] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[i*M+j] +
					A[a][8]*v_k1[i*M+j+1];
*/
	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u_k1[(i - 1) * M + j - 1] +
				A[a][6] * u_k1[(i - 1) * M + j + 1] +
				A[a][7] * u_k1[(i + 1) * M + j - 1] +
				A[a][8] * u_k1[(i + 1) * M + j + 1];
		}
	}

	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr + i + 1 - qq;

		a = M2 + i * M + j;

		B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
			A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
			A[a][5] * u_k1[(i - 1) * M + j - 1] +
			A[a][6] * u_k1[(i - 1) * M + j + 1] +
			//A[a][7]*u_k1[(i+1)*M+j-1] +
			A[a][8] * u_k1[(i + 1) * M + j + 1] +
			A[a][9] * u_k1[(i) * M + j - 1] +
			A[a][11] * u_k1[(i + 1) * M + j];
	}

	i = qq + w - 1;
	j = cntr + i + 1 - qq;

	a = M2 + i * M + j;

	B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
		A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
		A[a][5] * u_k1[(i - 1) * M + j - 1] +
		A[a][6] * u_k1[(i - 1) * M + j + 1] +
		A[a][7] * u_k1[(i + 1) * M + j - 1] +
		A[a][8] * u_k1[(i + 1) * M + j + 1];


	for (i = qq; i < qq + w - 1; i++)
	{
		j = cntr - i - 1 + qq;

		a = M2 + i * M + j;

		B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
			A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
			A[a][5] * u_k1[(i - 1) * M + j - 1] +
			A[a][6] * u_k1[(i - 1) * M + j + 1] +
			A[a][7] * u_k1[(i + 1) * M + j - 1] +
			A[a][10] * u_k1[(i) * M + j + 1] +
			A[a][11] * u_k1[(i + 1) * M + j];
		//A[a][8]*u_k1[(i+1)*M+j+1];
	}

	i = qq + w - 1;
	j = cntr - i - 1 + qq;

	a = M2 + i * M + j;

	B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
		A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
		A[a][5] * u_k1[(i - 1) * M + j - 1] +
		A[a][6] * u_k1[(i - 1) * M + j + 1] +
		A[a][7] * u_k1[(i + 1) * M + j - 1] +
		A[a][8] * u_k1[(i + 1) * M + j + 1];


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 2 - qq; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u_k1[(i - 1) * M + j - 1] +
				A[a][6] * u_k1[(i - 1) * M + j + 1] +
				A[a][7] * u_k1[(i + 1) * M + j - 1] +
				A[a][8] * u_k1[(i + 1) * M + j + 1];
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 2 + qq; j > 0; j--)
		{
			a = M2 + i * M + j;

			B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u_k1[(i - 1) * M + j - 1] +
				A[a][6] * u_k1[(i - 1) * M + j + 1] +
				A[a][7] * u_k1[(i + 1) * M + j - 1] +
				A[a][8] * u_k1[(i + 1) * M + j + 1];
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			B[a] = A[a][0] * v_k1[(i - 1) * M + j] + A[a][1] * v_k1[i * M + j - 1] + A[a][2] * v_k1[i * M + j] +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u_k1[(i - 1) * M + j - 1] +
				A[a][6] * u_k1[(i - 1) * M + j + 1] +
				A[a][7] * u_k1[(i + 1) * M + j - 1] +
				A[a][8] * u_k1[(i + 1) * M + j + 1];
		}
	}

	//Для Г1. l = n; m = 1,...,n-1;
	/*		i = N1;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] +
					A[a][3]*v_k1[i*M+j+1] +
					A[a][5]*u_k1[(i-1)*M+j-1] +
					A[a][6]*u_k1[(i-1)*M+j+1] +
					A[a][7]*u_k1[i*M+j-1] +
					A[a][8]*u_k1[i*M+j+1];
		}
*/
	//Для Г2. l = 1,...,n-1; m = 0;
	/*		j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] +
					A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[(i-1)*M+j] +
					A[a][6]*u_k1[(i-1)*M+j+1] +
					A[a][7]*u_k1[(i+1)*M+j] +
					A[a][8]*u_k1[(i+1)*M+j+1];
		}
*/
	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] +
					A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[i*M+j-1] +
					A[a][6]*u_k1[i*M+j+1] +
					A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8]*u_k1[(i+1)*M+j+1];
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*		j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] +
					A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[(i-1)*M+j-1] +
					A[a][6]*u_k1[(i-1)*M+j] +
					A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8]*u_k1[(i+1)*M+j];
		}
*/
	//Для S_00.
	/*		i = 0;
		j = 0;

			a = M2 + i*M + j;

			B[a] = A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] + A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[i*M+j] +
					A[a][6]*u_k1[i*M+j+1] +
					A[a][7]*u_k1[(i+1)*M+j] +
					A[a][8]*u_k1[(i+1)*M+j+1];
*/

	//Для S_0n.
	/*		i = 0;
		j = N;

			a = M2 + i*M + j;

			B[a] = A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] + A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[i*M+j-1] +
					A[a][6]*u_k1[i*M+j] +
					A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8]*u_k1[(i+1)*M+j];

*/
	//Для S_nn.
	/*		i = N1;
		j = N;

			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] +
					A[a][5]*u_k1[(i-1)*M+j-1] +
					A[a][6]*u_k1[(i-1)*M+j] +
					A[a][7]*u_k1[i*M+j-1] +
					A[a][8]*u_k1[i*M+j];
*/

	//Для S_n0.
	/*		i = N1;
		j = 0;

			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] +
					A[a][5]*u_k1[(i-1)*M+j] +
					A[a][6]*u_k1[(i-1)*M+j+1] +
					A[a][7]*u_k1[i*M+j] +
					A[a][8]*u_k1[i*M+j+1];
*/

	return 0;
}

//Метод Якоби
inline double motion_Jakobi(double* u2, double* u_k1, double* v2, double* v_k1)
{
	int i = 0;
	int j = 0;
	int a;

	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			u2[i * M + j] = u_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			u2[i * M + j] = u_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;

			u2[i * M + j] = u_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			u2[i * M + j] = u_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}


	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			v2[i * M + j] = v_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			v2[i * M + j] = v_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = M2 + i * M + j;

			v2[i * M + j] = v_k1[i * M + j] - D[a] * (B[a] - f[a]);
		}
	}

	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			v2[i * M + j] = v_k1[i * M + j] - D[a] * (B[a] - f[a]);
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
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] /*+ A[a][2]*u_k1[i*M+j]*/ + A[a][3] * u_k1[i * M + j + 1] +
					A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + j] +
					A[a][6] * v_k1[(i - 1) * M + j + 1] +
					A[a][7] * v_k1[(i + 1) * M + j] +
					A[a][8] * v_k1[(i + 1) * M + j + 1];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < N))
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * M + (j + 1)];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + j - 1] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + j - 1] +
					A[a][6] * v_k1[(i - 1) * M + j] +
					A[a][7] * v_k1[(i + 1) * M + j - 1] +
					A[a][8] * v_k1[(i + 1) * M + j];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	for (i = qq; i < qq + w - 1; i++)
	{
		//for(j = cntr-i-2+qq ; j > 0; j--)
		for (j = 1; j <= cntr - i - 1 + qq; j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] /*+ A[a][2]*u_k1[i*M+j]*/ + A[a][3] * u_k1[i * M + j + 1] +
					A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + j] +
					A[a][6] * v_k1[(i - 1) * M + j + 1] +
					A[a][7] * v_k1[(i + 1) * M + j] +
					A[a][8] * v_k1[(i + 1) * M + j + 1];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < cntr - i - 1 + qq))
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * M + (j + 1)];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == cntr - i - 1 + qq)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
					A[a][10] * v_k1[(i) * M + (j + 1)] +
					A[a][11] * v_k1[(i + 1) * M + (j)];
				//A[a][8]*v_k1[(i+1)*M+(j+1)];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = qq + w - 1;
	for (j = 1; j <= cntr - i - 1 + qq; j++)
	{
		a = i * M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * u2[(i - 1) * M + j] /*+ A[a][2]*u_k1[i*M+j]*/ + A[a][3] * u_k1[i * M + j + 1] +
				A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + j] +
				A[a][6] * v_k1[(i - 1) * M + j + 1] +
				A[a][7] * v_k1[(i + 1) * M + j] +
				A[a][8] * v_k1[(i + 1) * M + j + 1];

			u2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > 0) && (j < cntr - i - 1 + qq))
		{
			B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];

			u2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == cntr - i - 1 + qq)
		{
			B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];

			u2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq; i < qq + w - 1; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			if (j == cntr + i + 1 - qq)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
					//A[a][7]*v_k1[(i+1)*M+(j-1)] +
					A[a][8] * v_k1[(i + 1) * M + (j + 1)] +
					A[a][9] * v_k1[(i) * M + (j - 1)] +
					A[a][11] * v_k1[(i + 1) * M + (j)];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > cntr + i + 1 - qq) && (j < N))
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * M + (j + 1)];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + j - 1] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + j - 1] +
					A[a][6] * v_k1[(i - 1) * M + j] +
					A[a][7] * v_k1[(i + 1) * M + j - 1] +
					A[a][8] * v_k1[(i + 1) * M + j];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = qq + w - 1;
	for (j = cntr + i + 1 - qq; j < (M - 1); j++)
	{
		a = i * M + j;

		if (j == cntr + i + 1 - qq)
		{
			B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];

			u2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > cntr + i + 1 - qq) && (j < N))
		{
			B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
				A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
				A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
				A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
				A[a][8] * v_k1[(i + 1) * M + (j + 1)];

			u2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == N)
		{
			B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + j - 1] /*+ A[a][2]*u_k1[i*M+j]*/ +
				A[a][4] * u_k1[(i + 1) * M + j] +
				A[a][5] * v_k1[(i - 1) * M + j - 1] +
				A[a][6] * v_k1[(i - 1) * M + j] +
				A[a][7] * v_k1[(i + 1) * M + j - 1] +
				A[a][8] * v_k1[(i + 1) * M + j];

			u2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] /*+ A[a][2]*u_k1[i*M+j]*/ + A[a][3] * u_k1[i * M + j + 1] +
					A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + j] +
					A[a][6] * v_k1[(i - 1) * M + j + 1] +
					A[a][7] * v_k1[(i + 1) * M + j] +
					A[a][8] * v_k1[(i + 1) * M + j + 1];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < N))
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + (j - 1)] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][3] * u_k1[i * M + (j + 1)] + A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + (j - 1)] +
					A[a][6] * v_k1[(i - 1) * M + (j + 1)] +
					A[a][7] * v_k1[(i + 1) * M + (j - 1)] +
					A[a][8] * v_k1[(i + 1) * M + (j + 1)];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * u2[(i - 1) * M + j] + A[a][1] * u2[i * M + j - 1] /*+ A[a][2]*u_k1[i*M+j]*/ +
					A[a][4] * u_k1[(i + 1) * M + j] +
					A[a][5] * v_k1[(i - 1) * M + j - 1] +
					A[a][6] * v_k1[(i - 1) * M + j] +
					A[a][7] * v_k1[(i + 1) * M + j - 1] +
					A[a][8] * v_k1[(i + 1) * M + j];

				u2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 0; j < (M); j++)
		{
			a = i*M + j;

			if(j == 0)
			{
			    B[a] = A[a][0]*u2[(i-1)*M+j] //+ A[a][2]*u_k1[i*M+j]
                    + A[a][3]*u_k1[i*M+j+1] +
					A[a][5]*v_k1[(i-1)*M+j] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[i*M+j] +
					A[a][8]*v_k1[i*M+j+1];

                u2[i*M + j] = D[a]*(f[a] - B[a]);
			}

            if((j > 0) && (j < N))
			{
                B[a] = A[a][0]*u2[(i-1)*M+j] + A[a][1]*u2[i*M+j-1]
                //+ A[a][2]*u_k1[i*M+j]
                + A[a][3]*u_k1[i*M+j+1] +
					A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[i*M+j-1] +
					A[a][8]*v_k1[i*M+j+1];

                u2[i*M + j] = D[a]*(f[a] - B[a]);
			}

			if(j == N)
			{
			    B[a] = A[a][0]*u2[(i-1)*M+j] + A[a][1]*u2[i*M+j-1] //+ A[a][2]*u_k1[i*M+j]
			    + A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j] +
					A[a][7]*v_k1[i*M+j-1] +
					A[a][8]*v_k1[i*M+j];

                u2[i*M + j] = D[a]*(f[a] - B[a]);
			}
		}*/

	//Для Г2. l = 1,...,n-1; m = 0;
	/*		j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] +
					A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[(i-1)*M+j] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[(i+1)*M+j] +
					A[a][8]*v_k1[(i+1)*M+j+1];
		}
*/
	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = i*M + j;

			B[a] = A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] +
					A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[i*M+j-1] +
					A[a][6]*v_k1[i*M+j+1] +
					A[a][7]*v_k1[(i+1)*M+j-1] +
					A[a][8]*v_k1[(i+1)*M+j+1];
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*		j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] +
					A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j] +
					A[a][7]*v_k1[(i+1)*M+j-1] +
					A[a][8]*v_k1[(i+1)*M+j];
		}
*/
	//Для S_00.
	/*		i = 0;
		j = 0;

			a = i*M + j;

			B[a] = A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] + A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[i*M+j] +
					A[a][6]*v_k1[i*M+j+1] +
					A[a][7]*v_k1[(i+1)*M+j] +
					A[a][8]*v_k1[(i+1)*M+j+1];
*/

	//Для S_0n.
	/*		i = 0;
		j = N;

			a = i*M + j;

			B[a] = A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] + A[a][4]*u_k1[(i+1)*M+j] +
					A[a][5]*v_k1[i*M+j-1] +
					A[a][6]*v_k1[i*M+j] +
					A[a][7]*v_k1[(i+1)*M+j-1] +
					A[a][8]*v_k1[(i+1)*M+j];

*/
	//Для S_nn.
	/*		i = N1;
		j = N;

			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][1]*u_k1[i*M+j-1] + A[a][2]*u_k1[i*M+j] +
					A[a][5]*v_k1[(i-1)*M+j-1] +
					A[a][6]*v_k1[(i-1)*M+j] +
					A[a][7]*v_k1[i*M+j-1] +
					A[a][8]*v_k1[i*M+j];

*/
	//Для S_n0.
	/*		i = N1;
		j = 0;

			a = i*M + j;

			B[a] = A[a][0]*u_k1[(i-1)*M+j] + A[a][2]*u_k1[i*M+j] + A[a][3]*u_k1[i*M+j+1] +
					A[a][5]*v_k1[(i-1)*M+j] +
					A[a][6]*v_k1[(i-1)*M+j+1] +
					A[a][7]*v_k1[i*M+j] +
					A[a][8]*v_k1[i*M+j+1];
*/
	///////////////////////////////////////////////////Уравнение для v

	//Для внутренних узлов. l,m = 1,...,n-1
	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] /*+ A[a][2]*v_k1[i*M+j]*/ + A[a][3] * v_k1[i * M + j + 1] +
					A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < N))
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}


			if (j == N)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	for (i = qq; i < qq + w - 1; i++)
	{
		//for(j = cntr-i-2+qq ; j > 0; j--)
		for (j = 1; j <= cntr - i - 1 + qq; j++)
		{
			a = M2 + i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] /*+ A[a][2]*v_k1[i*M+j]*/ + A[a][3] * v_k1[i * M + j + 1] +
					A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < cntr - i - 1 + qq))
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == cntr - i - 1 + qq)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][10] * u2[(i) * M + j + 1] +
					A[a][11] * u2[(i + 1) * M + j];
				//A[a][8]*u_k1[(i+1)*M+j+1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = qq + w - 1;
	for (j = 1; j <= cntr - i - 1 + qq; j++)
	{
		a = M2 + i * M + j;

		if (j == 0)
		{
			B[a] = A[a][0] * v2[(i - 1) * M + j] /*+ A[a][2]*v_k1[i*M+j]*/ + A[a][3] * v_k1[i * M + j + 1] +
				A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u2[(i - 1) * M + j] +
				A[a][6] * u2[(i - 1) * M + j + 1] +
				A[a][7] * u2[(i + 1) * M + j] +
				A[a][8] * u2[(i + 1) * M + j + 1];

			v2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > 0) && (j < cntr - i - 1 + qq))
		{
			B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u2[(i - 1) * M + j - 1] +
				A[a][6] * u2[(i - 1) * M + j + 1] +
				A[a][7] * u2[(i + 1) * M + j - 1] +
				A[a][8] * u2[(i + 1) * M + j + 1];

			v2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == cntr - i - 1 + qq)
		{
			B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u2[(i - 1) * M + j - 1] +
				A[a][6] * u2[(i - 1) * M + j + 1] +
				A[a][7] * u2[(i + 1) * M + j - 1] +
				A[a][8] * u2[(i + 1) * M + j + 1];

			v2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq; i < qq + w - 1; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			if (j == cntr + i + 1 - qq)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					//A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8] * u2[(i + 1) * M + j + 1] +
					A[a][9] * u2[(i) * M + j - 1] +
					A[a][11] * u2[(i + 1) * M + j];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > cntr + i + 1 - qq) && (j < N))
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}


	i = qq + w - 1;
	for (j = cntr + i + 1 - qq; j < (M - 1); j++)
	{
		a = M2 + i * M + j;

		if (j == cntr + i + 1 - qq)
		{
			B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u2[(i - 1) * M + j - 1] +
				A[a][6] * u2[(i - 1) * M + j + 1] +
				A[a][7] * u2[(i + 1) * M + j - 1] +
				A[a][8] * u2[(i + 1) * M + j + 1];

			v2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if ((j > cntr + i + 1 - qq) && (j < N))
		{
			B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
				A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u2[(i - 1) * M + j - 1] +
				A[a][6] * u2[(i - 1) * M + j + 1] +
				A[a][7] * u2[(i + 1) * M + j - 1] +
				A[a][8] * u2[(i + 1) * M + j + 1];

			v2[i * M + j] = D[a] * (f[a] - B[a]);
		}

		if (j == N)
		{
			B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
				A[a][4] * v_k1[(i + 1) * M + j] +
				A[a][5] * u2[(i - 1) * M + j - 1] +
				A[a][6] * u2[(i - 1) * M + j] +
				A[a][7] * u2[(i + 1) * M + j - 1] +
				A[a][8] * u2[(i + 1) * M + j];

			v2[i * M + j] = D[a] * (f[a] - B[a]);
		}
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = M2 + i * M + j;

			if (j == 0)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] /*+ A[a][2]*v_k1[i*M+j]*/ + A[a][3] * v_k1[i * M + j + 1] +
					A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if ((j > 0) && (j < N))
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][3] * v_k1[i * M + j + 1] + A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j + 1] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j + 1];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}

			if (j == N)
			{
				B[a] = A[a][0] * v2[(i - 1) * M + j] + A[a][1] * v2[i * M + j - 1] /*+ A[a][2]*v_k1[i*M+j]*/ +
					A[a][4] * v_k1[(i + 1) * M + j] +
					A[a][5] * u2[(i - 1) * M + j - 1] +
					A[a][6] * u2[(i - 1) * M + j] +
					A[a][7] * u2[(i + 1) * M + j - 1] +
					A[a][8] * u2[(i + 1) * M + j];

				v2[i * M + j] = D[a] * (f[a] - B[a]);
			}
		}
	}

	//Для Г1. l = n; m = 1,...,n-1;
	/*i = N1;
		for(j = 0; j < (M); j++)
		{
			a = M2 + i*M + j;

			if(j == 0)
			{
			    B[a] = A[a][0]*v2[(i-1)*M+j] //+ A[a][2]*v_k1[i*M+j]
			     + A[a][3]*v_k1[i*M+j+1] +
					A[a][5]*u2[(i-1)*M+j] +
					A[a][6]*u2[(i-1)*M+j+1] +
					A[a][7]*u2[i*M+j] +
					A[a][8]*u2[i*M+j+1];

                v2[i*M + j] = D[a]*(f[a] - B[a]);
			}

            if((j > 0) && (j < N))
            {
                B[a] = A[a][0]*v2[(i-1)*M+j] + A[a][1]*v2[i*M+j-1] //+ A[a][2]*v_k1[i*M+j]
                + A[a][3]*v_k1[i*M+j+1] +
					A[a][5]*u2[(i-1)*M+j-1] +
					A[a][6]*u2[(i-1)*M+j+1] +
					A[a][7]*u2[i*M+j-1] +
					A[a][8]*u2[i*M+j+1];

                v2[i*M + j] = D[a]*(f[a] - B[a]);
            }

            if(j == N)
            {
                B[a] = A[a][0]*v2[(i-1)*M+j] + A[a][1]*v2[i*M+j-1] //+ A[a][2]*v_k1[i*M+j]
                + A[a][5]*u2[(i-1)*M+j-1] +
					A[a][6]*u2[(i-1)*M+j] +
					A[a][7]*u2[i*M+j-1] +
					A[a][8]*u2[i*M+j];

                v2[i*M + j] = D[a]*(f[a] - B[a]);
            }

		}*/

	//Для Г2. l = 1,...,n-1; m = 0;
	/*		j = 0;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] +
					A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[(i-1)*M+j] +
					A[a][6]*u_k1[(i-1)*M+j+1] +
					A[a][7]*u_k1[(i+1)*M+j] +
					A[a][8]*u_k1[(i+1)*M+j+1];
		}
*/
	//Для Г3. l = 0; m = 1,...,n-1;
	/*		i = 0;
		for(j = 1; j < (M-1); j++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] +
					A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[i*M+j-1] +
					A[a][6]*u_k1[i*M+j+1] +
					A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8]*u_k1[(i+1)*M+j+1];
		}
*/
	//Для Г4. l = 1,...,n-1; m = n;
	/*		j = N;
		for(i = 1; i < (M1-1); i++)
		{
			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] +
					A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[(i-1)*M+j-1] +
					A[a][6]*u_k1[(i-1)*M+j] +
					A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8]*u_k1[(i+1)*M+j];
		}
*/
	//Для S_00.
	/*		i = 0;
		j = 0;

			a = M2 + i*M + j;

			B[a] = A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] + A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[i*M+j] +
					A[a][6]*u_k1[i*M+j+1] +
					A[a][7]*u_k1[(i+1)*M+j] +
					A[a][8]*u_k1[(i+1)*M+j+1];
*/

	//Для S_0n.
	/*		i = 0;
		j = N;

			a = M2 + i*M + j;

			B[a] = A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] + A[a][4]*v_k1[(i+1)*M+j] +
					A[a][5]*u_k1[i*M+j-1] +
					A[a][6]*u_k1[i*M+j] +
					A[a][7]*u_k1[(i+1)*M+j-1] +
					A[a][8]*u_k1[(i+1)*M+j];

*/
	//Для S_nn.
	/*		i = N1;
		j = N;

			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][1]*v_k1[i*M+j-1] + A[a][2]*v_k1[i*M+j] +
					A[a][5]*u_k1[(i-1)*M+j-1] +
					A[a][6]*u_k1[(i-1)*M+j] +
					A[a][7]*u_k1[i*M+j-1] +
					A[a][8]*u_k1[i*M+j];

*/
	//Для S_n0.
	/*		i = N1;
		j = 0;

			a = M2 + i*M + j;

			B[a] = A[a][0]*v_k1[(i-1)*M+j] + A[a][2]*v_k1[i*M+j] + A[a][3]*v_k1[i*M+j+1] +
					A[a][5]*u_k1[(i-1)*M+j] +
					A[a][6]*u_k1[(i-1)*M+j+1] +
					A[a][7]*u_k1[i*M+j] +
					A[a][8]*u_k1[i*M+j+1];

*/

	/*fprintf(out,"B[ ] = \n",NULL);
	for(i=0;i<=N;i++)
		{
			fprintf(out,"%.5f ",B[i]);
		}
		fprintf(out,"\n\n ",NULL);*/
	return 0;
}


inline double motion(double* sigma_k1, double* sigma_k, double* u_k, double* v_k, double* u_k1, double* v_k1, double* u2, double* v2, double* e_k)
{
	int i = 0;
	int j = 0;
	int a;
	int bl = 1;
	int c_u;
	int c_v;

	/*---------------------------------------------*/

	motion_A(sigma_k1, e_k);
	motion_D();

	motion_f(sigma_k, sigma_k1, u_k, v_k, e_k);

	s_m = 0;
	//c = 0;
	while (bl)
	{
		motion_B(u_k1, v_k1);
		motion_Jakobi(u2, u_k1, v2, v_k1);

		//motion_Zeidel(u_k1, v_k1, u2, v2);

		c_u = 0;
		c_v = 0;

		//	if(s >= 1)
		//	{

		for (i = 1; i < qq; i++)
		{
			for (j = 1; j < (M - 1); j++)
			{
				a = i * M + j;
				if ((fabs(u_k1[a] - u2[a]) <= epsilon))
				{
					c_u += 1;
					//bl = 0;
				}

				if ((fabs(v_k1[a] - v2[a]) <= epsilon))
				{
					c_v += 1;
					//bl = 0;
				}
			}
		}

		for (i = qq; i < qq + w; i++)
		{
			for (j = cntr + i + 1 - qq; j < (M - 1); j++)
			{
				a = i * M + j;
				if ((fabs(u_k1[a] - u2[a]) <= epsilon))
				{
					c_u += 1;
					//bl = 0;
				}

				if ((fabs(v_k1[a] - v2[a]) <= epsilon))
				{
					c_v += 1;
					//bl = 0;
				}
			}
		}

		for (i = qq; i < qq + w; i++)
		{
			for (j = cntr - i - 1 + qq; j > 0; j--)
			{
				a = i * M + j;
				if ((fabs(u_k1[a] - u2[a]) <= epsilon))
				{
					c_u += 1;
					//bl = 0;
				}

				if ((fabs(v_k1[a] - v2[a]) <= epsilon))
				{
					c_v += 1;
					//bl = 0;
				}
			}
		}

		for (i = qq + w; i < (M1 - 1); i++)
		{
			for (j = 1; j < (M - 1); j++)
			{
				a = i * M + j;
				if ((fabs(u_k1[a] - u2[a]) <= epsilon))
				{
					c_u += 1;
					//bl = 0;
				}

				if ((fabs(v_k1[a] - v2[a]) <= epsilon))
				{
					c_v += 1;
					//bl = 0;
				}
			}
		}


		if ((c_u == (N1 - 1) * (N - 1) - (2 + (q - 1) * 2) / 2 * (q)) && (c_v >= (N1 - 1) * (N - 1) - (2 + (q - 1) * 2) / 2 * (q)))
		{
			bl = 0;
		}

		else if (s_m > 20)
		{
			//fprintf(out,"\nk = %i   t = %f   h = %f   s = %i\n\n",k,tau,h,s1);
			bl = 0;
		}

		else
		{
			for (i = 1; i < qq; i++)
			{
				for (j = 1; j < (M - 1); j++)
				{
					a = i * M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}

			for (i = qq; i < qq + w; i++)
			{
				for (j = cntr + i + 1 - qq; j < (M - 1); j++)
				{
					a = i * M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}

			for (i = qq; i < qq + w; i++)
			{
				for (j = cntr - i - 1 + qq; j > 0; j--)
				{
					a = i * M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}

			for (i = qq + w; i < (M1 - 1); i++)
			{
				for (j = 1; j < (M - 1); j++)
				{
					a = i * M + j;
					u_k1[a] = u2[a];
					v_k1[a] = v2[a];
				}
			}
		}

		s_m += 1;
		//if (s_m >= 1)
		//    bl = 0;
	}
	return 0;
}
