inline double continuity(double* Sigma_k, double* Sigma_k1, double* u_k, double* v_k)
{
	int i = 0, j = 0, a = 0;

	//Для внутренних узлов

	for (i = 1; i < qq; i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

			Sigma_k1[a] = (SigmaX_k[a] / tau) / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr + i + 1 - qq; j < (M - 1); j++)
		{
			a = i * M + j;

			//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

			Sigma_k1[a] = (SigmaX_k[a] / tau) / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}


	for (i = qq; i < qq + w; i++)
	{
		for (j = cntr - i - 1 + qq; j > 0; j--)
		{
			a = i * M + j;

			//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

			Sigma_k1[a] = (SigmaX_k[a] / tau) / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}


	for (i = qq + w; i < (M1 - 1); i++)
	{
		for (j = 1; j < (M - 1); j++)
		{
			a = i * M + j;

			//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

			Sigma_k1[a] = (SigmaX_k[a] / tau) / (1 / tau + (u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
				+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (4 * hy));
		}
	}

	//Для Г1. l = n; m = 1,...,n-1;


	/* Sigma_k1[a] = ( SigmaX_k[a]/(tau) )
			/( 1/(tau) + (u_k[i*M+j]-u_k[(i-1)*M+j])/(2*hx) + (v_k[i*M+j+1]-v_k[i*M+j-1])/(4*hy) );
	*/

	//Для Г2. l = w,...,n-1; m = 0;


	/*		Sigma_k1[a] = (SigmaX_k[a]/(2*tau))/( 1/(2*tau) + (u_k[(i+1)*M+j]-u_k[(i-1)*M+j])/(8*hx)
				+ (v_k[i*M+j+1]-v_k[i*M+j])/(4*hy) );

	*/


	//Для Г3. l = 0; m = q,...,n-1;
	/*	i = 0;
	for(j = 1; j < (M-1); j++)
	{

			a = i*M + j;

			//SigmaX_k[a] = Sigma_k[a] - (h*i*u_k[i*M+j] + h*j*v_k[i*M+j])/sqrt(j*j*h*h)*tau* (Sigma_k[i*M+j-1]-Sigma_k[a])
			//		/( (j-1)*h - j*h );

			SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

			Sigma_k1[a] = (SigmaX_k[a]/(2*tau))/( 1/(2*tau) + (u_k[(i+1)*M+j]-u_k[i*M+j])/(4*h)
				+ (v_k[i*M+j+1]-v_k[i*M+j-1])/(8*h) );

	}
*/

	//Для Г4. l = 1,...,n-1; m = n;


	/*		Sigma_k1[a] = (SigmaX_k[a]/(2*tau))/( 1/(2*tau) + (u_k[(i+1)*M+j]-u_k[(i-1)*M+j])/(8*hx)
				+ (v_k[i*M+j]-v_k[i*M+j-1])/(4*hy) );

	*/


	//Для Г5. l = w-1; m = 1,...,q-2;
	i = qq + w - 1;
	for (j = cntr - q + 2; j < (cntr + q - 1); j++)
	{
		a = i * M + j;

		//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

		Sigma_k1[a] = (SigmaX_k[a] / (2 * tau)) / (1 / (2 * tau) + (u_k[(i + 1) * M + j] - u_k[i * M + j]) / (4 * hx)
			+ (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (8 * hy));
	}

	//Для Г6.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr + i - qq;

		a = i * M + j;

		//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

		Sigma_k1[a] = (SigmaX_k[a] * (1 / (4 * tau) + 1 / (4 * tau))) / (1 / (4 * tau) + 1 / (4 * tau) + (u_k[(i) * M + j] - u_k[(i - 1) * M + j]) / (8 * hx)
			- u_k[(i - 1) * M + j] / (16 * hx) + (v_k[i * M + j + 1] - v_k[i * M + j]) / (8 * hy) + v_k[i * M + j + 1] / (16 * hy));
	}

	//Для Г7.
	for (i = qq + 1; i < qq + w - 1; i++)
	{
		j = cntr - i + qq;

		a = i * M + j;

		//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

		Sigma_k1[a] = (SigmaX_k[a] * (1 / (4 * tau) + 1 / (4 * tau))) / (1 / (4 * tau) + 1 / (4 * tau) + (u_k[(i) * M + j] - u_k[(i - 1) * M + j]) / (8 * hx)
			- u_k[(i - 1) * M + j] / (16 * hx) + (v_k[i * M + j] - v_k[i * M + j - 1]) / (8 * hy) - v_k[i * M + j - 1] / (16 * hy));
	}


	//Для S_00.
	/*	i = 0;
	j = 0;

		a = i*M + j;

		SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

		Sigma_k1[a] = (SigmaX_k[a]/(4*tau))/( 1/(4*tau) + (u_k[(i+1)*M+j]-u_k[i*M+j])/(8*h)
			+ (v_k[i*M+j+1]-v_k[i*M+j])/(8*h) );

*/
	//Для S_0q-1.
	/*	i = 0;
	j = q-1;

		a = i*M + j;

		SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

		Sigma_k1[a] = (SigmaX_k[a]/(4*tau))/( 1/(4*tau) + (u_k[(i+1)*M+j]-u_k[i*M+j])/(8*h)
			+ (v_k[i*M+j+1]-v_k[i*M+j])/(8*h) ); */

	//Для S_0n.
	/*	i = 0;
	j = N;

		a = i*M + j;

		//SigmaX_k[a] = Sigma_k[a] - (h*i*u_k[i*M+j] + h*j*v_k[i*M+j])/sqrt(i*i*h*h + j*j*h*h)*tau* (Sigma_k[i*M+j-1]-Sigma_k[a])
		//			/( sqrt((j-1)*(j-1)*h*h) - sqrt(j*j*h*h) );

		SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

		Sigma_k1[a] = (SigmaX_k[a]/(4*tau))/( 1/(4*tau) + (u_k[(i+1)*M+j]-u_k[i*M+j])/(8*h)
			+ (v_k[i*M+j]-v_k[i*M+j-1])/(8*h) );

*/
	//Для S_nn.
	/*i = N1;
	j = N;

		a = i*M + j;



		/*Sigma_k1[a] = (SigmaX_k[a]/(4*tau))/( 1/(4*tau) + (u_k[i*M+j]-u_k[(i-1)*M+j])/(8*hx)
			+ (v_k[i*M+j]-v_k[i*M+j-1])/(8*hy) );*/


	//Для S_n0.
	/*i = N1;
	j = 0;

		a = i*M + j;



		/*Sigma_k1[a] = (SigmaX_k[a]/(4*tau))/( 1/(4*tau) + (u_k[i*M+j]-u_k[(i-1)*M+j])/(8*hx)
			+ (v_k[i*M+j+1]-v_k[i*M+j])/(8*hy) );*/


	//Для S_w-1q-1.
	i = qq + w - 1;
	j = cntr + i - qq;

	a = i * M + j;

	//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

	Sigma_k1[a] = (SigmaX_k[a] * (3 / (4 * tau) + 1 / (8 * tau))) / (3 / (4 * tau) + 1 / (8 * tau) + (2 * u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j] - u_k[i * M + j]) / (8 * hx)
		+ (u_k[(i) * M + j] - u_k[(i - 1) * M + j]) / (16 * hx) + (2 * v_k[i * M + j + 1] - v_k[i * M + j - 1] - v_k[i * M + j]) / (8 * hy) + v_k[i * M + j] / (16 * hy));

	//Для S.
	i = qq + w - 1;
	j = cntr - i + qq;

	a = i * M + j;

	//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

	Sigma_k1[a] = (SigmaX_k[a] * (3 / (4 * tau) + 1 / (8 * tau))) / (3 / (4 * tau) + 1 / (8 * tau) + (2 * u_k[(i + 1) * M + j] - u_k[(i - 1) * M + j] - u_k[i * M + j]) / (8 * hx)
		+ (u_k[(i) * M + j] - u_k[(i - 1) * M + j]) / (16 * hx) + (v_k[i * M + j + 1] - 2 * v_k[i * M + j - 1] + v_k[i * M + j]) / (8 * hy) - v_k[i * M + j] / (16 * hy));


	//Для S_qq_0
	i = qq;
	j = cntr;

	a = i * M + j;

	//SigmaX_k[a] = Sigma_k[a] - trajectory(i, j, Sigma_k, u_k[a], v_k[a]);

	Sigma_k1[a] = (SigmaX_k[a] * (1 / (4 * tau) + 1 / (2 * tau))) / ((1 / (4 * tau) + 1 / (2 * tau)) + (u_k[(i) * M + j] - u_k[(i - 1) * M + j]) / (4 * hx)
		+ (u_k[(i + 1) * M + j] - u_k[i * M + j]) / (8 * hx) + (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (8 * hy) + (v_k[i * M + j + 1] - v_k[i * M + j - 1]) / (16 * hy));


	return 0;
}
