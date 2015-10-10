inline double trajectory(int i, int j, double* arr, double u_k, double v_k)
{
	double result = 0;

	if (u_k == 0 && v_k == 0)
	{
		result = 0;
	}

	if (u_k > 0 && v_k > 0)
	{
		result = u_k * tau * ((arr[(i - 1) * M + j] - arr[i * M + j]) / ((i - 1) * hx - i * hx))
			+ v_k * tau * ((arr[i * M + (j - 1)] - arr[i * M + j]) / ((j - 1) * hy - j * hy));
	}


	if (u_k < 0 && v_k > 0)
	{
		result = u_k * tau * ((arr[(i + 1) * M + j] - arr[i * M + j]) / ((i + 1) * hx - i * hx))
			+ v_k * tau * ((arr[i * M + (j - 1)] - arr[i * M + j]) / ((j - 1) * hy - j * hy));
	}


	if (u_k < 0 && v_k < 0)
	{
		result = u_k * tau * ((arr[(i + 1) * M + j] - arr[i * M + j]) / ((i + 1) * hx - i * hx))
			+ v_k * tau * ((arr[i * M + (j + 1)] - arr[i * M + j]) / ((j + 1) * hy - j * hy));
	}

	if (u_k > 0 && v_k < 0)
	{
		result = u_k * tau * ((arr[(i - 1) * M + j] - arr[i * M + j]) / ((i - 1) * hx - i * hx))
			+ v_k * tau * ((arr[i * M + (j + 1)] - arr[i * M + j]) / ((j + 1) * hy - j * hy));
	}


	if (u_k > 0 && v_k == 0)
	{
		result = u_k * tau * ((arr[(i - 1) * M + j] - arr[i * M + j]) / ((i - 1) * hx - i * hx));
	}


	if (u_k == 0 && v_k > 0)
	{
		result = v_k * tau * ((arr[i * M + (j - 1)] - arr[i * M + j]) / ((j - 1) * hy - j * hy));
	}


	if (u_k < 0 && v_k == 0)
	{
		result = u_k * tau * ((arr[(i + 1) * M + j] - arr[i * M + j]) / ((i + 1) * hx - i * hx));
	}


	if (u_k == 0 && v_k < 0)
	{
		result = v_k * tau * ((arr[i * M + (j + 1)] - arr[i * M + j]) / ((j + 1) * hy - j * hy));
	}

	return result;
}
