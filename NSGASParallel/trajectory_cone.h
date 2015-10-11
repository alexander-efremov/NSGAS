inline double trajectory(int i, int j, double* arr, double u_k, double v_k, const int m)
{
	double result = 0;
	if (u_k == 0 && v_k == 0)
	{
		result = 0;
	}
	int idx = i * m + j;
	int idx2 = i * m + (j - 1);
	if (u_k > 0 && v_k > 0)
	{
		result = u_k * tau * ((arr[(i - 1) * m + j] - arr[idx]) / ((i - 1) * hx - i * hx))
			+ v_k * tau * ((arr[idx2] - arr[idx]) / ((j - 1) * hy - j * hy));
	}

	if (u_k < 0 && v_k > 0)
	{
		result = u_k * tau * ((arr[(i + 1) * m + j] - arr[idx]) / ((i + 1) * hx - i * hx))
			+ v_k * tau * ((arr[idx2] - arr[idx]) / ((j - 1) * hy - j * hy));
	}

	if (u_k < 0 && v_k < 0)
	{
		result = u_k * tau * ((arr[(i + 1) * m + j] - arr[idx]) / ((i + 1) * hx - i * hx))
			+ v_k * tau * ((arr[i * m + (j + 1)] - arr[idx]) / ((j + 1) * hy - j * hy));
	}

	if (u_k > 0 && v_k < 0)
	{
		result = u_k * tau * ((arr[(i - 1) * m + j] - arr[idx]) / ((i - 1) * hx - i * hx))
			+ v_k * tau * ((arr[i * m + (j + 1)] - arr[idx]) / ((j + 1) * hy - j * hy));
	}

	if (u_k > 0 && v_k == 0)
	{
		result = u_k * tau * ((arr[(i - 1) * m + j] - arr[idx]) / ((i - 1) * hx - i * hx));
	}

	if (u_k == 0 && v_k > 0)
	{
		result = v_k * tau * ((arr[idx2] - arr[idx]) / ((j - 1) * hy - j * hy));
	}

	if (u_k < 0 && v_k == 0)
	{
		result = u_k * tau * ((arr[(i + 1) * m + j] - arr[idx]) / ((i + 1) * hx - i * hx));
	}

	if (u_k == 0 && v_k < 0)
	{
		result = v_k * tau * ((arr[i * m + (j + 1)] - arr[idx]) / ((j + 1) * hy - j * hy));
	}
	return result;
}
