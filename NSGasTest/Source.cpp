//// m2 = M2
//inline void init_arrays(const int array_element_count, const int param_array_element_count)
//{
//	A = new double*[2 * array_element_count];
//	for (int i = 0; i < 2 * array_element_count; ++i)
//	{
//		A[i] = new double[param_array_element_count];
//	}
//	B = new double[2 * array_element_count];
//	D = new double[2 * array_element_count];
//	f = new double[2 * array_element_count];
//	sigma_k = new double[array_element_count];
//	sigma_k1 = new double[array_element_count];
//	u_k = new double[array_element_count];
//	u_k1 = new double[array_element_count];
//	v_k = new double[array_element_count];
//	v_k1 = new double[array_element_count];
//	u2 = new double[array_element_count];
//	v2 = new double[array_element_count];
//	e_k = new double[array_element_count];
//	e_k1 = new double[array_element_count];
//	e2 = new double[array_element_count];
//	T = new double[array_element_count];
//	sigma_kk = new double[array_element_count];
//	u_kk = new double[array_element_count];
//	v_kk = new double[array_element_count];
//	e_kk = new double[array_element_count];
//	sigmaX_k = new double[array_element_count];
//	uX_k = new double[array_element_count];
//	vY_k = new double[array_element_count];
//	eR_k = new double[array_element_count];
//	for (int i = 0; i < 2 * array_element_count; ++i)
//	{
//		std::fill_n(A[i], param_array_element_count, 0.);
//	}
//	std::fill_n(B, 2 * array_element_count, 0.);
//	std::fill_n(D, 2 * array_element_count, 0.);
//	std::fill_n(f, 2 * array_element_count, 0.);
//	std::fill_n(sigma_k, array_element_count, 0.);
//	std::fill_n(sigma_k1, array_element_count, 0.);
//	std::fill_n(u_k, array_element_count, 0.);
//	std::fill_n(v_k, array_element_count, 0.);
//	std::fill_n(u_k1, array_element_count, 0.);
//	std::fill_n(v_k1, array_element_count, 0.);
//	std::fill_n(u2, array_element_count, 0.);
//	std::fill_n(v2, array_element_count, 0.);
//	std::fill_n(e_k, array_element_count, 0.);
//	std::fill_n(e_k1, array_element_count, 0.);
//	std::fill_n(e2, array_element_count, 0.);
//	std::fill_n(T, array_element_count, 0.);
//	std::fill_n(sigma_kk, array_element_count, 0.);
//	std::fill_n(u_kk, array_element_count, 0.);
//	std::fill_n(v_kk, array_element_count, 0.);
//	std::fill_n(e_kk, array_element_count, 0.);
//}
//
//// m2 = M2
//inline void clear_memory(const int array_element_count)
//{
//	for (int i = 0; i < 2 * array_element_count; i++)
//	{
//		delete [] A[i];
//	}
//	delete [] A;
//
//	delete [] D;
//	delete [] B;
//	delete [] f;
//	delete [] sigma_k;
//	delete [] sigma_k1;
//	delete [] u_k;
//	delete [] u_k1;
//	delete [] v_k;
//	delete [] v_k1;
//	delete [] u2;
//	delete [] v2;
//	delete [] e_k;
//	delete [] e_k1;
//	delete [] e2;
//	delete [] T;
//	delete [] sigma_kk;
//	delete [] u_kk;
//	delete [] v_kk;
//	delete [] e_kk;
//	delete [] sigmaX_k;
//	delete [] uX_k;
//	delete [] vY_k;
//	delete [] eR_k;
//}

//
//FILE* fout = nullptr;
//FILE* fout_itr = nullptr;
//FILE* fdensity = nullptr;
//FILE* fdensity_new = nullptr;
//FILE* fvelocity = nullptr;
//FILE* ftemperature = nullptr;
//FILE* fpressure = nullptr;
//if (need_print)
//{
//	fout = fopen("out_p.txt", "w");
//	fout_itr = fopen("out_itr_p.txt", "w");
//	fdensity = fopen("density_p.dat", "w");
//	fdensity_new = fopen("density-new_p.dat", "w");
//	fvelocity = fopen("velocity_p.dat", "w");
//	ftemperature = fopen("temperature_p.dat", "w");
//	fpressure = fopen("pressure_p.dat", "w");
//	print_file_header(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr, tau, hx, hy, N);
//}
//const double gamma = 1.4;
//const int time_steps_nbr = 1; // time_steps_nbr - количество шагов по времени
//init_arrays(M2, 12);
//
//set_initial_boundary_conditions(gamma, qq, w, M, M1, M2, Mah2);
//StartTimer();
//for (int current_time_step = 1; current_time_step <= time_steps_nbr; current_time_step++)
//{
//	int s_end = 0;
//	int s_m = 0;
//	int s_e = 0;
//	int s_itr;
//	prepare_to_iterate(M, M1, qq, w, cntr);
//	s_itr = interate_over_nonlinearity(gamma, qq, M, M1, w, cntr, N, q, s_m, s_e, s_end);
//	if (need_print)
//		print_to_file(gamma, s_m, s_e, current_time_step, s_itr, s_end, tau, hx, hy, M, M1, N, Mah2, fout, fdensity, fdensity_new, fvelocity, ftemperature, fpressure, fout_itr);
//}
//double time = GetTimer();
//if (need_print)
//{
//	print_new_line(fout, fdensity, fvelocity, ftemperature, fpressure);
//	close_files(fout, fdensity, fvelocity, ftemperature, fpressure, fout_itr);
//}
//clear_memory(M2);
//return time;