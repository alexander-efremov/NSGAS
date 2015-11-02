#ifndef __CALC_H__
#define __CALC_H__

#include <math.h>
#include <malloc.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 1
#define omp_get_num_threads() 1
#endif

#ifdef __MIC__
#define ALIGN 64
#elif __AVX__
#define ALIGN 32
#elif __AVX2__
#define ALIGN 32
#elif __SSE2__
#define ALIGN 16
#elif __SSE3__
#define ALIGN 16
#elif __SSE4__
#define ALIGN 16
#elif __SSE4_1__
#define ALIGN 16
#elif __SSE4_2__
#define ALIGN 16
#else
#define ALIGN 16
#endif

typedef double float_type;
// указывая __restrict мы подсказываем компилятору, что области памяти не пересекаются
// то есть нет явления memory aliasing
// тогда компилятор может применить автоматическую векторизацию цикла
// таким образом уничтожаются зависимости типа ANTI и FLOW
typedef const float_type* __restrict const cnst_arr_t;
typedef float_type* __restrict const cnst_ptr_arr_t;
typedef float_type** __restrict const cnst_ptr_2d_arr_t;

#pragma omp declare simd uniform(gamma_mah2) linear(value)
float_type Mu(float_type gamma_mah2, float_type value)
{	
	// 0.8 = omega
	return pow(gamma_mah2 * value * value, 0.8);
}

#pragma omp declare simd uniform(gamma_value) linear(sigma_k_value, e_k_value)
float_type P(float_type gamma_value, float_type sigma_k_value, float_type e_k_value)
{
	return (gamma_value - 1) * sigma_k_value * sigma_k_value * e_k_value * e_k_value;
}

#endif // __CALC_H__