#ifndef __CALC_H__
#define __CALC_H__

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 1
#define omp_get_num_threads() 1
#endif

#ifdef __GNUC__
#define __pure
/*почему то на GCC не работает и ломает результаты расчета*/
//#define __pure __attribute__((const))
// for memory alignment
#define _ALIGN(N)  __attribute__((aligned(N))) 
#elif __INTEL_COMPILER
#define __pure __declspec(const)
// for memory alignment
#define _ALIGN(N)  __declspec(align(N))
#elif __NVCC__
#define __pure
#else
#define __pure 
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
typedef const float_type* __restrict const cnst_arr_t;
typedef float_type* __restrict const cnst_ptr_arr_t;
typedef float_type** __restrict const cnst_ptr_2d_arr_t;

inline __pure float_type Mu(float_type gamma_mah2, float_type value)
{
	const float_type omega = 0.8;
	return pow(gamma_mah2 * value * value, omega);
}

#pragma omp declare simd
inline __pure float_type P(float_type gamma_value, float_type sigma_k_value, float_type e_k_value)
{
	return (gamma_value - 1) * sigma_k_value * sigma_k_value * e_k_value * e_k_value;
}

#endif // __CALC_H__