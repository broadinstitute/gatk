#if defined(__x86_64__)
#define SIMD_ENGINE avx
#define SIMD_ENGINE_AVX

#include "template.h"

#include "define-float.h"
#include "shift_template.cc"
#include "pairhmm-template-kernel.cc"

#include "define-double.h"
#include "shift_template.cc"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_avxd<double>(testcase* tc, double* nextlog);
template float compute_full_prob_avxs<float>(testcase* tc, float* nextlog);
#endif
