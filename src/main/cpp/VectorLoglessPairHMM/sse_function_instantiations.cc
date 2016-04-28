#if defined(__POWER8_VECTOR__)
#define SIMD_ENGINE sse
#define SIMD_ENGINE_SSE

#include "template.h"

#include "define-sse-float.h"
#include "shift_template.cc"
#include "pairhmm-template-kernel.cc"

#include "define-sse-double.h"
#include "shift_template.cc"
#include "pairhmm-template-kernel.cc"

template double compute_full_prob_ssed<double>(testcase* tc, double* nextlog);
template float compute_full_prob_sses<float>(testcase* tc, float* nextlog);
#endif
