#include "headers.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;

//static members from ConvertChar
uint8_t ConvertChar::conversionTable[255];
//Global function pointers in utils.h
float (*g_compute_full_prob_float)(testcase *tc, float* before_last_log) = 0;
double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log) = 0;
#if defined(__POWER8_VECTOR__)
unsigned long g_max_num_threads = omp_get_num_procs() *3/8; // 3 threads per core
#endif
//Static members in ContextBase
template<>
bool ContextBase<double>::staticMembersInitializedFlag = false;
template<>
double ContextBase<double>::jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE] = { };
template<>
double ContextBase<double>::matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] = { };
template<>
bool ContextBase<float>::staticMembersInitializedFlag = false;
template<>
float ContextBase<float>::jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE] = { };
template<>
float ContextBase<float>::matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] = { };
template<> bool ContextBase<double>::staticMembersInitializedFlag1 = false;
template<> bool ContextBase<float>::staticMembersInitializedFlag1 = false;
template<> double ContextBase<double>::ph2pr[128] = {0.0};
template<> float ContextBase<float>::ph2pr[128] = {0.0F};
template<> double ContextBase<double>::INITIAL_CONSTANT = 0.0;
template<> float ContextBase<float>::INITIAL_CONSTANT = 0.0F;
template<> double ContextBase<double>::LOG10_INITIAL_CONSTANT = 0.0;
template<> float ContextBase<float>::LOG10_INITIAL_CONSTANT = 0.0F;
template<> double ContextBase<double>::RESULT_THRESHOLD = 0.0;
template<> float ContextBase<float>::RESULT_THRESHOLD = 0.0F;


void initialize_function_pointers()
{
#if defined(__x86_64__)
  g_compute_full_prob_float = compute_full_prob_avxs<float>;
  g_compute_full_prob_double = compute_full_prob_avxd<double>;
#elif defined(__POWER8_VECTOR__)
  g_compute_full_prob_float = compute_full_prob_sses<float>;
  g_compute_full_prob_double = compute_full_prob_ssed<double>;
  char *s = getenv("PHMM_N_THREADS");
  if (s) {
    char *endp;
    long l = strtol(s, &endp, 10);
    if (endp && *endp == 0) g_max_num_threads = l;
  }
#else
#error "unsupported platform"
#endif
}
