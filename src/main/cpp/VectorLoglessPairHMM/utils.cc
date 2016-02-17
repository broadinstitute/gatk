#include "headers.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;

//static members from ConvertChar
uint8_t ConvertChar::conversionTable[255];
//Global function pointers in utils.h
float (*g_compute_full_prob_float)(testcase *tc, float* before_last_log) = 0;
double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log) = 0;
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


void initialize_function_pointers()
{
  g_compute_full_prob_float = compute_full_prob_avxs<float>;
  g_compute_full_prob_double = compute_full_prob_avxd<double>;
}
