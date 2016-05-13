#include "headers.h"
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;

//static members from ConvertChar
uint8_t ConvertChar::conversionTable[255];
//Global function pointers in utils.h
float (*g_compute_full_prob_float)(testcase *tc, float* before_last_log) = 0;
double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log) = 0;


void initialize_function_pointers()
{
  g_compute_full_prob_float = compute_full_prob_avxs<float>;
  g_compute_full_prob_double = compute_full_prob_avxd<double>;
}
