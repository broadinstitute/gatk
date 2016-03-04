#ifndef PAIRHMM_UTIL_H
#define PAIRHMM_UTIL_H

#include "common_data_structure.h"

template<class T>
std::string to_string(T obj)
{
  std::stringstream ss;
  std::string ret_string;
  ss.clear();
  ss << std::scientific << obj;
  ss >> ret_string;
  ss.clear();
  return ret_string;
}
void debug_dump(std::string filename, std::string s, bool to_append, bool add_newline=true);

int read_mod_testcase(std::ifstream& fptr, testcase* tc, bool reformat=false);

extern float (*g_compute_full_prob_float)(testcase *tc, float *before_last_log);
extern double (*g_compute_full_prob_double)(testcase *tc, double* before_last_log);
template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log=0);
template<class NUMBER>
NUMBER compute_full_prob_avxd(testcase *tc, NUMBER *before_last_log=0);
template<class NUMBER>
NUMBER compute_full_prob_avxs(testcase *tc, NUMBER *before_last_log=0);

void initialize_function_pointers();

#endif
