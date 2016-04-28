#ifndef COMMON_HEADERS_H
#define COMMON_HEADERS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <ctype.h>

#include <sys/time.h>

#if defined(__x86_64__)
#include <immintrin.h>
#include <emmintrin.h>
#elif defined(__POWER8_VECTOR__)
#include <omp.h>
#endif

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fenv.h>

#endif
