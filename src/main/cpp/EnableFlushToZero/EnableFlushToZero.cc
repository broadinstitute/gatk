#include <stdio.h>
#include <xmmintrin.h>

class EnableFlushToZero {
 public:
  EnableFlushToZero();
};

EnableFlushToZero::EnableFlushToZero() {
  printf("INFO JNI - Enable Flush To Zero mode.\n");
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
}

EnableFlushToZero ftz;
