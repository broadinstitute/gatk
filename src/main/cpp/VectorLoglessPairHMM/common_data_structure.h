#ifndef COMMON_DATA_STRUCTURE_H
#define COMMON_DATA_STRUCTURE_H

#include "headers.h"

#define CAT(X,Y) X##Y
#define CONCAT(X,Y) CAT(X,Y)

#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5

//#define MROWS  500
//#define MCOLS  1000


#define MAX_QUAL 254
#define MAX_JACOBIAN_TOLERANCE 8.0
#define JACOBIAN_LOG_TABLE_STEP 0.0001
#define JACOBIAN_LOG_TABLE_INV_STEP (1.0 / JACOBIAN_LOG_TABLE_STEP)
#define MAXN 70000
#define LOG10_CACHE_SIZE  (4*MAXN)  // we need to be able to go up to 2*(2N) when calculating some of the coefficients
#define JACOBIAN_LOG_TABLE_SIZE ((int) (MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1)

template<class NUMBER>
struct ContextBase
{
  public:
    NUMBER ph2pr[128];
    NUMBER INITIAL_CONSTANT;
    NUMBER LOG10_INITIAL_CONSTANT;
    NUMBER RESULT_THRESHOLD; 

    static bool staticMembersInitializedFlag;
    static NUMBER jacobianLogTable[JACOBIAN_LOG_TABLE_SIZE];
    static NUMBER matchToMatchProb[((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1];

    static void initializeStaticMembers()
    {
      if(!staticMembersInitializedFlag)
      {
        //Order of calls important - Jacobian first, then MatchToMatch
        initializeJacobianLogTable();
        initializeMatchToMatchProb();
        staticMembersInitializedFlag = true;
      }
    }

    static void deleteStaticMembers()
    {
      if(staticMembersInitializedFlag)
      {
        staticMembersInitializedFlag = false;
      }
    }

    //Called only once during library load - don't bother to optimize with single precision fp
    static void initializeJacobianLogTable()
    {
      for (int k = 0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
        jacobianLogTable[k] = (NUMBER)(log10(1.0 + pow(10.0, -((double) k) * JACOBIAN_LOG_TABLE_STEP)));
      }
    }

    //Called only once per library load - don't bother optimizing with single fp
    static void initializeMatchToMatchProb()
    {
      double LN10 = log(10);
      double INV_LN10 = 1.0/LN10;
      for (int i = 0, offset = 0; i <= MAX_QUAL; offset += ++i)
        for (int j = 0; j <= i; j++) {
          double log10Sum = approximateLog10SumLog10(-0.1*i, -0.1*j);
          double matchToMatchLog10 =
            log1p(-std::min(1.0,pow(10,log10Sum))) * INV_LN10;
          matchToMatchProb[offset + j] = (NUMBER)(pow(10,matchToMatchLog10));
        }
    }
    //Called during computation - use single precision where possible
    static int fastRound(NUMBER d) {
      return (d > ((NUMBER)0.0)) ? (int) (d + ((NUMBER)0.5)) : (int) (d - ((NUMBER)0.5));
    }
    //Called during computation - use single precision where possible
    static NUMBER approximateLog10SumLog10(NUMBER small, NUMBER big) {
      // make sure small is really the smaller value
      if (small > big) {
        NUMBER t = big;
        big = small;
        small = t;
      }

      if (std::isinf(small) == -1 || std::isinf(big) == -1)
        return big;

      NUMBER diff = big - small;
      if (diff >= ((NUMBER)MAX_JACOBIAN_TOLERANCE))
        return big;

      // OK, so |y-x| < tol: we use the following identity then:
      // we need to compute log10(10^x + 10^y)
      // By Jacobian logarithm identity, this is equal to
      // max(x,y) + log10(1+10^-abs(x-y))
      // we compute the second term as a table lookup with integer quantization
      // we have pre-stored correction for 0,0.1,0.2,... 10.0
      int ind = fastRound((NUMBER)(diff * ((NUMBER)JACOBIAN_LOG_TABLE_INV_STEP))); // hard rounding
      return big + jacobianLogTable[ind];
    }
};

template<class NUMBER>
struct Context : public ContextBase<NUMBER>
{};

template<>
struct Context<double> : public ContextBase<double>
{
  Context():ContextBase<double>()
  {
    for (int x = 0; x < 128; x++)
      ph2pr[x] = pow(10.0, -((double)x) / 10.0);

    INITIAL_CONSTANT = ldexp(1.0, 1020.0);
    LOG10_INITIAL_CONSTANT = log10(INITIAL_CONSTANT);
    RESULT_THRESHOLD = 0.0;
  }

  double LOG10(double v){ return log10(v); }
  inline double POW(double b, double e) { return pow(b,e); }

  static double _(double n){ return n; }
  static double _(float n){ return ((double) n); }
};

template<>
struct Context<float> : public ContextBase<float>
{
  Context() : ContextBase<float>()
  {
    for (int x = 0; x < 128; x++)
    {
      ph2pr[x] = powf(10.f, -((float)x) / 10.f);
    }

    INITIAL_CONSTANT = ldexpf(1.f, 120.f);
    LOG10_INITIAL_CONSTANT = log10f(INITIAL_CONSTANT);
    RESULT_THRESHOLD = ldexpf(1.f, -110.f);
  }

  float LOG10(float v){ return log10f(v); }
  inline float POW(float b, float e) { return powf(b,e); }

  static float _(double n){ return ((float) n); }
  static float _(float n){ return n; }
};

#define SET_MATCH_TO_MATCH_PROB(output, insQual, delQual)                       \
{                                                                               \
  int minQual = delQual;                                                        \
  int maxQual = insQual;                                                        \
  if (insQual <= delQual)                                                       \
  {                                                                             \
    minQual = insQual;                                                          \
    maxQual = delQual;                                                          \
  }                                                                             \
  (output) = (MAX_QUAL < maxQual) ?                                             \
  ((NUMBER)1.0) - ctx.POW(((NUMBER)10), ctx.approximateLog10SumLog10(((NUMBER)-0.1)*minQual, ((NUMBER)-0.1)*maxQual))       \
  : ctx.matchToMatchProb[((maxQual * (maxQual + 1)) >> 1) + minQual];           \
}



typedef struct
{
        int rslen, haplen;
	/*int *q, *i, *d, *c;*/
	/*int q[MROWS], i[MROWS], d[MROWS], c[MROWS];*/
	char *q, *i, *d, *c;
        char *hap, *rs;
	int *ihap;
	int *irs;
} testcase;

#define MIN_ACCEPTED 1e-28f
#define NUM_DISTINCT_CHARS 5
#define AMBIG_CHAR 4

class ConvertChar {

  static uint8_t conversionTable[255] ;

public:

  static void init() {
    assert (NUM_DISTINCT_CHARS == 5) ;
    assert (AMBIG_CHAR == 4) ;

    conversionTable['A'] = 0 ;
    conversionTable['C'] = 1 ;
    conversionTable['T'] = 2 ;
    conversionTable['G'] = 3 ;
    conversionTable['N'] = 4 ;
  }

  static inline uint8_t get(uint8_t input) {
    return conversionTable[input] ;
  }

};

int normalize(char c);
int read_testcase(testcase *tc, FILE* ifp=0);

#endif
