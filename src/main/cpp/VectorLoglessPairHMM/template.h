#ifndef TEMPLATES_H_
#define TEMPLATES_H_

#include "headers.h"


#define ALIGNED __attribute__((aligned(32)))

typedef union __attribute__((aligned(32))) {
        ALIGNED __m256 ALIGNED d;
        ALIGNED __m128i ALIGNED s[2];
        ALIGNED float  ALIGNED f[8];
        ALIGNED __m256i ALIGNED i;
} ALIGNED mix_F ALIGNED;

typedef union __attribute__((aligned(32))) {
        ALIGNED __m128 ALIGNED d;
        ALIGNED __m64 ALIGNED s[2];
        ALIGNED float  ALIGNED f[4];
        ALIGNED __m128i ALIGNED i;
} ALIGNED mix_F128 ALIGNED;

typedef union ALIGNED {
  __m128i vec ;
  __m128 vecf ;
  uint32_t masks[4] ;
} MaskVec_F ;

typedef union ALIGNED {
  __m64 vec ;
  __m64 vecf ;
  uint32_t masks[2] ;
} MaskVec_F128 ;

typedef union ALIGNED
{
        ALIGNED __m128i ALIGNED i;
        ALIGNED __m128 ALIGNED f;
} ALIGNED IF_128f ALIGNED;

typedef union ALIGNED
{
        ALIGNED int    ALIGNED i;
        ALIGNED float  ALIGNED f;
} ALIGNED IF_32 ALIGNED;

typedef union __attribute__((aligned(32))) {
        ALIGNED __m256d ALIGNED d;
        ALIGNED __m128i ALIGNED s[2];
        ALIGNED double  ALIGNED f[4];
        ALIGNED __m256i ALIGNED i;
} ALIGNED mix_D ALIGNED;

typedef union __attribute__((aligned(32))) {
        ALIGNED __m128d ALIGNED d;
        ALIGNED __m64 ALIGNED s[2];
        ALIGNED double  ALIGNED f[2];
        ALIGNED __m128i ALIGNED i;
} ALIGNED mix_D128 ALIGNED;

typedef union ALIGNED {
  __m128i vec ;
  __m128d vecf ;
  uint64_t masks[2] ;
} MaskVec_D ;

typedef union ALIGNED {
  __m64 vec ;
  __m64 vecf ;
  uint64_t masks[1] ;
} MaskVec_D128 ;

typedef union ALIGNED
{
        ALIGNED __m128i ALIGNED i;
        ALIGNED __m128d ALIGNED f;
} ALIGNED IF_128d ALIGNED;

typedef union ALIGNED
{
        ALIGNED int64_t ALIGNED i;
        ALIGNED double  ALIGNED f;
} ALIGNED IF_64 ALIGNED;


#include "common_data_structure.h"

#endif


