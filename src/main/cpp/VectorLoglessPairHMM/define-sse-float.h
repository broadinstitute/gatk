#ifdef PRECISION
#undef PRECISION
#undef MAIN_TYPE
#undef MAIN_TYPE_SIZE
#undef UNION_TYPE
#undef IF_128
#undef IF_MAIN_TYPE
#undef SHIFT_CONST1
#undef SHIFT_CONST2
#undef SHIFT_CONST3
#undef _128_TYPE
#undef SIMD_TYPE
#undef AVX_LENGTH
#undef HAP_TYPE
#undef MASK_TYPE
#undef MASK_ALL_ONES

#undef VEC_EXTRACT_UNIT
#undef VEC_INSERT_UNIT
#undef SET_VEC_ZERO
#undef VEC_OR
#undef VEC_ADD
#undef VEC_SUB
#undef VEC_MUL
#undef VEC_DIV
#undef VEC_BLEND
#undef VEC_BLENDV
#undef VEC_CAST_256_128
#undef VEC_EXTRACT_128
#undef VEC_EXTRACT_UNIT
#undef VEC_SET1_VAL128
#undef VEC_MOVE
#undef VEC_CAST_128_256
#undef VEC_INSERT_VAL
#undef VEC_CVT_128_256
#undef VEC_SET1_VAL
#undef VEC_POPCVT_CHAR
#undef VEC_LDPOPCVT_CHAR
#undef VEC_CMP_EQ
#undef VEC_SET_LSE
#undef SHIFT_HAP
#undef MASK_VEC
#undef VEC_SSE_TO_AVX
#undef VEC_SHIFT_LEFT_1BIT
#undef MASK_ALL_ONES
#undef COMPARE_VECS
#undef _256_INT_TYPE
#undef BITMASK_VEC
#endif

#define SSE
#define PRECISION s

#define MAIN_TYPE float
#define MAIN_TYPE_SIZE 32
#define UNION_TYPE mix_F128
#define IF_128 IF_128f
#define IF_MAIN_TYPE IF_32
#define SHIFT_CONST1 3
#define SHIFT_CONST2 4
#define SHIFT_CONST3 0
#define _128_TYPE __m128
#define SIMD_TYPE __m128
#define _256_INT_TYPE __m128i
#define AVX_LENGTH 4
//#define MAVX_COUNT  (MROWS+3)/AVX_LENGTH
#define HAP_TYPE UNION_TYPE
#define MASK_TYPE uint32_t
#define MASK_ALL_ONES 0xFFFFFFFF
#define MASK_VEC MaskVec_F

#define VEC_EXTRACT_UNIT(__v1, __im)            \
    _mm_extract_epi32(__v1, __im)

#define VEC_INSERT_UNIT(__v1,__ins,__im)        \
    _mm_insert_epi32(__v1,__ins,__im)

#define VEC_OR(__v1, __v2)                      \
    _mm_or_ps(__v1, __v2)

#define VEC_ADD(__v1, __v2)                     \
    _mm_add_ps(__v1, __v2)

#define VEC_SUB(__v1, __v2)                     \
    _mm_sub_ps(__v1, __v2)

#define VEC_MUL(__v1, __v2)                     \
    _mm_mul_ps(__v1, __v2)

#define VEC_DIV(__v1, __v2)                     \
    _mm_div_ps(__v1, __v2)

#define VEC_CMP_EQ(__v1, __v2)                  \
    _mm_cmpeq_ps(__v1, __v2)

#define VEC_BLEND(__v1, __v2, __mask)           \
    _mm_blend_ps(__v1, __v2, __mask)

#define VEC_BLENDV(__v1, __v2, __maskV)         \
    _mm_blendv_ps(__v1, __v2, __maskV)

#define SHIFT_HAP(__v1, __val)                  \
    _vector_shift_lastsses(__v1, __val.f)

#define VEC_CVT_128_256(__v1)                   \
    _mm_cvtepi32_ps(__v1.i)

#define VEC_SET1_VAL(__val)                     \
    _mm_set1_ps(__val)

#define VEC_POPCVT_CHAR(__ch)                   \
    _mm_cvtepi32_ps(_mm_set1_epi32(__ch))

#define VEC_SET_LSE(__val)                      \
    _mm_set_ps(zero, zero, zero, __val);

#define VEC_LDPOPCVT_CHAR(__addr)               \
    _mm_cvtepi32_ps(_mm_loadu_si128((__m128i const *)__addr))

#define VEC_SSE_TO_AVX(__vsLow, __vsHigh, __vdst)       \
    __vdst = _mm_cvtpi32x2_ps(__vsLow, __vsHigh)

#define VEC_SHIFT_LEFT_1BIT(__vs)               \
    __vs = _mm_slli_epi32(__vs, 1)

class BitMaskVec_sse_float {

    MASK_VEC combined_ ;

    public:
    inline MASK_TYPE& getLowEntry(int index) {
        return combined_.masks[index] ;
    }
    inline MASK_TYPE& getHighEntry(int index) {
        return combined_.masks[AVX_LENGTH/2+index] ;
    }

    inline const SIMD_TYPE& getCombinedMask() {
        return combined_.vecf ;
    }

    inline void shift_left_1bit() {
        VEC_SHIFT_LEFT_1BIT(combined_.vec) ;
    }

} ;

#define BITMASK_VEC BitMaskVec_sse_float
