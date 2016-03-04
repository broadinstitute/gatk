#include <iostream>

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

#define PRECISION s

#define MAIN_TYPE float
#define MAIN_TYPE_SIZE 32
#define UNION_TYPE mix_F
#define IF_128 IF_128f
#define IF_MAIN_TYPE IF_32
#define SHIFT_CONST1 12
#define SHIFT_CONST2 3
#define SHIFT_CONST3 4
#define _128_TYPE __m128
#define SIMD_TYPE __m256
#define _256_INT_TYPE __m256i
#define AVX_LENGTH 8
#define HAP_TYPE UNION_TYPE
#define MASK_TYPE uint32_t
#define MASK_ALL_ONES 0xFFFFFFFF
#define MASK_VEC MaskVec_F

#define SET_VEC_ZERO(__vec)                     \
    __vec= _mm256_setzero_ps()

#define VEC_OR(__v1, __v2)                      \
    _mm256_or_ps(__v1, __v2)

#define VEC_ADD(__v1, __v2)                     \
    _mm256_add_ps(__v1, __v2)

#define VEC_SUB(__v1, __v2)                     \
    _mm256_sub_ps(__v1, __v2)

#define VEC_MUL(__v1, __v2)                     \
    _mm256_mul_ps(__v1, __v2)

#define VEC_DIV(__v1, __v2)                     \
    _mm256_div_ps(__v1, __v2)

#define VEC_BLEND(__v1, __v2, __mask)           \
    _mm256_blend_ps(__v1, __v2, __mask)

#define VEC_BLENDV(__v1, __v2, __maskV)         \
    _mm256_blendv_ps(__v1, __v2, __maskV)

#define VEC_CAST_256_128(__v1)                  \
    _mm256_castps256_ps128 (__v1)

#define VEC_EXTRACT_128(__v1, __im)             \
    _mm256_extractf128_ps (__v1, __im)

#define VEC_EXTRACT_UNIT(__v1, __im)            \
    _mm_extract_epi32(__v1, __im)

#define VEC_SET1_VAL128(__val)                  \
    _mm_set1_ps(__val)

#define VEC_MOVE(__v1, __val)                   \
    _mm_move_ss(__v1, __val)

#define VEC_CAST_128_256(__v1)                  \
    _mm256_castps128_ps256(__v1)

#define VEC_INSERT_VAL(__v1, __val, __pos)      \
    _mm256_insertf128_ps(__v1, __val, __pos)

#define VEC_CVT_128_256(__v1)                   \
    _mm256_cvtepi32_ps(__v1.i)

#define VEC_SET1_VAL(__val)                     \
    _mm256_set1_ps(__val)

#define VEC_POPCVT_CHAR(__ch)                   \
    _mm256_cvtepi32_ps(_mm256_set1_epi32(__ch))

#define VEC_LDPOPCVT_CHAR(__addr)               \
    _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i const *)__addr))

#define VEC_CMP_EQ(__v1, __v2)                  \
    _mm256_cmp_ps(__v1, __v2, _CMP_EQ_OQ)

#define VEC_SET_LSE(__val)                      \
    _mm256_set_ps(zero, zero, zero, zero, zero, zero, zero, __val);

#define SHIFT_HAP(__v1, __val)                  \
    _vector_shift_lastavxs(__v1, __val.f);

#define VEC_SSE_TO_AVX(__vsLow, __vsHigh, __vdst)       \
    __vdst = _mm256_castps128_ps256(__vsLow) ;            \
__vdst = _mm256_insertf128_ps(__vdst, __vsHigh, 1) ;

#define VEC_SHIFT_LEFT_1BIT(__vs)               \
    __vs = _mm_slli_epi32(__vs, 1)

#define COMPARE_VECS(__v1, __v2, __first, __last) {                     \
    float* ptr1 = (float*) (&__v1) ;                                    \
    float* ptr2 = (float*) (&__v2) ;                                    \
    for (int ei=__first; ei <= __last; ++ei) {                          \
        if (ptr1[ei] != ptr2[ei]) {                                       \
            std::cout << "Float Mismatch at " << ei << ": "                 \
            << ptr1[ei] << " vs. " << ptr2[ei] << std::endl ;     \
            exit(0) ;                                                       \
        }                                                                 \
    }                                                                   \
}

class BitMaskVec_float {

    MASK_VEC low_, high_ ;
    SIMD_TYPE combined_ ;

    public:

    inline MASK_TYPE& getLowEntry(int index) {
        return low_.masks[index] ;
    }
    inline MASK_TYPE& getHighEntry(int index) {
        return high_.masks[index] ;
    }

    inline const SIMD_TYPE& getCombinedMask() {
        VEC_SSE_TO_AVX(low_.vecf, high_.vecf, combined_) ;
        return combined_ ;
    }

    inline void shift_left_1bit() {
        VEC_SHIFT_LEFT_1BIT(low_.vec) ;
        VEC_SHIFT_LEFT_1BIT(high_.vec) ;
    }

} ;

#define BITMASK_VEC BitMaskVec_float
