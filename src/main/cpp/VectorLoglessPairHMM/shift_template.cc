#ifdef PRECISION

inline void CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn, MAIN_TYPE &shiftOut)
{
    IF_128 xlow , xhigh;
    /* cast x to xlow */
    xlow.f = VEC_CAST_256_128(x.d);
    /* extract x,1 to xhigh */
    xhigh.f = VEC_EXTRACT_128(x.d, 1);
    /* extract xlow[3] */
    IF_128 shiftOutL128;
    shiftOutL128.i = _mm_srli_si128(xlow.i, SHIFT_CONST1);
    /* extract xhigh[3] */
    IF_MAIN_TYPE shiftOutH;
    shiftOutH.i = VEC_EXTRACT_UNIT(xhigh.i, SHIFT_CONST2);
    shiftOut = shiftOutH.f;
    /* shift xlow */
    xlow.i = _mm_slli_si128 (xlow.i, SHIFT_CONST3);
    /* shift xhigh */
    xhigh.i = _mm_slli_si128 (xhigh.i, SHIFT_CONST3);
    /*movss shiftIn to xlow[0] */
    _128_TYPE shiftIn128 = VEC_SET1_VAL128(shiftIn);
    xlow.f = VEC_MOVE(xlow.f , shiftIn128);
    /*movss xlow[3] to xhigh[0] */
    xhigh.f = VEC_MOVE(xhigh.f, shiftOutL128.f);
    /* cast xlow to x */
    x.d = VEC_CAST_128_256(xlow.f);
    /* insert xhigh to x,1 */
    x.d = VEC_INSERT_VAL(x.d, xhigh.f, 1);
}


inline void CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION) (UNION_TYPE &x, MAIN_TYPE shiftIn)
{
    IF_128 xlow , xhigh;
    /* cast x to xlow */
    xlow.f = VEC_CAST_256_128(x.d);
    /* extract x,1 to xhigh */
    xhigh.f = VEC_EXTRACT_128(x.d, 1);
    /* extract xlow[3] */
    IF_128 shiftOutL128;
    shiftOutL128.i = _mm_srli_si128(xlow.i, SHIFT_CONST1);
    /* shift xlow */
    xlow.i = _mm_slli_si128 (xlow.i, SHIFT_CONST3);
    /* shift xhigh */
    xhigh.i = _mm_slli_si128 (xhigh.i, SHIFT_CONST3);
    /*movss shiftIn to xlow[0] */
    _128_TYPE shiftIn128 = VEC_SET1_VAL128(shiftIn);
    xlow.f = VEC_MOVE(xlow.f , shiftIn128);
    /*movss xlow[3] to xhigh[0] */
    xhigh.f = VEC_MOVE(xhigh.f, shiftOutL128.f);
    /* cast xlow to x */
    x.d = VEC_CAST_128_256(xlow.f);
    /* insert xhigh to x,1 */
    x.d = VEC_INSERT_VAL(x.d, xhigh.f, 1);
}

#endif
