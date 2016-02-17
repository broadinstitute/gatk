#ifdef PRECISION

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>


void CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(const testcase& tc, int COLS, int numMaskVecs, MASK_TYPE (*maskArr)[NUM_DISTINCT_CHARS]) {

    const int maskBitCnt = MAIN_TYPE_SIZE ;

    for (int vi=0; vi < numMaskVecs; ++vi) {
        for (int rs=0; rs < NUM_DISTINCT_CHARS; ++rs) {
            maskArr[vi][rs] = 0 ;
        }
        maskArr[vi][AMBIG_CHAR] = MASK_ALL_ONES ;
    }

    for (int col=1; col < COLS; ++col) {
        int mIndex = (col-1) / maskBitCnt ;
        int mOffset = (col-1) % maskBitCnt ;
        MASK_TYPE bitMask = ((MASK_TYPE)0x1) << (maskBitCnt-1-mOffset) ;

        char hapChar = ConvertChar::get(tc.hap[col-1]);

        if (hapChar == AMBIG_CHAR) {
            for (int ci=0; ci < NUM_DISTINCT_CHARS; ++ci)
                maskArr[mIndex][ci] |= bitMask ;
        }

        maskArr[mIndex][hapChar] |= bitMask ;
        // bit corresponding to col 1 will be the MSB of the mask 0
        // bit corresponding to col 2 will be the MSB-1 of the mask 0
        // ...
        // bit corresponding to col 32 will be the LSB of the mask 0
        // bit corresponding to col 33 will be the MSB of the mask 1
        // ...
    }

}

void CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(const testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess) {

    for (int ri=0; ri < numRowsToProcess; ++ri) {
        rsArr[ri] = ConvertChar::get(tc.rs[ri+beginRowIndex-1]) ;
    }

    for (int ei=0; ei < AVX_LENGTH; ++ei) {
        lastMaskShiftOut[ei] = 0 ;
    }
}

#define SET_MASK_WORD(__dstMask, __srcMask, __lastShiftOut, __shiftBy, __maskBitCnt){ \
    MASK_TYPE __bitMask = (((MASK_TYPE)0x1) << __shiftBy) - 1 ;            \
    MASK_TYPE __nextShiftOut = (__srcMask & __bitMask) << (__maskBitCnt - __shiftBy) ; \
    __dstMask = (__srcMask >> __shiftBy) | __lastShiftOut ;        \
    __lastShiftOut = __nextShiftOut ;                    \
}


void CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)(int maskIndex, BITMASK_VEC& bitMaskVec, MASK_TYPE (*maskArr) [NUM_DISTINCT_CHARS], char* rsArr, MASK_TYPE* lastMaskShiftOut, int maskBitCnt) {

    for (int ei=0; ei < AVX_LENGTH/2; ++ei) {
        SET_MASK_WORD(bitMaskVec.getLowEntry(ei), maskArr[maskIndex][rsArr[ei]],
                lastMaskShiftOut[ei], ei, maskBitCnt) ;

        int ei2 = ei + AVX_LENGTH/2 ; // the second entry index
        SET_MASK_WORD(bitMaskVec.getHighEntry(ei), maskArr[maskIndex][rsArr[ei2]],
                lastMaskShiftOut[ei2], ei2, maskBitCnt) ;
    }

}


inline void CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (BITMASK_VEC& bitMaskVec, SIMD_TYPE& distm, SIMD_TYPE& _1_distm, SIMD_TYPE& distmChosen) {

    distmChosen = VEC_BLENDV(distm, _1_distm, bitMaskVec.getCombinedMask()) ;

    bitMaskVec.shift_left_1bit() ;
}

/*
 * This function:
 * 1- Intializes probability values p_MM, p_XX, P_YY, p_MX, p_GAPM and pack them into vectors
 * 2- Precompute parts of "distm" which only depeneds on a row number and pack it into vector
 */
 
template<class NUMBER> void CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)(int ROWS, int COLS, NUMBER* shiftOutM, NUMBER *shiftOutX, NUMBER *shiftOutY, Context<NUMBER> ctx, testcase *tc,  SIMD_TYPE *p_MM, SIMD_TYPE *p_GAPM, SIMD_TYPE *p_MX, SIMD_TYPE *p_XX, SIMD_TYPE *p_MY, SIMD_TYPE *p_YY, SIMD_TYPE *distm1D)
{
    NUMBER zero = ctx._(0.0);
    NUMBER init_Y = ctx.INITIAL_CONSTANT / (tc->haplen);
    for (int s=0;s<ROWS+COLS+AVX_LENGTH;s++)
    {
        shiftOutM[s] = zero;
        shiftOutX[s] = zero;
        shiftOutY[s] = init_Y;
    }

    NUMBER *ptr_p_MM = (NUMBER *)p_MM;
    NUMBER *ptr_p_XX = (NUMBER *)p_XX;
    NUMBER *ptr_p_YY = (NUMBER *)p_YY;
    NUMBER *ptr_p_MX = (NUMBER *)p_MX;
    NUMBER *ptr_p_MY = (NUMBER *)p_MY;
    NUMBER *ptr_p_GAPM = (NUMBER *)p_GAPM;

    *ptr_p_MM = ctx._(0.0);
    *ptr_p_XX = ctx._(0.0);
    *ptr_p_YY = ctx._(0.0);
    *ptr_p_MX = ctx._(0.0);
    *ptr_p_MY = ctx._(0.0);
    *ptr_p_GAPM = ctx._(0.0);

    for (int r = 1; r < ROWS; r++)
    {
        int _i = tc->i[r-1] & 127;
        int _d = tc->d[r-1] & 127;
        int _c = tc->c[r-1] & 127;

        //*(ptr_p_MM+r-1) = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];
        SET_MATCH_TO_MATCH_PROB(*(ptr_p_MM+r-1), _i, _d);
        *(ptr_p_GAPM+r-1) = ctx._(1.0) - ctx.ph2pr[_c];
        *(ptr_p_MX+r-1) = ctx.ph2pr[_i];
        *(ptr_p_XX+r-1) = ctx.ph2pr[_c];
        *(ptr_p_MY+r-1) = ctx.ph2pr[_d];
        *(ptr_p_YY+r-1) = ctx.ph2pr[_c];
    }

    NUMBER *ptr_distm1D = (NUMBER *)distm1D;
    for (int r = 1; r < ROWS; r++)
    {
        int _q = tc->q[r-1] & 127;
        ptr_distm1D[r-1] = ctx.ph2pr[_q];
    }
}

/*
 * This function handles pre-stripe computation:
 * 1- Retrieve probaility vectors from memory 
 * 2- Initialize M, X, Y vectors with all 0's (for the first stripe) and shifting the last row from previous stripe for the rest 
 */

template<class NUMBER> inline void CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(
        int stripeIdx, Context<NUMBER> ctx, testcase *tc, SIMD_TYPE &pGAPM, SIMD_TYPE &pMM, SIMD_TYPE &pMX, SIMD_TYPE &pXX, SIMD_TYPE &pMY, SIMD_TYPE &pYY,
        SIMD_TYPE &rs, UNION_TYPE &rsN, SIMD_TYPE &distm, SIMD_TYPE &_1_distm,  SIMD_TYPE *distm1D, SIMD_TYPE N_packed256, SIMD_TYPE *p_MM , SIMD_TYPE *p_GAPM ,
        SIMD_TYPE *p_MX, SIMD_TYPE *p_XX , SIMD_TYPE *p_MY, SIMD_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1,
        UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER* shiftOutX, NUMBER* shiftOutM)
{
    int i = stripeIdx;
    pGAPM = p_GAPM[i];
    pMM   = p_MM[i];
    pMX   = p_MX[i];
    pXX   = p_XX[i];
    pMY   = p_MY[i];
    pYY   = p_YY[i];

    NUMBER zero = ctx._(0.0);
    NUMBER init_Y = ctx.INITIAL_CONSTANT / (tc->haplen);
    UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
    UNION_TYPE packed3;  packed3.d = VEC_SET1_VAL(3.0);

    distm = distm1D[i];
    _1_distm = VEC_SUB(packed1.d, distm);

    distm = VEC_DIV(distm, packed3.d);

    /* initialize M_t_2, M_t_1, X_t_2, X_t_1, Y_t_2, Y_t_1 */
    M_t_2.d = VEC_SET1_VAL(zero);
    X_t_2.d = VEC_SET1_VAL(zero);

    if (i==0) {
        M_t_1.d = VEC_SET1_VAL(zero);
        X_t_1.d = VEC_SET1_VAL(zero);
        Y_t_2.d = VEC_SET_LSE(init_Y);
        Y_t_1.d = VEC_SET1_VAL(zero);
    }
    else {
        X_t_1.d = VEC_SET_LSE(shiftOutX[AVX_LENGTH]);
        M_t_1.d = VEC_SET_LSE(shiftOutM[AVX_LENGTH]);
        Y_t_2.d = VEC_SET1_VAL(zero);
        Y_t_1.d = VEC_SET1_VAL(zero);
    }
    M_t_1_y = M_t_1;
}

/*
 *  This function is the main compute kernel to compute M, X and Y
 */

inline void CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(UNION_TYPE &M_t, UNION_TYPE &X_t, UNION_TYPE &Y_t, UNION_TYPE &M_t_y,
        UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1,
        SIMD_TYPE pMM, SIMD_TYPE pGAPM, SIMD_TYPE pMX, SIMD_TYPE pXX, SIMD_TYPE pMY, SIMD_TYPE pYY, SIMD_TYPE distmSel)
{
    /* Compute M_t <= distm * (p_MM*M_t_2 + p_GAPM*X_t_2 + p_GAPM*Y_t_2) */
    M_t.d = VEC_MUL(VEC_ADD(VEC_ADD(VEC_MUL(M_t_2.d, pMM), VEC_MUL(X_t_2.d, pGAPM)), VEC_MUL(Y_t_2.d, pGAPM)), distmSel);
    //M_t.d = VEC_MUL( VEC_ADD(VEC_MUL(M_t_2.d, pMM), VEC_MUL(VEC_ADD(X_t_2.d, Y_t_2.d), pGAPM)), distmSel);

    M_t_y = M_t;

    /* Compute X_t */
    X_t.d = VEC_ADD(VEC_MUL(M_t_1.d, pMX) , VEC_MUL(X_t_1.d, pXX));

    /* Compute Y_t */
    Y_t.d = VEC_ADD(VEC_MUL(M_t_1_y.d, pMY) , VEC_MUL(Y_t_1.d, pYY));
}

/*
 * This is the main compute function. It operates on the matrix in s stripe manner.
 * The stripe height is determined by the SIMD engine type. 
 * Stripe height: "AVX float": 8, "AVX double": 4
 * For each stripe the operations are anti-diagonal based. 
 * Each anti-diagonal (M_t, Y_t, X_t) depends on the two previous anti-diagonals (M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, Y_t_1).
 * Each stripe (except the fist one) depends on the last row of the previous stripe.
 * The last stripe computation handles the addition of the last row of M and X, that's the reason for loop spliting.
 */

template<class NUMBER> NUMBER CONCAT(CONCAT(compute_full_prob_,SIMD_ENGINE), PRECISION) (testcase *tc, NUMBER *before_last_log = NULL)
{
    int ROWS = tc->rslen + 1;
    int COLS = tc->haplen + 1;
    int MAVX_COUNT = (ROWS+AVX_LENGTH-1)/AVX_LENGTH;

    /* Probaility arrays */
    SIMD_TYPE p_MM   [MAVX_COUNT], p_GAPM [MAVX_COUNT], p_MX   [MAVX_COUNT];
    SIMD_TYPE p_XX   [MAVX_COUNT], p_MY   [MAVX_COUNT], p_YY   [MAVX_COUNT];

    /* For distm precomputation */
    SIMD_TYPE distm1D[MAVX_COUNT];

    /* Carries the values from each stripe to the next stripe */
    NUMBER shiftOutM[ROWS+COLS+AVX_LENGTH], shiftOutX[ROWS+COLS+AVX_LENGTH], shiftOutY[ROWS+COLS+AVX_LENGTH];

    /* The vector to keep the anti-diagonals of M, X, Y*/
    /* Current: M_t, X_t, Y_t */
    /* Previous: M_t_1, X_t_1, Y_t_1 */
    /* Previous to previous: M_t_2, X_t_2, Y_t_2 */ 
    UNION_TYPE  M_t, M_t_1, M_t_2, X_t, X_t_1, X_t_2, Y_t, Y_t_1, Y_t_2, M_t_y, M_t_1_y;

    /* Probality vectors */
    SIMD_TYPE pGAPM, pMM, pMX, pXX, pMY, pYY;

    struct timeval start, end;
    NUMBER result_avx2;
    Context<NUMBER> ctx;
    UNION_TYPE rs , rsN;
    HAP_TYPE hap;
    SIMD_TYPE distmSel, distmChosen ;
    SIMD_TYPE distm, _1_distm;

    int r, c;
    NUMBER zero = ctx._(0.0);
    UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
    SIMD_TYPE N_packed256 = VEC_POPCVT_CHAR('N');
    NUMBER init_Y = ctx.INITIAL_CONSTANT / (tc->haplen);
    int remainingRows = (ROWS-1) % AVX_LENGTH;
    int stripe_cnt = ((ROWS-1) / AVX_LENGTH) + (remainingRows!=0);

    const int maskBitCnt = MAIN_TYPE_SIZE ;
    const int numMaskVecs = (COLS+ROWS+maskBitCnt-1)/maskBitCnt ; // ceil function

    /* Mask precomputation for distm*/
    MASK_TYPE maskArr[numMaskVecs][NUM_DISTINCT_CHARS] ;
    CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(*tc, COLS, numMaskVecs, maskArr) ;

    char rsArr[AVX_LENGTH] ;
    MASK_TYPE lastMaskShiftOut[AVX_LENGTH] ;

    /* Precompute initialization for probabilities and shift vector*/
    CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)<NUMBER>(ROWS, COLS, shiftOutM, shiftOutX, shiftOutY,
            ctx, tc, p_MM, p_GAPM, p_MX, p_XX, p_MY, p_YY, distm1D);

    for (int i=0;i<stripe_cnt-1;i++)
    {
        //STRIPE_INITIALIZATION
        CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM);
        CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(*tc, rsArr, lastMaskShiftOut, i*AVX_LENGTH+1, AVX_LENGTH) ;
        // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors

        BITMASK_VEC bitMaskVec ;

        for (int begin_d=1;begin_d<COLS+AVX_LENGTH;begin_d+=MAIN_TYPE_SIZE)
        {
            int numMaskBitsToProcess = std::min(MAIN_TYPE_SIZE, COLS+AVX_LENGTH-begin_d) ;
            CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)((begin_d-1)/MAIN_TYPE_SIZE, bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt) ;

            for (int mbi=0; mbi < numMaskBitsToProcess; ++mbi) {
                CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (bitMaskVec, distm, _1_distm, distmChosen) ;
                int ShiftIdx = begin_d + mbi + AVX_LENGTH;

                CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                        pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(M_t, shiftOutM[ShiftIdx], shiftOutM[begin_d+mbi]);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(X_t, shiftOutX[ShiftIdx], shiftOutX[begin_d+mbi]);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(Y_t_1, shiftOutY[ShiftIdx], shiftOutY[begin_d+mbi]);

                M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;
            }
        }
    }

    int i = stripe_cnt-1;
    {
        //STRIPE_INITIALIZATION
        CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM);

        if (remainingRows==0) remainingRows=AVX_LENGTH;
        CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(*tc, rsArr, lastMaskShiftOut, i*AVX_LENGTH+1, remainingRows) ;

        SIMD_TYPE sumM, sumX;
        sumM = VEC_SET1_VAL(zero);
        sumX = VEC_SET1_VAL(zero);

        // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors
        BITMASK_VEC bitMaskVec ;

        for (int begin_d=1;begin_d<COLS+remainingRows-1;begin_d+=MAIN_TYPE_SIZE)
        {
            int numMaskBitsToProcess = std::min(MAIN_TYPE_SIZE, COLS+remainingRows-1-begin_d) ;
            CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE),PRECISION)((begin_d-1)/MAIN_TYPE_SIZE, bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt) ;

            for (int mbi=0; mbi < numMaskBitsToProcess; ++mbi) {

                CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (bitMaskVec, distm, _1_distm, distmChosen) ;
                int ShiftIdx = begin_d + mbi +AVX_LENGTH;

                CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                        pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                sumM  = VEC_ADD(sumM, M_t.d);
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(M_t, shiftOutM[ShiftIdx]);

                sumX  = VEC_ADD(sumX, X_t.d);
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(X_t, shiftOutX[ShiftIdx]);

                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(Y_t_1, shiftOutY[ShiftIdx]);

                M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;

            }
        }
        UNION_TYPE sumMX;
        sumMX.d = VEC_ADD(sumM, sumX);
        result_avx2 = sumMX.f[remainingRows-1];
    }
    return result_avx2;
}

#endif

