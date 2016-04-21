#include <altivec.h>

typedef vector unsigned int __m128i;
typedef vector float __m128;
typedef vector double __m128d;
#define _mm_add_ps(a,b) vec_add(a,b)
#define _mm_add_pd(a,b) vec_add(a,b)
#define _mm_sub_ps(a,b) vec_sub(a,b)
#define _mm_sub_pd(a,b) vec_sub(a,b)
#define _mm_div_ps(a,b) vec_div(a,b)
#define _mm_div_pd(a,b) vec_div(a,b)
#define _mm_mul_ps(a,b) vec_mul(a,b)
#define _mm_mul_pd(a,b) vec_mul(a,b)
// r := (b==0)? a32_0 : ((b==1) ? a32_1 : ((b==2) ? a32_2 : a32_3))
#define _mm_extract_epi32(a,b) vec_extract(a,b<=3?b:3)
// r := (b==0)? a64_0 : a64_1
#define _mm_extract_epi64(a,b) vec_extract((vector unsigned long long)a,b<=1?b:1)
// r0 := (c==0) ? b : a32_0
// r1 := (c==1) ? b : a32_1
// r2 := (c==2) ? b : a32_2
// r3 := (c==3) ? b : a32_3
#define _mm_insert_epi32(a,b,c) (c<=3 ? vec_insert(b,a,c) : a)
// r0 := (c==0) ? b : a64_0
// r1 := (c==1) ? b : a64_1
#define _mm_insert_epi64(a,b,c) (c<=1 ? (__m128i)vec_insert(b,(vector unsigned long long)a,c) : a)
// r0 := a32_0 << b
// r1 := a32_1 << b
// r2 := a32_2 << b
// r3 := a32_3 << b
#define _mm_slli_epi32(a,b) vec_sl(a,((__m128i){b,b,b,b}))
// r0 := a64_0 << b
// r1 := a64_1 << b
#define _mm_slli_epi64(a,b) (__m128i)vec_sl((vector unsigned long long)a,((vector unsigned long long){b,b}))
// r0 := r1 := r2 := r3 := a32
#define _mm_set1_ps(a) ((__m128){a,a,a,a})// vec_splat((vector float){a}, 0) // http://lists.freebsd.org/pipermail/freebsd-ppc/2014-July/007110.html 
// r0 := r1 := a64
#define _mm_set1_pd(a) ((__m128d){a,a})// vec_splat((vector double){a}, 0)
// r0 := a32
// r1 := a32
// r2 := a32
// r3 := a32
#define _mm_set1_epi32(a) ((__m128i){a,a,a,a})//vec_splat((__m128i){a}, 0)
// r0 := d32
// r1 := c32
// r2 := b32
// r3 := a32
#define _mm_set_ps(a,b,c,d) ((__m128){d,c,b,a})
// r0 := b64
// r1 := a64
#define _mm_set_pd(a,b) ((__m128d){b,a})
// r0 := (float)a32_0
// r1 := (float)a32_1
// r2 := (float)a32_2
// r3 := (float)a32_3
#define _mm_cvtepi32_ps(a) vec_ctf(a,0) // http://www.filibeto.org/unix/macos/lib/dev/documentation/Performance/Conceptual/Accelerate_sse_migration/Accelerate_sse_migration.pdf
// r0 := (double)a64_0
// r1 := (double)a64_1
#if defined(vec_ctd)
#define _mm_cvtepi32_pd(a) vec_ctd(a,0) /* will be supported by a future version of gcc */
#else
#if __GNUC__ > 4
#define _mm_cvtepi32_pd(a) ((__m128d){(double)vec_extract(a,0),(double)vec_extract(a,1)})
#else
// Some g++ cannot compile this because of "sorry, unimplemented: unexpected AST of kind compound_literal_expr"
#if !defined(_MM_CVTEPI32_PD)
#define _MM_CVTEPI32_PD
extern "C" __inline __m128d _mm_cvtepi32_pd(vector unsigned int a) {
  double const a0 = (double)vec_extract((vector unsigned long long)a,0);
  double const a1 = (double)vec_extract((vector unsigned long long)a,1);
  return (__m128d){a0,a1};
}
#endif
#endif
#endif
// __m128i r := __m128i a << (int b * 8)
#define _mm_slli_si128(a,b) \
  (b==4?vec_perm(a,(__m128i){0,0,0,0},(vector unsigned char){ 16,16,16,16, 0,1,2,3, 4,5,6,7, 8,9,10,11 }): \
   (b==8?vec_perm(a,(__m128i){0,0,0,0},(vector unsigned char){ 16,16,16,16, 16,16,16,16, 0,1,2,3, 4,5,6,7}): \
	(__m128i){(exit(-1),0U)}))
// r0 := (c32_0 & 0x80000000) ? b32_0 : a32_0
// r1 := (c32_1 & 0x80000000) ? b32_1 : a32_1
// r2 := (c32_2 & 0x80000000) ? b32_2 : a32_2
// r3 := (c32_3 & 0x80000000) ? b32_3 : a32_3
//#define _mm_blendv_ps(a,b,c) vec_sel(a,b,(vector unsigned int){((unsigned int)c[0]&0x80000000)?0xFFFFFFFF:0,((unsigned int)c[1]&0x80000000)?0xFFFFFFFF:0,((unsigned int)c[2]&0x80000000)?0xFFFFFFFF:0,((unsigned int)c[3]&0x80000000)?0xFFFFFFFF:0})
#define _mm_blendv_ps(a,b,c) vec_sel(a,b,vec_sra((vector signed int)c, (vector unsigned int){31,31,31,31}))
// r0 := (c64_0 & 0x80000000) ? b64_0 : a64_0
// r1 := (c64_1 & 0x80000000) ? b64_1 : a64_1
//#define _mm_blendv_pd(a,b,c) vec_sel(a,b,(vector unsigned long long){(c[0]&0x8000000000000000UL)?0xFFFFFFFFFFFFFFFFUL:0,(c[1]&0x8000000000000000UL)?0xFFFFFFFFFFFFFFFFUL:0})
#define _mm_blendv_pd(a,b,c) vec_sel(a,b,vec_sra((vector signed long long)c, (vector unsigned long long){63,63}))
//#if !defined(_MM_BLENDV_PD)
//#define _MM_BLENDV_PD
//extern "C" __inline vector double _mm_blendv_pd(vector double a, vector double b, vector double c) {
//  vector signed long long c1 = (vector signed long long)c;
//  vector unsigned long long c2 = (vector unsigned long long){63,63};
//  vector signed long long d = vec_sra(c1, c2);
//  return vec_sel(a, b, d);
//}
//#endif
