#!/usr/bin/env bash
# Rebuild Intel GKL's PairHMM (and gkl_utils) for arm64 (Apple Silicon) by mapping the AVX
# intrinsics to NEON with SIMDe, and install it to mavenLocal as gkl:0.9.1-arm64.
#
# IMPORTANT — correctness vs speed:
#   The native PairHMM (AVX on x86, NEON here) is a *vectorized float* kernel. It is NOT
#   bit-identical to the pure-Java PairHMM: GATK's own VectorPairHMMUnitTest compares it to the
#   reference only within tolerance 1e-5. On Apple Silicon GATK uses GKL's AVX (256-bit) path
#   (GKL excludes AVX-512 on __APPLE__), which we lower to NEON via SIMDe.
#   => Use this jar ONLY if you want native PairHMM SPEED and can accept results that match the
#      x86 result within ~1e-5 (i.e. NOT bit-identical). For bit-identical scores, do NOT use it;
#      let GATK fall back to the pure-Java PairHMM (the default when no arm64 GKL is present).
#   Enable with: -DgklVersion=0.9.1-arm64  (and -pairHMM AVX_LOGLESS_CACHING to force native).
#
# Requirements: JAVA_HOME=arm64 JDK 17, Xcode CLT (clang), git, curl, the upstream x86 gkl-0.9.1.jar
# available in the Gradle cache (GATK must have been built at least once).
set -euo pipefail

VER=0.9.1; ARM_VER="${VER}-arm64"
SIMDE_TAG="${SIMDE_TAG:-v0.8.2}"
WORK="${WORK:-$HOME/gatk-arm64-natives}"
: "${JAVA_HOME:?set JAVA_HOME to an arm64 JDK 17}"
export PATH="$JAVA_HOME/bin:$PATH"

mkdir -p "$WORK"; cd "$WORK"
[ -d gkl ]   || git clone --depth 1 --branch "$VER" https://github.com/Intel-HLS/GKL.git gkl
[ -d simde ] || git clone --depth 1 --branch "$SIMDE_TAG" https://github.com/simd-everywhere/simde.git simde
SIMDE="$WORK/simde"
cd gkl/src/main/native

# SIMDe shim: x86 intrinsic headers -> NEON
mkdir -p shim
printf '#include <simde/x86/avx512.h>\n#include <simde/x86/svml.h>\n' > shim/x86intrin.h
printf '#include <simde/x86/avx512.h>\n#include <simde/x86/svml.h>\n' > shim/immintrin.h

# Patch common/avx.h: on non-x86 skip cpuid and report the AVX (not AVX512) path (idempotent)
grep -q 'SIMDe/NEON' common/avx.h || perl -0pi -e '
  s/(#define COMMON_AVX_H\n)/$1\n#if !defined(__x86_64__) && !defined(__i386__)\n#include <x86intrin.h>\n#include <stdint.h>\ninline int  check_xcr0_ymm()      { return 1; }\ninline int  check_xcr0_zmm()      { return 0; }\ninline bool is_avx_supported()    { return true; }\ninline bool is_sse_supported()    { return true; }\ninline bool is_avx2_supported()   { return true; }\ninline bool is_avx512_supported() { return false; }\n#else\n/;
  s/(\n#endif\s*)$/\n#endif \/\/ !x86 (SIMDe\/NEON) vs x86\n#endif\n/;
' common/avx.h

CXXFLAGS="-std=c++14 -O3 -fPIC -shared -Ishim -I$SIMDE -Icommon -Ipairhmm -Iutils \
  -I$JAVA_HOME/include -I$JAVA_HOME/include/darwin -DSIMDE_ENABLE_NATIVE_ALIASES -march=armv8-a+simd"

# Apple path = AVX (256) only; no avx512_impl.cc
clang++ $CXXFLAGS pairhmm/IntelPairHmm.cc pairhmm/avx_impl.cc pairhmm/pairhmm_common.cc -o /tmp/libgkl_pairhmm.dylib
clang++ $CXXFLAGS utils/utils.cc -o /tmp/libgkl_utils.dylib
file /tmp/libgkl_pairhmm.dylib | grep -q arm64
file /tmp/libgkl_utils.dylib   | grep -q arm64

# Build the arm64 jar = upstream x86 jar with the two mac dylibs swapped to arm64
GKLJAR=$(find "$HOME/.gradle/caches" -name "gkl-${VER}.jar" | grep -vE "sources|javadoc" | head -1)
[ -n "$GKLJAR" ] || { echo "upstream gkl-${VER}.jar not in Gradle cache; build GATK once first"; exit 1; }
rm -rf /tmp/gklpatch && mkdir -p /tmp/gklpatch/com/intel/gkl/native
cp /tmp/libgkl_pairhmm.dylib /tmp/libgkl_utils.dylib /tmp/gklpatch/com/intel/gkl/native/
cp "$GKLJAR" "/tmp/gkl-${ARM_VER}.jar"
( cd /tmp/gklpatch && jar uf "/tmp/gkl-${ARM_VER}.jar" \
    com/intel/gkl/native/libgkl_pairhmm.dylib com/intel/gkl/native/libgkl_utils.dylib )

# Install to mavenLocal (reuse the upstream POM, version bumped)
DEST="$HOME/.m2/repository/com/intel/gkl/gkl/${ARM_VER}"
mkdir -p "$DEST"
cp "/tmp/gkl-${ARM_VER}.jar" "$DEST/gkl-${ARM_VER}.jar"
curl -fsSL "https://repo1.maven.org/maven2/com/intel/gkl/gkl/${VER}/gkl-${VER}.pom" \
  | sed "s#<version>${VER}</version>#<version>${ARM_VER}</version>#" > "$DEST/gkl-${ARM_VER}.pom"
echo "Installed com.intel.gkl:gkl:${ARM_VER} (NEON PairHMM) to mavenLocal"
