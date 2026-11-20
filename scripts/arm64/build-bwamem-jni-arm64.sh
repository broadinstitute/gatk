#!/usr/bin/env bash
# Rebuild gatk-bwamem-jni for arm64 (Apple Silicon) and install it to the local Maven repo
# as version "<ver>-arm64", so GATK can consume it via -DbwaMemJni.version=1.0.4-arm64
# (or -DuseArm64Natives=true). See docs/arm64/native-libs.md.
#
# Strategy: bwa's only SIMD code is the integer Smith-Waterman in bwa/ksw.c (SSE2). We map the
# SSE2 intrinsics to NEON with sse2neon (integer ops -> bit-identical alignments). The JNI
# loader is patched to accept os.arch=aarch64/arm64. Validated by the repo's own BwaMemIndexTest
# (exact CIGAR/SNV assertions generated on x86) passing on arm64.
#
# Requirements: JAVA_HOME pointing at an arm64 JDK 17, Xcode CLT (clang), make, git, curl.
set -euo pipefail

VER="${VER:-1.0.4}"
ARM_VER="${VER}-arm64"
SSE2NEON_TAG="${SSE2NEON_TAG:-v1.8.0}"
WORK="${WORK:-$HOME/gatk-arm64-natives}"
M2="$HOME/.m2/repository/org/broadinstitute/gatk-bwamem-jni/${ARM_VER}"
: "${JAVA_HOME:?set JAVA_HOME to an arm64 JDK 17}"
export PATH="$JAVA_HOME/bin:$PATH"

mkdir -p "$WORK"; cd "$WORK"
[ -d bwamem ] || git clone --depth 1 --branch "$VER" https://github.com/broadinstitute/gatk-bwamem-jni.git bwamem
cd bwamem/src/main/c

# sse2neon shim: make <emmintrin.h>/<xmmintrin.h> resolve to NEON
mkdir -p sse2neon
[ -f sse2neon/sse2neon.h ] || curl -fsSL -o sse2neon/sse2neon.h \
  "https://raw.githubusercontent.com/DLTcollab/sse2neon/${SSE2NEON_TAG}/sse2neon.h"
printf '#include "sse2neon.h"\n' > sse2neon/emmintrin.h
printf '#include "sse2neon.h"\n' > sse2neon/xmmintrin.h

# Patch Makefile: clang + arm64 flags (idempotent)
grep -q 'sse2neon' Makefile || perl -0pi -e '
  s/^CC=gcc$/CC=cc/m;
  s{(UNAME := \$\(shell uname\))}{$1\nARCH := \$(shell uname -m)\nifneq (,\$(filter \$(ARCH),arm64 aarch64))\nCFLAGS += -I\$(CURDIR)/sse2neon -march=armv8-a+simd -flax-vector-conversions -DSSE2NEON_SUPPRESS_WARNINGS\nendif}m;
' Makefile

# Patch JNI loader to accept aarch64/arm64 (idempotent)
LOADER="../java/org/broadinstitute/hellbender/utils/bwa/BwaMemIndex.java"
grep -q 'aarch64' "$LOADER" || perl -0pi -e '
  s/(if \( !"x86_64"\.equals\(osArch\) && !"amd64"\.equals\(osArch\))( \))/$1\n                                && !"aarch64".equals(osArch) && !"arm64".equals(osArch)$2/;
' "$LOADER"

make clean >/dev/null 2>&1 || true
make all
file libbwa.Darwin.dylib | grep -q arm64 || { echo "ERROR: not arm64"; exit 1; }

# Assemble jar: compiled classes + dylib at archive root
cd "$WORK/bwamem"
rm -rf out jarstage && mkdir out jarstage
javac -d out $(find src/main/java -name '*.java')
cp -r out/* jarstage/
cp src/main/c/libbwa.Darwin.dylib jarstage/
( cd jarstage && jar cf "/tmp/gatk-bwamem-jni-${ARM_VER}.jar" . )

# Install to mavenLocal
mkdir -p "$M2"
cp "/tmp/gatk-bwamem-jni-${ARM_VER}.jar" "$M2/gatk-bwamem-jni-${ARM_VER}.jar"
cat > "$M2/gatk-bwamem-jni-${ARM_VER}.pom" <<POM
<?xml version="1.0" encoding="UTF-8"?>
<project xmlns="http://maven.apache.org/POM/4.0.0">
  <modelVersion>4.0.0</modelVersion>
  <groupId>org.broadinstitute</groupId>
  <artifactId>gatk-bwamem-jni</artifactId>
  <version>${ARM_VER}</version>
  <packaging>jar</packaging>
</project>
POM
echo "Installed org.broadinstitute:gatk-bwamem-jni:${ARM_VER} to mavenLocal"
