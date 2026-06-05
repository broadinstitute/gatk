#!/usr/bin/env bash
# Attempt to rebuild the jhdf5 JNI (libjhdf5.2.11.0.dylib) for arm64 by compiling the HDF5 1.10.11
# Java JNI sources (renamed to the ncsa.hdf.hdf5lib namespace) against an arm64 HDF5 1.10.
#
# STATUS: PARTIAL / NOT VALIDATED. The resulting dylib loads on Apple Silicon, resolves all 1249
# ncsa JNI symbols, and CREATES HDF5 files, but HDF5LibraryUnitTest still fails 10/12: group
# create + flush raise "Inappropriate type". Root cause: the bundled jarhdf5-2.11.0.jar (CISD
# jhdf5 2.11.0) Java classes were built against an OLDER HDF5 than 1.10.11, so some deprecated-API
# / enum behaviors differ at runtime. For a CORRECT arm64 native, build the CISD jhdf5 2.11.0
# source against ITS matching HDF5 version (it emits exactly these ncsa symbols), then repackage.
# Until then GATK holds hdf5-java-bindings at the upstream x86 jar (default); on aarch64 the CNV/
# HDF5 tools (CollectReadCounts, DenoiseReadCounts, ...) are not native yet. See docs/arm64/native-libs.md.
#
# This script documents the (reproducible) groundwork; it does NOT install a validated jar.
set -euo pipefail
: "${JAVA_HOME:?set JAVA_HOME to an arm64 JDK 17}"
WORK="${WORK:-$HOME/gatk-arm64-natives}"; mkdir -p "$WORK"; cd "$WORK"

brew list hdf5@1.10 >/dev/null 2>&1 || brew install hdf5@1.10
brew list libaec    >/dev/null 2>&1 || brew install libaec
H510=/opt/homebrew/opt/hdf5@1.10
SZ=$(ls /opt/homebrew/opt/libaec/lib/libsz*.dylib | head -1)

[ -d hdf5src ] || { mkdir hdf5src; curl -fsSL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.11/src/hdf5-1.10.11.tar.gz | tar xz -C hdf5src; }
SRC="$WORK/hdf5src/hdf5-1.10.11"
rm -rf hdf5jni && mkdir hdf5jni
cp "$SRC"/java/src/jni/*.c "$SRC"/java/src/jni/*.h hdf5jni/
cd hdf5jni

# Rename hdf.hdf5lib -> ncsa.hdf.hdf5lib in JNI function names and class-path string literals
perl -0pi -e 's/Java_hdf_hdf5lib/Java_ncsa_hdf_hdf5lib/g; s{"hdf/hdf5lib}{"ncsa/hdf/hdf5lib}g;' *.c *.h
# Add 3 legacy constants present in jhdf5 2.11.0 but dropped from HDF5 1.10.11
perl -0pi -e 's/(\nH5_GCC_DIAG_ON\("missing-prototypes"\))/\nJNIEXPORT jint JNICALL Java_ncsa_hdf_hdf5lib_HDF5Constants_H5D_1CHUNK_1BTREE(JNIEnv *e, jclass c){return 0;}\nJNIEXPORT jint JNICALL Java_ncsa_hdf_hdf5lib_HDF5Constants_H5F_1ACC_1DEBUG(JNIEnv *e, jclass c){return 0;}\nJNIEXPORT jint JNICALL Java_ncsa_hdf_hdf5lib_HDF5Constants_H5F_1SCOPE_1DOWN(JNIEnv *e, jclass c){return 2;}\n$1/' h5Constants.c

for c in *.c; do
  clang -c -O2 -fPIC -I. -I"$SRC/src" -I"$H510/include" \
    -I"$JAVA_HOME/include" -I"$JAVA_HOME/include/darwin" "$c" -o "${c%.c}.o"
done
clang -dynamiclib -install_name libjhdf5.2.11.0.dylib -o libjhdf5.2.11.0.dylib \
  *.o "$H510/lib/libhdf5.a" "$SZ" -lz -lm
echo "Built libjhdf5.2.11.0.dylib ($(nm -gU libjhdf5.2.11.0.dylib | grep -c Java_ncsa) ncsa symbols)."
echo "VALIDATE before use:  ./gradlew test --tests '*HDF5LibraryUnitTest' -Dhdf5Bindings.version=1.2.0-hdf5_2.11.0-arm64"
echo "Currently 10/12 fail (version mismatch) -> do NOT install as default."
