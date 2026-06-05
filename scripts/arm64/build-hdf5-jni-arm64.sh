#!/usr/bin/env bash
# Rebuild the jhdf5 JNI (libjhdf5.2.11.0.dylib) for arm64 and install hdf5-java-bindings:
# 1.2.0-hdf5_2.11.0-arm64 to mavenLocal. Consumed by GATK via -DuseArm64Natives (or
# -Dhdf5Bindings.version=1.2.0-hdf5_2.11.0-arm64).
#
# The upstream jar's Java (ncsa.hdf.hdf5lib) was built against HDF5 1.8.14, so we build the JNI
# against HDF5 1.8.14 (NOT 1.10 -- a 1.10 native fails at runtime with "Bad value"/"Inappropriate
# type"). Recipe:
#   - build HDF5 1.8.14 static C lib for arm64,
#   - take the HDF5 1.10.11 Java JNI .c (only modern copy available), rename hdf.hdf5lib->ncsa,
#   - add HD* macros + map the 1.10 *2 object/ref APIs to their 1.8 forms,
#   - keep ONLY the JNI functions the x86 jar exports (drop 1.10-only extras),
#   - add a variable-length-string H5DwriteString (GATK writes VL strings; the only jar function
#     the 1.10 source doesn't provide),
#   - compile against 1.8.14 + link statically (self-contained dylib).
# Validated: HDF5LibraryUnitTest 12/12 + HDF5SimpleCountCollectionUnitTest pass on arm64.
#
# Requirements: JAVA_HOME=arm64 JDK17, Xcode CLT, git, curl. Run GATK once so the upstream
# hdf5-java-bindings jar is in the Gradle cache (used as the symbol allow-list).
set -euo pipefail
: "${JAVA_HOME:?set JAVA_HOME to an arm64 JDK 17}"; export PATH="$JAVA_HOME/bin:$PATH"
WORK="${WORK:-$HOME/gatk-arm64-natives}"; mkdir -p "$WORK"; cd "$WORK"
ARMVER=1.2.0-hdf5_2.11.0-arm64

# 1) HDF5 1.8.14 static lib (arm64)
[ -f hdf5-1.8.14/src/.libs/libhdf5.a ] || {
  [ -d hdf5-1.8.14 ] || { curl -fsSL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz | tar xz; }
  cp hdf5-1.10.11/bin/config.guess hdf5-1.8.14/bin/ 2>/dev/null || true   # arm64-aware (from 1.10.11 tree)
  cp hdf5-1.10.11/bin/config.sub   hdf5-1.8.14/bin/ 2>/dev/null || true
  ( cd hdf5-1.8.14 && CFLAGS="-O2 -Wno-implicit-function-declaration -Wno-error=implicit-int -Wno-int-conversion" \
      ./configure --enable-static --disable-shared --disable-fortran --disable-cxx --disable-hl --disable-tests --disable-tools >/dev/null \
      && make -C src -j4 >/dev/null )
}
[ -d hdf5-1.10.11/java/src/jni ] || { curl -fsSL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.11/src/hdf5-1.10.11.tar.gz | tar xz; }

# 2) JNI sources: 1.10.11 renamed to ncsa
rm -rf hdf5jni-arm64 && mkdir hdf5jni-arm64
cp hdf5-1.10.11/java/src/jni/*.c hdf5-1.10.11/java/src/jni/*.h hdf5jni-arm64/
cd hdf5jni-arm64
perl -0pi -e 's/Java_hdf_hdf5lib/Java_ncsa_hdf_hdf5lib/g; s{"hdf/hdf5lib}{"ncsa/hdf/hdf5lib}g;' *.c *.h
# legacy constants dropped from 1.10.11
perl -0pi -e 's/(\nH5_GCC_DIAG_ON\("missing-prototypes"\))/\nJNIEXPORT jint JNICALL Java_ncsa_hdf_hdf5lib_HDF5Constants_H5D_1CHUNK_1BTREE(JNIEnv *e, jclass c){return 0;}\nJNIEXPORT jint JNICALL Java_ncsa_hdf_hdf5lib_HDF5Constants_H5F_1ACC_1DEBUG(JNIEnv *e, jclass c){return 0;}\nJNIEXPORT jint JNICALL Java_ncsa_hdf_hdf5lib_HDF5Constants_H5F_1SCOPE_1DOWN(JNIEnv *e, jclass c){return 2;}\n$1/' h5Constants.c

# 3) compat header: HD* macros + 1.10 *2 APIs -> 1.8 forms
cat > hdcompat.h <<'CC'
#ifndef HDCOMPAT_H
#define HDCOMPAT_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define HDprintf printf
#define HDsscanf sscanf
#define HDsprintf sprintf
#define HDsnprintf snprintf
#define HDstrndup strndup
#ifndef H5O_INFO_ALL
#define H5O_INFO_ALL 0
#endif
#define H5Oget_info2(a,b,c)                  H5Oget_info((a),(b))
#define H5Oget_info_by_name2(a,b,c,d,e)      H5Oget_info_by_name((a),(b),(c),(e))
#define H5Oget_info_by_idx2(a,b,c,d,e,f,g,h) H5Oget_info_by_idx((a),(b),(c),(d),(e),(f),(h))
#define H5Ovisit2(a,b,c,d,e,f)               H5Ovisit((a),(b),(c),(d),(e))
#define H5Ovisit_by_name2(a,b,c,d,e,f,g,h)   H5Ovisit_by_name((a),(b),(c),(d),(e),(f),(h))
#define H5Rdereference2(a,b,c,d)             H5Rdereference((a),(c),(d))
#endif
CC

# 4) keep only the functions the upstream x86 jar exports (symbol allow-list)
XJAR=$(find "$HOME/.gradle/caches" -name "hdf5-java-bindings-1.2.0-hdf5_2.11.0.jar" | grep -vE "sources|javadoc|arm64" | head -1)
rm -rf /tmp/xlib && mkdir /tmp/xlib && ( cd /tmp/xlib && unzip -oq "$XJAR" org/broadinstitute/hdf5/libjhdf5.2.11.0.dylib )
nm -gU /tmp/xlib/org/broadinstitute/hdf5/libjhdf5.2.11.0.dylib | grep -oE "Java_[A-Za-z0-9_]+" | sort -u > /tmp/target_syms.txt
python3 - <<'PY'
import re, pathlib
target=set(l.strip() for l in open('/tmp/target_syms.txt'))
for cf in pathlib.Path('.').glob('*.c'):
    s=cf.read_text(); out=[]; last=0
    pat=re.compile(r'JNIEXPORT\b[^\n]*\bJNICALL\s*\n?\s*(Java_[A-Za-z0-9_]+)\s*\(', re.M)
    while True:
        m=pat.search(s,last)
        if not m: out.append(s[last:]); break
        name=m.group(1); ds=m.start(); p=s.find(')',m.end()-1); b=s.find('{',p)
        d=0; j=b
        while j<len(s):
            if s[j]=='{':d+=1
            elif s[j]=='}':
                d-=1
                if d==0:break
            j+=1
        e=j+1
        out.append(s[last:e] if name in target else s[last:ds]); last=e
    cf.write_text(''.join(out))
PY

# 5) variable-length-string H5DwriteString (GATK writes VL strings)
cat >> h5dImp.c <<'CC'

JNIEXPORT jint JNICALL
Java_ncsa_hdf_hdf5lib_H5_H5DwriteString(JNIEnv *env, jclass clss, jint dataset_id, jint mem_type_id,
                                        jint mem_space_id, jint file_space_id, jint xfer_plist_id,
                                        jobjectArray j_buf)
{
    herr_t status = FAIL; jsize n = 0, i; char **wdata = NULL; UNUSED(clss);
    if (NULL == j_buf) H5_NULL_ARGUMENT_ERROR(ENVONLY, "H5DwriteString: buffer is NULL");
    if ((n = ENVPTR->GetArrayLength(ENVONLY, j_buf)) <= 0) { CHECK_JNI_EXCEPTION(ENVONLY, JNI_TRUE);
        H5_BAD_ARGUMENT_ERROR(ENVONLY, "H5DwriteString: length <= 0"); }
    if (NULL == (wdata = (char **)calloc((size_t)n, sizeof(char *)))) H5_OUT_OF_MEMORY_ERROR(ENVONLY, "H5DwriteString: oom");
    for (i = 0; i < n; i++) {
        jstring obj = (jstring)ENVPTR->GetObjectArrayElement(ENVONLY, (jobjectArray)j_buf, i);
        if (NULL != obj) { const char *utf8 = NULL;
            PIN_JAVA_STRING(ENVONLY, obj, utf8, NULL, "H5DwriteString: not pinned");
            if (NULL != utf8) wdata[i] = strdup(utf8);
            UNPIN_JAVA_STRING(ENVONLY, obj, utf8); ENVPTR->DeleteLocalRef(ENVONLY, obj); }
    }
    if ((status = H5Dwrite((hid_t)dataset_id,(hid_t)mem_type_id,(hid_t)mem_space_id,(hid_t)file_space_id,(hid_t)xfer_plist_id, wdata)) < 0)
        H5_LIBRARY_ERROR(ENVONLY);
done:
    if (wdata) { for (i = 0; i < n; i++) if (wdata[i]) free(wdata[i]); free(wdata); }
    return (jint)status;
}
CC

# 6) compile + link (static HDF5 1.8.14 -> self-contained)
H8SRC="$WORK/hdf5-1.8.14/src"
for c in *.c; do clang -c -O2 -fPIC -include hdcompat.h -I. -I"$H8SRC" -I"$JAVA_HOME/include" -I"$JAVA_HOME/include/darwin" "$c" -o "${c%.c}.o"; done
clang -dynamiclib -install_name libjhdf5.2.11.0.dylib -o libjhdf5.2.11.0.dylib *.o "$H8SRC/.libs/libhdf5.a" -lz -lm
file libjhdf5.2.11.0.dylib | grep -q arm64

# 7) package + install to mavenLocal
DEST="$HOME/.m2/repository/org/broadinstitute/hdf5-java-bindings/$ARMVER"; mkdir -p "$DEST"
cp "$XJAR" "/tmp/hdf5-$ARMVER.jar"
rm -rf /tmp/hp && mkdir -p /tmp/hp/org/broadinstitute/hdf5
cp libjhdf5.2.11.0.dylib /tmp/hp/org/broadinstitute/hdf5/
( cd /tmp/hp && jar uf "/tmp/hdf5-$ARMVER.jar" org/broadinstitute/hdf5/libjhdf5.2.11.0.dylib )
cp "/tmp/hdf5-$ARMVER.jar" "$DEST/hdf5-java-bindings-$ARMVER.jar"
curl -fsSL "https://repo1.maven.org/maven2/org/broadinstitute/hdf5-java-bindings/1.2.0-hdf5_2.11.0/hdf5-java-bindings-1.2.0-hdf5_2.11.0.pom" \
  | sed "s#<version>1.2.0-hdf5_2.11.0</version>#<version>$ARMVER</version>#" > "$DEST/hdf5-java-bindings-$ARMVER.pom"
echo "Installed org.broadinstitute:hdf5-java-bindings:$ARMVER (HDF5 1.8.14, arm64) to mavenLocal"
