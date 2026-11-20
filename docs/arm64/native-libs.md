# Rebuilding GATK's native libraries for arm64 (Apple Silicon)

Each library below is built from a forked upstream repo, compiled on an Apple Silicon mac,
and published to `mavenLocal()` as a fat jar (the new arm64 dylib added alongside the existing
x86 dylibs) under a `-arm64` version. GATK then selects those versions via
`-DuseArm64Natives=true` (see [README](README.md)). Because the GKL/bwa/fermi jars locate
their native by `os.arch` at runtime, no GATK Java change is needed once the jar contains an
arm64 dylib.

Pin exact upstream tags + the SIMDe / sse2neon / HDF5 / protobuf versions you build against so
the `-arm64` jars are reproducible.

## Intel GKL — `com.intel.gkl:gkl` (Intel-HLS/GKL, CMake + JNI) — empirical NEON build exists; default is Java

**The NEON PairHMM was built and measured** (`scripts/arm64/build-gkl-arm64.sh`): GKL's AVX
(256-bit) PairHMM + `gkl_utils` are lowered to NEON via SIMDe and installed as `gkl:0.9.1-arm64`
in mavenLocal. It loads and runs on Apple Silicon and matches the x86 reference within GATK's
PairHMM tolerance of **1e-5** (`VectorPairHMMUnitTest`). **It is NOT bit-identical** to the
pure-Java PairHMM — exactly like the x86 AVX kernel, a vectorized float kernel differs from the
Java double computation in the low bits.

So there are two mutually-exclusive modes for the score path:
- **Bit-identical scores (default):** no arm64 GKL → GATK falls back to the pure-Java PairHMM,
  which is bit-identical across architectures (Java 17 strictfp; `PairHMMUnitTest` 25928/25928 on
  arm64). Slower (no SIMD).
- **Native speed (opt-in):** `-DgklVersion=0.9.1-arm64` (+ `-pairHMM AVX_LOGLESS_CACHING`) uses the
  NEON PairHMM. Faster, but results match x86 only within ~1e-5 (NOT bit-identical).

The rest of this section explains why the other GKL kernels are not pursued:

- **PairHMM (AVX2/AVX512, floating point)** drives the genotype likelihoods/scores. A NEON/SIMDe
  build will not be bit-for-bit identical to x86 AVX (FMA/reduction ordering), so it could not be
  used anyway. The **pure-Java PairHMM** is bit-identical across architectures (Java 17 is
  strictfp) — *verified: `PairHMMUnitTest` passes 25928/25928 on arm64*. GATK's default
  `FASTEST_AVAILABLE` automatically falls back to it when GKL's x86 native fails to load on arm64.
- **Compression (ISA-L)** is hand-written x86 assembly (won't port) and produces different
  compressed *bytes* per codec anyway; decompressed content is identical, and the JDK
  deflater/inflater (`--use-jdk-deflater/--use-jdk-inflater`, or the automatic fallback) is the
  bit-identical-content path.
- **SmithWaterman (AVX2/AVX512, integer)** is the *only* kernel that could be both bit-identical
  and faster on NEON (via SIMDe). It is tracked as a **future optional optimization** — build just
  that kernel with SIMDe, gate on bit-exactness vs the Java aligner. Not required for correctness
  (the Java aligner is bit-identical and already used as the fallback).

So on aarch64 `gklVersion` stays `0.9.1` (the x86 jar; its native simply fails to load → Java
fallback). **Do not force AVX on aarch64** (`-pairHMM AVX_LOGLESS_CACHING`,
`--smith-waterman AVX_ENABLED`) — those bypass the fallback and will error.

Verify: `HaplotypeCaller` on arm64 with defaults runs via the Java PairHMM (bit-identical scores);
`PairHMMUnitTest` 25928/25928 on arm64.

## bwa-mem JNI — `org.broadinstitute:gatk-bwamem-jni` & fermi-lite JNI — `gatk-fermilite-jni`
No Java fallback — these must produce working arm64 dylibs. Use **`sse2neon.h`** (define
`__SSE2__`, replace `<emmintrin.h>`/`<xmmintrin.h>` includes, `-march=armv8-a+simd`) — the same
approach Bioconda uses to build bwa/fermi-lite for osx-arm64. Patch the Makefile arch branch
and the JNI loader to recognize `aarch64`.

Verify: `BwaMemIndexImageCreator` on a small fasta; an SV smoke test (`FindBreakpointEvidenceSpark`)
or `FilterAlignmentArtifacts` for fermi-lite.

## GenomicsDB — `org.genomicsdb:genomicsdb:1.5.5` — ALREADY NATIVE (no action)
The published `1.5.5` jar bundles `libtiledbgenomicsdb.dylib` as a **universal binary
(x86_64 + arm64)** (`lipo -archs` → `x86_64 arm64`). So GenomicsDB runs natively on Apple
Silicon with the existing dependency — no rebuild, no version change.

Verified on arm64: `GenomicsDBIntegrationTest` — the native lib loads and `GenomicsDBImport`
runs; the 2 LFS-independent tests pass. (The other 5 fail only because the git-lfs reference
`human_g1k_v37.20.21.fasta` wasn't pulled — environmental, not a data difference.)

## HDF5 jhdf5 — `org.broadinstitute:hdf5-java-bindings` — DONE (`1.2.0-hdf5_2.11.0-arm64`)

**Fixed and validated** (`scripts/arm64/build-hdf5-jni-arm64.sh`). The key was the HDF5 version:
the upstream jar's `ncsa.hdf.hdf5lib` Java was built against **HDF5 1.8.14** (read from the x86
dylib: `strings ... | grep "HDF5 Version"` → `1.8.14`). Building against 1.10.x produced runtime
"Bad value"/"Inappropriate type" errors (the 1.8→1.10 API/enum drift). The working build:

1. compile **HDF5 1.8.14** as a static arm64 C lib;
2. take the HDF5 1.10.11 Java JNI `.c` (only modern copy), rename `hdf.hdf5lib`→`ncsa.hdf.hdf5lib`;
3. add HD* macros and map the 1.10 `*2` object/reference APIs to their 1.8 forms;
4. keep ONLY the JNI functions the x86 jar exports (drop 1.10-only extras);
5. add a **variable-length-string `H5DwriteString`** (GATK writes VL strings via
   `H5Tset_size(type, H5T_VARIABLE)`; it's the one jar function the 1.10 source lacks — a flat
   buffer segfaults in `H5T_convert`, so it must build a `char**` pointer array);
6. compile against 1.8.14 and link statically → self-contained `libjhdf5.2.11.0.dylib` (arm64).

Validated on arm64: **`HDF5LibraryUnitTest` 12/12** and **`HDF5SimpleCountCollectionUnitTest`**
pass (also via `-DuseArm64Natives`). CNV/HDF5 tools (`CollectReadCounts`, `DenoiseReadCounts`,
`CreateReadCountPanelOfNormals`, germline-CNV read-count I/O) now run natively on Apple Silicon.


This is the one native lib not yet rebuilt. The bindings jar bundles a **prebuilt** x86
`libjhdf5.2.11.0.dylib` (SIS/CISD jhdf5 2.11.0) — the repo has no C sources to recompile, so an
arm64 build means building **SIS jhdf5 2.11.0 from upstream** against an arm64 HDF5 C library and
repackaging. Risk: jhdf5 2.11.0 targets an HDF5 1.10/1.12-era C API, while Homebrew's current
`hdf5` is a much newer 2.x — likely API/ABI mismatch, so a matching older HDF5 (1.10.x/1.12.x)
must be built for arm64 first. Held at x86 for now (`hdf5Bindings.version` default), so on arm64
the native fails to load and the **CNV/HDF5 tools** (`CollectReadCounts`, `DenoiseReadCounts`,
`CreateReadCountPanelOfNormals`, GermlineCNV read-count I/O) do not run natively yet.

### Groundwork done + precise blocker
- The x86 `libjhdf5.2.11.0.dylib` exports **1041 `Java_ncsa_hdf_hdf5lib_*` JNI symbols** with HDF5
  **statically embedded** (`otool -L` shows no external HDF5 dep). Namespace is the pre-1.10
  `ncsa.hdf.hdf5lib` (classes live in the bundled `jarhdf5-2.11.0.jar`).
- An **arm64 HDF5 1.10 is available**: `brew install hdf5@1.10` →
  `/opt/homebrew/opt/hdf5@1.10/lib/libhdf5.dylib` (arm64).
- The HDF5 **1.10.11 source ships the JNI implementations** (`java/src/jni/*.c`, **1253**
  `Java_hdf_hdf5lib_*` functions — a *superset* of the 1041 needed).

### Why it's not a quick swap (the real risk)
The 1.10.11 JNI uses the **`hdf.hdf5lib`** namespace, not `ncsa.hdf.hdf5lib`, and is a newer API
than CISD jhdf5 2.11.0. A build therefore needs (a) renaming `hdf_hdf5lib`→`ncsa_hdf_hdf5lib` in
the JNI function names, and (b) **matching each function's JNI signature** to what
`jarhdf5-2.11.0.jar` declares — any drift surfaces only at runtime (UnsatisfiedLinkError or a
crash) when a CNV tool calls that specific `H5*` function. So it must be built signature-matched
(generate headers via `javac -h` from the jar's `ncsa.hdf.hdf5lib` classes, map to the 1.10.11
`.c` bodies) and **validated** before shipping.

### Recipe (to complete + validate)
1. `brew install hdf5@1.10` (done).
2. `javac -h` the `ncsa.hdf.hdf5lib.*` classes from `jarhdf5-2.11.0.jar` → exact JNI headers.
3. Take HDF5 1.10.11 `java/src/jni/*.c`, rename `hdf_hdf5lib`→`ncsa_hdf_hdf5lib`, compile against
   `hdf5@1.10` + the generated headers → `libjhdf5.2.11.0.dylib` (arm64, statically linking
   `libhdf5.a`). Resolve any missing/mismatched signatures against the jar.
4. Replace `org/broadinstitute/hdf5/libjhdf5.2.11.0.dylib` in the jar → install as
   `1.2.0-hdf5_2.11.0-arm64` → `-Dhdf5Bindings.version=1.2.0-hdf5_2.11.0-arm64`.
5. **Validate:** `CollectReadCounts` (writes `.hdf5`) → `DenoiseReadCounts`; diff vs x86 output.

Robust alternative: build the **CISD jhdf5 2.11.0** source directly (it produces exactly these
`ncsa` symbols), then repackage — heavier build but no signature-matching guesswork.

## MUMmer 4.0.0rc1 (bundled in GATK, not a jar)
Compile from source on arm64-mac, zip `nucmer`/`delta-filter`/`show-snps` with the same layout
as the existing x86 zip, and drop it at
`src/main/resources/org/broadinstitute/hellbender/utils/alignment/MUMmer-4.0.0rc1_mac_aarch64.zip`.
`MummerExecutor` already selects it on arm64 once present.

Verify: `CompareReferences --mode FULL_ALIGNMENT` on two small references.
