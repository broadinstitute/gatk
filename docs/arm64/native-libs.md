# Rebuilding GATK's native libraries for arm64 (Apple Silicon)

Each library below is built from a forked upstream repo, compiled on an Apple Silicon mac,
and published to `mavenLocal()` as a fat jar (the new arm64 dylib added alongside the existing
x86 dylibs) under a `-arm64` version. GATK then selects those versions via
`-DuseArm64Natives=true` (see [README](README.md)). Because the GKL/bwa/fermi jars locate
their native by `os.arch` at runtime, no GATK Java change is needed once the jar contains an
arm64 dylib.

Pin exact upstream tags + the SIMDe / sse2neon / HDF5 / protobuf versions you build against so
the `-arm64` jars are reproducible.

## Intel GKL — `com.intel.gkl:gkl` (Intel-HLS/GKL, CMake + JNI)
Three independent kernels:
- **PairHMM (AVX2/AVX512)** and **SmithWaterman (AVX2/AVX512)**: compile the intrinsic kernels
  with **[SIMDe](https://github.com/simd-everywhere/simde)** (`simde/x86/avx2.h`,
  `simde/x86/avx512.h`, `-DSIMDE_ENABLE_NATIVE_ALIASES`), which lowers AVX intrinsics to NEON.
  Drop `-mavx*`, add `-march=armv8-a+simd`. SmithWaterman (integer) is low numeric risk;
  PairHMM (FP, log-space) is the numerics watch-item, especially the AVX512→SIMDe path.
- **Compression (ISA-L, x86 asm)**: use ISA-L's `aarch64/` backend or **zlib-ng** (ARM NEON).
  Lowest risk overall — GKL's `IntelDeflaterFactory`/`IntelInflaterFactory` already fall back
  to the JDK, and GATK exposes `--use-jdk-deflater/--use-jdk-inflater`.
- Confirm GKL's `NativeLibraryLoader` maps `os.arch=aarch64` → the arm64 dylibs.

Verify: `HaplotypeCaller -pairHMM AVX_LOGLESS_CACHING` (forces native, throws if unavailable);
`--smith-waterman AVX_ENABLED`; `PrintReads` without `--use-jdk-*` and confirm the startup log
prints `Deflater: IntelDeflater`.

## bwa-mem JNI — `org.broadinstitute:gatk-bwamem-jni` & fermi-lite JNI — `gatk-fermilite-jni`
No Java fallback — these must produce working arm64 dylibs. Use **`sse2neon.h`** (define
`__SSE2__`, replace `<emmintrin.h>`/`<xmmintrin.h>` includes, `-march=armv8-a+simd`) — the same
approach Bioconda uses to build bwa/fermi-lite for osx-arm64. Patch the Makefile arch branch
and the JNI loader to recognize `aarch64`.

Verify: `BwaMemIndexImageCreator` on a small fasta; an SV smoke test (`FindBreakpointEvidenceSpark`)
or `FilterAlignmentArtifacts` for fermi-lite.

## GenomicsDB — `org.genomicsdb:genomicsdb`
**Try a version bump first** — recent releases ship osx-arm64 natives. Only build from source
(CMake; htslib/protobuf/openssl arm64; protobuf ABI must match `protobuf-java:3.25.5`) if no
released jar has the arm64 native.

Verify: `GenomicsDBImport` of small GVCFs → `GenotypeGVCFs`.

## HDF5 jhdf5 — `org.broadinstitute:hdf5-java-bindings`
Rebuild `libjhdf5.dylib` + HDF5 C lib for arm64 against a brew/conda arm64 HDF5; pin the HDF5
version to the `_2.11.0` suffix.

Verify: `CollectReadCounts` (writes `.hdf5`) → `DenoiseReadCounts`.

## MUMmer 4.0.0rc1 (bundled in GATK, not a jar)
Compile from source on arm64-mac, zip `nucmer`/`delta-filter`/`show-snps` with the same layout
as the existing x86 zip, and drop it at
`src/main/resources/org/broadinstitute/hellbender/utils/alignment/MUMmer-4.0.0rc1_mac_aarch64.zip`.
`MummerExecutor` already selects it on arm64 once present.

Verify: `CompareReferences --mode FULL_ALIGNMENT` on two small references.
