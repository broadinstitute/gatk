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

## HDF5 jhdf5 — `org.broadinstitute:hdf5-java-bindings` — REMAINING GAP (x86-only)
This is the one native lib not yet rebuilt. The bindings jar bundles a **prebuilt** x86
`libjhdf5.2.11.0.dylib` (SIS/CISD jhdf5 2.11.0) — the repo has no C sources to recompile, so an
arm64 build means building **SIS jhdf5 2.11.0 from upstream** against an arm64 HDF5 C library and
repackaging. Risk: jhdf5 2.11.0 targets an HDF5 1.10/1.12-era C API, while Homebrew's current
`hdf5` is a much newer 2.x — likely API/ABI mismatch, so a matching older HDF5 (1.10.x/1.12.x)
must be built for arm64 first. Held at x86 for now (`hdf5Bindings.version` default), so on arm64
the native fails to load and the **CNV/HDF5 tools** (`CollectReadCounts`, `DenoiseReadCounts`,
`CreateReadCountPanelOfNormals`, GermlineCNV read-count I/O) do not run natively yet.

Recipe (future): build HDF5 1.10.x/1.12.x for arm64 → build SIS jhdf5 2.11.0 JNI against it →
replace `org/broadinstitute/hdf5/libjhdf5.2.11.0.dylib` in the jar → install as
`1.2.0-hdf5_2.11.0-arm64` → `-Dhdf5Bindings.version=1.2.0-hdf5_2.11.0-arm64`.
Verify: `CollectReadCounts` (writes `.hdf5`) → `DenoiseReadCounts`.

## MUMmer 4.0.0rc1 (bundled in GATK, not a jar)
Compile from source on arm64-mac, zip `nucmer`/`delta-filter`/`show-snps` with the same layout
as the existing x86 zip, and drop it at
`src/main/resources/org/broadinstitute/hellbender/utils/alignment/MUMmer-4.0.0rc1_mac_aarch64.zip`.
`MummerExecutor` already selects it on arm64 once present.

Verify: `CompareReferences --mode FULL_ALIGNMENT` on two small references.
