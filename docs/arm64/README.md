# Running GATK natively on Apple Silicon (macOS arm64)

> Status: **in progress** (branch `arm64-native-port`). This document tracks the work to make
> GATK run natively on Apple Silicon (aarch64) instead of under x86_64 Rosetta emulation.

Today the [main README](../../README.md) instructs ARM-Mac users to install the **x86_64**
miniconda and rely on Rosetta. The goal of this effort is to remove that requirement: a
native arm64 JVM, native osx-arm64 Conda environment, NEON-accelerated native libraries, and
optional Apple-GPU (Metal) acceleration for the deep-learning tools.

## Correctness principle: bit-equivalence is the acceptance gate

Native arm64 results **must be bit-identical to the current (x86) version**. This gives one
decision rule for every accelerated path:

> Keep a native acceleration (NEON, GPU/MPS) **only if it is bit-identical** to the reference
> path. Otherwise use the deterministic path — the pure-Java PairHMM / Smith-Waterman, which
> run natively on arm64 and are architecture-independent, hence inherently bit-equivalent.

Implications:
- **Deterministic / integer / exact libraries** (bwa-mem, fermi-lite, MUMmer, GenomicsDB, HDF5):
  bit-equivalence is achievable — same algorithms; `sse2neon` maps the integer SSE ops to NEON.
  Validated by diffing arm64 output against x86 output on identical input.
- **Floating-point PairHMM** (Intel GKL AVX): a NEON/SIMDe build is unlikely to be bit-identical
  to x86 AVX (FMA/reduction ordering). NEON is enabled by default for SPEED; the bit-equivalent
  pure-Java path is available via `-DgklVersion=0.9.1`.
- **Compression** (deflate): compressed *bytes* differ per codec, but decompressed content is
  identical; bit-equivalence is defined on logical content (reads/variants), not the gzip stream.
- **GPU/MPS**: enabled only after confirming its output is bio-identical to the CPU output.

### IMPORTANT empirical correction: bit-identical floating-point scores are NOT fully achievable arm64↔x86

Measured on real data (`RampedHaplotypeCallerIntegrationTest`, which compares raw likelihood
matrices **exactly, no tolerance**): the arm64 output differs from the committed x86 expected at
the **last ULP** of a double — e.g. `-3.450078543453124` (x86) vs `-3.4500785434531234` (arm64),
~6e-16. **This happens even with the pure-Java PairHMM** (`-DgklVersion=0.9.1`), so it is **not**
a NEON artifact. Root cause: Java's `strictfp` only governs `+ - * /`; the transcendental
`Math.log/exp/pow` GATK uses in its likelihood code are **not** guaranteed identical across CPU
architectures (only `StrictMath.*` is, and GATK uses `Math.*`). So:
- **Bio-identical: YES** — these ~1e-16 differences never change a variant call or filter. All CNV
  tools (CollectReadCounts 6/6, DenoiseReadCounts 48/48, CreateReadCountPanelOfNormals 36/36) pass
  exactly on arm64; HaplotypeCaller VCF-level tests pass; only the internal exact-likelihood-dump
  ramp tests (2) show the ULP difference.
- **Bit-identical to x86: NOT achievable** for the FP caller on arm64 (independent of NEON vs Java)
  without rewriting GATK's hot paths to `StrictMath` — which would also change x86 output, i.e.
  redefine the reference. The realistic, correct bar for cross-architecture is **bio-identical**.

## Prerequisites (macOS Apple Silicon)

- A **native arm64 JDK 17** (GATK/Gradle require Java 17; the build fails under Java 18+).
  e.g. `brew install openjdk@17`, then point Gradle at it:
  `export JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home`.
  Confirm `java -XshowSettings:properties 2>&1 | grep os.arch` → `aarch64`.

## Status

- [x] **Track A — foundation** (aarch64 detection, build scaffolding, snappy, netlib):
      `./gradlew installDist` builds and `NativeUtilsUnitTest` passes on native arm64 JDK 17
      (`os.arch=aarch64`, `sysctl.proc_translated=0`). End-to-end verified: `gatk PrintReads`
      produced a valid BGZF BAM natively (pure-Java `JdkDeflater` path, no Rosetta).
- [~] **Track B — native library rebuilds** (scripts in `scripts/arm64/`, jars in mavenLocal,
      enabled with `-DuseArm64Natives=true`):
  - [x] **MUMmer** — bundled arm64 binaries, bit-identical (`MummerExecutorUnitTest` 4/4).
  - [x] **gatk-bwamem-jni** `1.0.4-arm64` (sse2neon) — bit-identical alignments (7/7).
  - [x] **gatk-fermilite-jni** `1.2.0-arm64` (sse2neon) — bit-identical assembly (2/2).
  - [x] **GKL** `0.9.1-arm64` — NEON PairHMM (SIMDe). **Default = speed** (per user); matches x86
        within 1e-5 (NOT bit-identical). Bit-identical Java path via `-DgklVersion=0.9.1`.
  - [x] **GenomicsDB** `1.5.5` — already a universal x86_64+arm64 binary; runs natively, no rebuild.
  - [x] **hdf5-java-bindings** `1.2.0-hdf5_2.11.0-arm64` — jhdf5 rebuilt against HDF5 1.8.14 (+ VL
        `H5DwriteString`); `HDF5LibraryUnitTest` 12/12 + `HDF5SimpleCountCollectionUnitTest` pass.
        CNV/HDF5 tools now native. **All native deps done.**
  Full native arm64 GATK assembles with `installDist -DuseArm64Natives=true`.

## Validation results (real LFS data, `-DuseArm64Natives`)

- **CNV + HaplotypeCaller packages: 8775 / 8786 tests pass** (99.87%). Zero hard native errors
  (the `libgkl_compression` "incompatible architecture" lines are benign — it falls back to
  `JdkDeflater`).
- **All 312 GATK tools load on arm64** (`./gatk <Tool> --help`). The ~120 that return a non-zero
  exit are Picard-style tools whose `--help` exits non-zero *by convention* — they print usage
  and load fine.
- **Speed**: HaplotypeCaller with NEON PairHMM vs pure-Java — ~**1.3–1.4× faster** (WES chr20:
  32.9s vs 43.4s; sparse-WGS 3 Mb: 18.3s vs 26.4s). NEON and Java call the **identical 2977
  variants** (CHROM/POS/REF/ALT/GT, 0 diffs) → bio-identical. (Still far from DRAGEN/FPGA — expected.)
- **Python (osx-arm64 conda)**: env builds (pysam 0.24.0, pyvcf3, OpenBLAS, cpu/MPS PyTorch);
  `gcnvkernel` imports and **PyMC NUTS sampling runs**. One fix required: PyTensor's bundled clang
  fails to link on arm64, so the env sets `PYTENSOR_FLAGS=cxx=/usr/bin/clang++` (system Xcode
  clang) — without it GermlineCNVCaller/DetermineGermlineContigPloidy fail to compile.

### The 11 failures — `ModelSegments` (MCMC), an inherent cross-arch limit (not a port bug)

`ModelSegments` is MCMC/Gibbs-based; its posterior means drift ~1% on arm64 vs the x86 expected
(e.g. `0.303401` vs `0.306041`, well above the 1e-5 test tolerance). Root cause: the same Java
`Math.log/exp/pow` ~1-ULP cross-architecture difference, **amplified by stochastic sampling**
(the DeepVariant phenomenon). The segmentation is **structurally bio-identical** (same segments);
only posterior estimates differ. Making it match the x86 expected exactly is **not achievable on
arm64** without either (a) converting GATK's copy-number hot paths to `StrictMath` (which would
also change x86 output, i.e. requires regenerating the expected files on a canonical platform), or
(b) running these steps under x86 Rosetta. This is a property of Java FP + MCMC, not of the port.
- [~] Track C — osx-arm64 conda env: arm64 template added; gradle auto-selects it and generates
      the arm64 yml. *Not yet created/run end-to-end (requires a full conda solve).*
- [~] Track D — PyTorch MPS: device auto-selection (CUDA→MPS→CPU) added to nvscorevariants.py
      and syntax-checked. *Not yet run (requires the conda env + model).*

## What blocks native arm64 today

| Dependency | Kind | arm64 status | Has Java fallback? |
|---|---|---|---|
| `com.intel.gkl:gkl` | JNI (AVX PairHMM / SmithWaterman, ISA-L deflate) | x86-only | yes (slower) |
| `org.broadinstitute:gatk-bwamem-jni` | JNI (bwa-mem, SSE2) | x86-only | **no** |
| `org.broadinstitute:gatk-fermilite-jni` | JNI (fermi-lite, SSE) | x86-only | **no** |
| `org.genomicsdb:genomicsdb` | JNI (C++/TileDB) | check latest | **no** |
| `org.broadinstitute:hdf5-java-bindings` | JNI (jhdf5) | x86-only | **no** |
| bundled MUMmer 4.0.0rc1 binaries | native exes | x86-only | **no** |
| `com.github.fommil.netlib:netlib-native_*-x86_64` | JNI (BLAS) | x86-only | yes (+ Accelerate) |
| `org.xerial.snappy:snappy-java:1.1.10.4` | JNI | bump to 1.1.11.x | yes |
| Conda env (`blas=mkl`, `pytorch=*mkl*`) | Python | MKL is Intel-only | n/a |

## Strategy

1. **In-repo foundation (this branch):** add aarch64 architecture detection, scaffold the
   build to optionally consume locally-rebuilt `-arm64` native jars (defaulting to the
   upstream x86 jars so the build stays green and Rosetta keeps working), bump snappy, drop
   the useless x86 netlib natives on arm64 (rely on Apple **Accelerate** via
   `dev.ludovic.netlib` from Spark 3.5), and produce a native osx-arm64 Conda environment.
2. **Native library rebuilds (external forks → `mavenLocal()`):** rebuild each JNI library
   for arm64 and publish a fat jar (arm64 dylib added next to the x86 ones) as `<ver>-arm64`.
   See [native-libs.md](native-libs.md).
3. **GPU:** see [gpu.md](gpu.md). On Apple Silicon there is **no CUDA** — the realistic GPU
   win is the **PyTorch MPS** backend for the deep-learning tools (`NVScoreVariants`,
   `TrainVariantAnnotationsModel`/`ScoreVariantAnnotations`). GPU acceleration of the core
   C/SIMD kernels (PairHMM, Smith-Waterman, bwa) would require writing Metal compute shaders
   from scratch — tracked as a research stretch goal, not part of this branch.

## Consuming locally-rebuilt arm64 native jars

`mavenLocal()` is already a repository in `build.gradle`. Once a native lib has been rebuilt
and `publishToMavenLocal`'d under its `-arm64` version, enable it for the GATK build with:

```
./gradlew installDist -DuseArm64Natives=true
```

Without that flag the build uses the upstream x86 artifacts (so the project always builds,
even before any native lib has been ported).

## Verifying you are truly native (not Rosetta)

```
sysctl sysctl.proc_translated          # 0 (or absent) = native; 1 = running under Rosetta
java -XshowSettings:properties 2>&1 | grep os.arch    # aarch64
file <extracted-temp-dylib>            # "Mach-O 64-bit ... arm64"
```
Activity Monitor's **Kind** column shows **Apple** for native processes, **Intel** under Rosetta.
