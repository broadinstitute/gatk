# Running GATK natively on Apple Silicon (macOS arm64)

> Status: **in progress** (branch `arm64-native-port`). This document tracks the work to make
> GATK run natively on Apple Silicon (aarch64) instead of under x86_64 Rosetta emulation.

Today the [main README](../../README.md) instructs ARM-Mac users to install the **x86_64**
miniconda and rely on Rosetta. The goal of this effort is to remove that requirement: a
native arm64 JVM, native osx-arm64 Conda environment, NEON-accelerated native libraries, and
optional Apple-GPU (Metal) acceleration for the deep-learning tools.

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
