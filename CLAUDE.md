# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

GATK4 (Genome Analysis Toolkit) — a large Java 17 / Gradle bioinformatics toolkit that unifies GATK and Picard tools under one command-line framework, with optional Apache Spark backends for distributed execution. Built on [htsjdk](https://github.com/samtools/htsjdk) for BAM/CRAM/VCF I/O. Some tools (CNV, variant scoring) have Python components run via a Conda environment.

## Build & run

```bash
./gradlew installDist      # FAST dev build; run tools locally from the clone (use this while iterating)
./gradlew installAll       # dev build, local + Spark cluster
./gradlew bundle           # full build -> build/gatk-VERSION.zip standalone distribution
./gradlew localJar         # build/libs/gatk-package-VERSION-local.jar (usable outside clone)
./gradlew sparkJar         # build/libs/gatk-package-VERSION-spark.jar (no Spark/Hadoop bundled)
./gradlew clean

./gatk --list              # list all tools
./gatk ToolName --help     # help for one tool
./gatk ToolName <args>     # run a tool (Spark tools end in "Spark" and accept "-- --spark-master ..." )
```

Always run tools through the `./gatk` wrapper, never `java -jar` directly — the wrapper sets critical JVM/system properties (e.g. htsjdk compression level). Add `org.gradle.daemon=true` to `~/.gradle/gradle.properties` to avoid gradle startup cost.

## Test

```bash
./gradlew test                                      # non-cloud unit + integration tests (default)
./gradlew test --tests *SomeTestClass               # single class
./gradlew test --tests *SomeTest.someMethod         # single method
./gradlew test --tests all.in.some.package*         # package
./gradlew jacocoTestReport                           # coverage -> build/reports/jacoco/test/html/index.html
```

- Test report: `build/reports/tests/test/index.html`. Framework is **TestNG** (not JUnit).
- The `TEST_TYPE` env var selects a disjoint partition: `unit`, `integration`, `cloud`, `spark`, `conda`, or `all`. See `testsettings.gradle` for the exact include/exclude groups (cloud, bucket, python, R, funcotatorValidation, variantcalling, spark). Cloud tests need `gcloud` auth plus `HELLBENDER_TEST_PROJECT` / `HELLBENDER_TEST_STAGING` / `HELLBENDER_TEST_INPUTS` env vars.
- `TEST_VERBOSITY=minimal` reduces output; `GATK_STACKTRACE_ON_USER_EXCEPTION=true` shows stack traces for `UserException`.

## Test data (git-lfs) — required before running tests

Large test files are stored in git-lfs. Run `git lfs install` then `git lfs pull` (≈5 GB) from the repo root. Build-time large files live in `src/main/resources/large/` (ML models, native libs — packaged into the jar, required to *build*); test-only large files in `src/test/resources/large/` (required to run tests). Any file placed under a `resources/large/` dir is auto-tracked by lfs — never run `git lfs track` manually.

## Architecture

**Entry point & tool discovery.** [Main.java](src/main/java/org/broadinstitute/hellbender/Main.java) is the dispatcher. It discovers tools at runtime by scanning packages (`getPackageList()`) for subclasses of `CommandLineProgram`. To add a tool you generally just write a `CommandLineProgram` subclass in a scanned package and annotate it with `@CommandLineProgramProperties` (assigning a program group from `cmdline/programgroups/`); it then appears in `--list` automatically. Picard tools are wrapped via `PicardCommandLineProgram` / `PicardCommandLineProgramExecutor`.

**The "walker" engine** ([engine/](src/main/java/org/broadinstitute/hellbender/engine/)) is the core abstraction. `GATKTool extends CommandLineProgram` defines the lifecycle: `onTraversalStart()` → `traverse()` → `onTraversalSuccess()` → `closeTool()`. Rather than implement `traverse()` directly, most tools extend a **walker base class** that handles input iteration and hands the tool one unit at a time via an `apply(...)` callback:
- `ReadWalker` — one read at a time
- `LocusWalker` — pileup per genomic locus
- `VariantWalker` / `MultiVariantWalker` — one variant at a time
- `FeatureWalker`, `IntervalWalker`, `AssemblyRegionWalker` (drives the assembly-based callers)

Each walker provides standard contexts: `ReadsContext`, `ReferenceContext`, `FeatureContext` (features come from `FeatureDataSource`/`FeatureManager`; all paths flow through `GATKPath`).

**Spark mirror.** [engine/spark/](src/main/java/org/broadinstitute/hellbender/engine/spark/) re-implements the same traversals for distributed execution: `GATKSparkTool` plus `ReadWalkerSpark`, `LocusWalkerSpark`, `VariantWalkerSpark`, `AssemblyRegionWalkerSpark`, with `SparkSharder` for partitioning. A tool typically has both a local form and a `...Spark` form.

**Pluggable command-line components.** Read filters ([engine/filters/](src/main/java/org/broadinstitute/hellbender/engine/filters/)), read transformers ([transformers/](src/main/java/org/broadinstitute/hellbender/transformers/)), and variant annotations are plugins selectable from the command line via the plugin descriptor framework in [cmdline/GATKPlugin/](src/main/java/org/broadinstitute/hellbender/cmdline/GATKPlugin/). Shared argument names live in `StandardArgumentDefinitions` / `ReadFilterArgumentDefinitions`.

**Tools** ([tools/](src/main/java/org/broadinstitute/hellbender/tools/)) are grouped by domain: `walkers/` (HaplotypeCaller, Mutect2, BQSR, GenotypeGVCFs, VQSR, annotators, …), `copynumber/` (CNV), `sv/` (structural variants), `funcotator/` (functional annotation), `dragstr/`, `genomicsdb/`, `spark/`, plus `examples/` which are the canonical templates for writing a new tool of each walker type.

**Python** ([src/main/python/org/broadinstitute/hellbender/](src/main/python/org/broadinstitute/hellbender/)): `gcnvkernel` (germline CNV model), `scorevariants` (NN variant scoring), `gatktool` (Java↔Python bridge). These run inside the `gatk` Conda environment. Create/update it locally with `./gradlew localDevCondaEnv`, then `source activate gatk`.

## Conventions (from README dev guidelines)

- Follows the [Google Java Style guide](https://google.github.io/styleguide/javaguide.html). No PRs that introduce compiler warnings.
- Every tool needs ≥1 end-to-end integration test checking expected output, plus unit tests for non-trivial utilities; target >90% coverage. Every bug fix needs a regression test.
- Validate public-method arguments and throw `IllegalArgumentException` on invalid input. Prefer `final`. Don't base logic on `toString()` output. Use log4j2 (`org.apache.logging.log4j.Logger`) for logging.
- Don't push to master — open a PR; rebase and squash on merge. Keep test data files <100 KB unless they belong in `resources/large/` (lfs).

## Docs & WDL generation

```bash
./gradlew gatkDoc          # -> build/docs/gatkdoc
./gradlew gatkWDLGen       # WDL wrappers -> build/docs/wdlGen (tools annotated for WDL gen only)
```
