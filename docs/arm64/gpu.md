# GPU acceleration on Apple Silicon — what is and isn't realistic

People often ask for "GPU like DRAGEN" on a Mac. Two clarifications up front:

- **DRAGEN is FPGA, not GPU.** The GPU-accelerated GATK is **NVIDIA Clara Parabricks**, which
  is **CUDA**.
- **Apple Silicon has no CUDA.** The GPU is programmed with **Metal**. So there is no drop-in
  GPU path for GATK on a Mac.

## Realistic now: PyTorch MPS for the deep-learning tools

GATK's neural-network tools run in PyTorch, and PyTorch ships an **MPS (Metal Performance
Shaders)** backend that uses the Apple GPU. Tools that benefit:

- `NVScoreVariants`
- `TrainVariantAnnotationsModel` / `ScoreVariantAnnotations`

Plan: in the `scorevariants` Python package, auto-select the device in priority order
**CUDA → MPS → CPU** (`torch.backends.mps.is_available()`), and stop pinning the MKL/CPU-only
PyTorch build in the Conda env. This is a bounded, real GPU win on Apple Silicon.

Caveats: some torch ops still fall back to CPU on MPS; set
`PYTORCH_ENABLE_MPS_FALLBACK=1` so unsupported ops don't hard-fail. Validate scoring output
parity against the CPU path.

## Not realistic here: GPU for the core variant-calling kernels

DRAGEN/Parabricks-style acceleration of **PairHMM, Smith-Waterman, and bwa** on a Mac GPU
would require rewriting those kernels as **Metal compute shaders** from scratch — there is no
existing Metal implementation to reuse, and the numeric validation burden is large. This is a
research-grade effort and is **out of scope** for the native-arm64 port. It is recorded here
as a possible future direction; the CPU NEON path (via SIMDe / hand-NEON in the rebuilt GKL)
is the supported acceleration for those kernels on Apple Silicon.

## gcnvkernel (PyMC / PyTensor)

No Metal backend exists; PyTensor would need a JAX-Metal path which is experimental and
fragile. These tools stay on CPU (NEON via the native BLAS / Accelerate).
