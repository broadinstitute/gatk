==============================
GATK gCNV computational kernel
==============================

`gcnvkernel` is a python module that provides the computational backend for GATK
germline copy number variant (gCNV) tools and workflows.

This module implements inference schemes for read-depth profile denoising, germline
integer copy number variation discovery, germline contig ploidy determination, associated
I/O methods, and helper CLI scripts. `gcnvkernel` additionally provides general-purpose
inference schemas built on the top of `PyMC3` and `theano`.

The module is organized as follows::

    `gcnvkernel.inference`: Classes related to the general task of variational inference,
        such as convergence tracking, deterministic annealing protocol, and structured
        stochastic optimizers.

    `gcnvkernel.io`: I/O routines for saving and loading various data structures, such as
        genomic intervals, read count data, global and sample-specific posteriors, and
        sample metadata.

    `gcnvkernel.models`: `PyMC3` model declarations, `theano` symbolic operations (e.g.
        forward-backward algorithm for HMMs), and custom probability distributions.

    `gcnvkernel.preprocess`: Routines for filtering interval lists.

    `gcnvkernel.postprocess`: Classes for post-processing gcnvkernel output (e.g.
        Viterbi segmentation and segment quality calculation).

    `gcnvkernel.structs`: Classes that represents genomic intervals, interval annotations,
        and sample metadata.

    `gcnvkernel.tasks`: Inference task executors for different models and their use cases.
        This module additionally provides general-purpose templates for setting up arbitrary
        ADVI inference tasks for models with mixed discrete and continuous latent variables.

    `gcnvkernel.utils`: Contains utility classes for CLI scripts, an implementation of
        recursive least squares (RLS), and a number of useful mathematical functions.

Even though `gcnvkernel` can be used as a standalone tool, its intended use case is to be
privately called by GATK CLI tools. As such, the outputs of `gcnvkernel` are considered
_intermediate_ and to be further processed and properly formatted by GATK. Advanced users
may choose to use all or parts of `gcnvkernel` methods directly, either using the python CLI
scripts provided in `src/main/resources/org/broadinstitute/hellbender/tools/copynumber`, 
or in interactive python environments for exploratory analyses.

