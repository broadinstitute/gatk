package org.broadinstitute.hellbender.utils.samples;

/**
 *
 */
public enum PedigreeValidationType {
    /**
     * Require if a pedigree file is provided at all samples in the VCF or BAM files have a corresponding
     * entry in the pedigree file(s).
     */
    STRICT,

    /**
     * Do not enforce any overlap between the VCF/BAM samples and the pedigree data
     * */
    SILENT
}
