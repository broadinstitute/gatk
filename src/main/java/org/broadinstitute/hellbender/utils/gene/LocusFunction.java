package org.broadinstitute.hellbender.utils.gene;

/**
 * Describes the behavior of a locus relative to a gene.  Note that these are enumerated in a specific order:
 * (INTERGENIC, INTRONIC, UTR, CODING, RIBOSOMAL) so that e.g. if a locus is CODING in one gene and UTR in another,
 * we count it as CODING.
 */
public enum LocusFunction {
    INTERGENIC, INTRONIC, UTR, CODING, RIBOSOMAL
}
