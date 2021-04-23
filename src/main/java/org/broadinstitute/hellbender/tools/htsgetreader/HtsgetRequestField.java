package org.broadinstitute.hellbender.tools.htsgetreader;

/**
 * Fields which can be used to filter a htsget request as defined by the spec
 */
public enum HtsgetRequestField {
    QNAME,
    FLAG,
    RNAME,
    POS,
    MAPQ,
    CIGAR,
    RNEXT,
    PNEXT,
    TLEN,
    SEQ,
    QUAL
}
