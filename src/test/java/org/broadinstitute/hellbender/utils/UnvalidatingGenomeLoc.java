package org.broadinstitute.hellbender.utils;

/**
 * GenomeLocs are very useful objects to keep track of genomic locations and perform set operations
 * with them.
 *
 * However, GenomeLocs are bound to strict validation through the GenomeLocParser and cannot
 * be created easily for small tasks that do not require the rigors of the GenomeLocParser validation
 *
 * UnvalidatingGenomeLoc is a simple utility to create GenomeLocs without going through the parser.
 *
 * WARNING: SHOULD BE USED ONLY BY EXPERT USERS WHO KNOW WHAT THEY ARE DOING!
 */
@SuppressWarnings("serial")
public final class UnvalidatingGenomeLoc extends GenomeLoc {

    public UnvalidatingGenomeLoc(String contigName, int contigIndex, int start, int stop) {
        super(contigName, contigIndex, start, stop);
    }
}
