package org.broadinstitute.hellbender.utils;

/**
 * In broad terms, each sequencing platform can be classified by whether it flows nucleotides in some order
 * such that homopolymers get sequenced in a single event (ie 454 or Ion) or it reads each position in the
 * sequence one at a time, regardless of base composition (Illumina or Solid).  This information is primarily
 * useful in the BQSR process
 */
public enum SequencerFlowClass {
    DISCRETE,
    FLOW,
    OTHER //Catch-all for unknown platforms, as well as relics that GATK doesn't handle well (Capillary, Helicos)
}
