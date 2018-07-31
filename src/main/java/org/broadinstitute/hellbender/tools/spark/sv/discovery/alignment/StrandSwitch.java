package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

/**
 * For symbolizing the change of strand from one alignment to the next
 * of an assembly contig.
 */
public enum StrandSwitch {
    NO_SWITCH, FORWARD_TO_REVERSE, REVERSE_TO_FORWARD;
}
