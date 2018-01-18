package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.engine.GATKTool;

/**
 * Tool to extract an interval from the given reference and create a new FASTA file that contains only data
 * from that interval.
 *
 * If only a contig is specified, a FASTA file containing the whole contig will be created.
 * If a contig and a start position are specified, a FASTA file containing all bases starting from the given start position to the end of the contig will be created.
 * If a contig and an end position are specified, a FASTA file containing all bases starting from the beginning of the contig to the end position will be created.
 * If a contig, start, and end position are specified, a FASTA file containing all bases from the start to the end position on the given contig will be created.
 *
 * Inputs:
 *     Contig (required)
 *     Start Position (1-based, inclusive) (optional)
 *     End Position (1-based, inclusive) (optional)
 *
 * Outputs:
 *     FASTA file
 *     FASTA index
 *     FASTA sequence dictionary
 */
public class ExtractSequence extends GATKTool {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public void traverse() {

    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
