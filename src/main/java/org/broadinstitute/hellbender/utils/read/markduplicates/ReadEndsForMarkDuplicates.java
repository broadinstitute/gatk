package org.broadinstitute.hellbender.utils.read.markduplicates;

/**
 * Little struct-like class to hold read pair (and fragment) end data for MarkDuplicatesWithMateCigar
 *
 * @author Nils Homer
 */
public final class ReadEndsForMarkDuplicates extends ReadEnds {
    /*
    What do we need to store you ask?  Well, we need to store:
       - byte: orientation
       - short: libraryId, readGroup, tile, x, y, score
       - int: read1ReferenceIndex, read1Coordinate, read2ReferenceIndex, read2Coordinate
       - long: read1IndexInFile, read2IndexInFile
     */
    public static final int SIZE_OF = (1 * 1) + (5 * 2) + (4 * 4) + (8 * 2) + 1
            + 8 + // last 8 == reference overhead
            13; // This is determined experimentally with JProfiler

    public short score = 0;
    public long read1IndexInFile = -1;
    public long read2IndexInFile = -1;
}