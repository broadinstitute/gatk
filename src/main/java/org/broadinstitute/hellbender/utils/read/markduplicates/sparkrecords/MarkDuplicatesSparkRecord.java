package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

import java.util.Map;

/**
 * A common interface for the data types that represent reads for mark duplicates spark. This is done to satisfy limitations
 * in spark as it is often faster to perform a map on an RDD than to split it off into multiple streams.
 *
 * The implementing classes are specific to the needs of various code paths in MarkDuplicatesSpark.
 */
public abstract class MarkDuplicatesSparkRecord {
    protected final int partitionIndex;
    protected final String name;
    MarkDuplicatesSparkRecord(int partitionIndex, String name) {
        this.name = name;
        this.partitionIndex = partitionIndex;
    }

    // Required abstract methods
    public abstract Type getType();
    // NOTE: these keys are typically stored as a transient field to prevent serialization for performances purposes, this is only guaranteed to give the correct result right after construction
    public abstract ReadsKey key();


    // A fragment containing only one read without a mapped mate
    public static Fragment newFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap) {
        return new Fragment (first, header, partitionIndex, scoringStrategy, headerLibraryMap);
    }

    // An optimization for reducing the serialized data passed around when indicating that there was a mapped read at a location
    public static EmptyFragment newEmptyFragment(GATKRead read, SAMFileHeader header, Map<String, Byte> headerLibraryMap) {
        return new EmptyFragment(read, header, headerLibraryMap);
    }

    // An object representing a pair of primary and secondary reads with a particular span for duplicate marking
    public static Pair newPair(GATKRead first, GATKRead second, SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap) {
        return new Pair(first, second, header, partitionIndex, scoringStrategy, headerLibraryMap);
    }

    // An object representing a read or group of reads that we want to pass through the tool without being duplicate marked
    public static Passthrough getPassthrough(GATKRead read, int partitionIndex) {
        return new Passthrough(read, partitionIndex);
    }


    public int getPartitionIndex(){
      return partitionIndex;
    }

    public String getName() {
      return name;
    }

    public enum Type {
        FRAGMENT, PAIR, PASSTHROUGH, EMPTY_FRAGMENT
    }
}
