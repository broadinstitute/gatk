package org.broadinstitute.hellbender.tools.spark.utils;

/**
 * Interface for classes that provide ways to add long primitives and test for set membership
 */
public interface QueryableLongSet {

    boolean add(final long val); //TODO needed for experimental kmer filtering approach
    boolean contains(final long val);
    boolean containsAll(final long[] vals);

}
