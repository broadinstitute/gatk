package org.broadinstitute.hellbender.engine.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceTwoBitSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

/**
 * Internal interface to load a reference sequence.
 */
public interface ReferenceSource {

    SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException;

    default ReferenceBases getReferenceBases(final SimpleInterval interval) throws IOException {
        return getReferenceBases(interval);
    }

    default SAMSequenceDictionary getReferenceSequenceDictionary() throws IOException {
        return getReferenceSequenceDictionary(null);
    }

    /**
     * Returns whether this reference source can be used with Spark broadcast.
     * Currently, only {@link ReferenceTwoBitSource} is compatible with the Spark broadcast implementation.
     */
    default boolean isCompatibleWithSparkBroadcast(){
        return this instanceof ReferenceTwoBitSource;
    }
}
