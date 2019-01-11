package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

/**
 * Internal interface to load a reference sequence.
 */
public interface ReferenceSparkSource {
    ReferenceBases getReferenceBases(SimpleInterval interval) throws IOException;
    SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException;

    /**
     * Returns whether this reference source can be used with Spark broadcast.
     * Currently, only {@link ReferenceTwoBitSparkSource} is compatible with the Spark broadcast implementation.
     */
    default public boolean isCompatibleWithSparkBroadcast(){
        return this instanceof ReferenceTwoBitSparkSource;
    }
}
