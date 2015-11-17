package org.broadinstitute.hellbender.engine.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceTwoBitSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

/**
 * Internal interface to load a reference sequence.
 */
public interface ReferenceSource {
    ReferenceBases getReferenceBases(PipelineOptions pipelineOptions, SimpleInterval interval) throws IOException;
    SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException;

    /**
     * Returns whether this reference source can be used with Spark broadcast.
     * Currently, only {@link ReferenceTwoBitSource} is compatible with the Spark broadcast implementation.
     */
    default public boolean isCompatibleWithSparkBroadcast(){
        return this instanceof ReferenceTwoBitSource;
    }
}
