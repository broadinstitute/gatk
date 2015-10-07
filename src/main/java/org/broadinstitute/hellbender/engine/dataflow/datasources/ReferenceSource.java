package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;

/**
 * Internal interface to load a reference sequence.
 */
public interface ReferenceSource {
    ReferenceBases getReferenceBases(PipelineOptions pipelineOptions, SimpleInterval interval) throws IOException;
    SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException;
}
