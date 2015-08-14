package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.reference.HadoopReferenceSequenceFileFactory;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.apache.hadoop.fs.Path;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;
import java.util.Map;

public class ReferenceHadoopSource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    private final String referencePath;

    public ReferenceHadoopSource(final String referencePath) {
        this.referencePath = referencePath;
    }

    public ReferenceBases getReferenceBases(SimpleInterval interval, PipelineOptions pipelineOptions) throws IOException {
        ReferenceSequenceFile referenceSequenceFile = HadoopReferenceSequenceFileFactory.getReferenceSequenceFile(new Path(referencePath));
        ReferenceSequence sequence = referenceSequenceFile.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
        return new ReferenceBases(sequence.getBases(), interval);
    }

}
