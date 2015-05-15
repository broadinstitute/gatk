package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

public class ReadsPreprocessingPipeline extends DataflowCommandLineProgram {

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @Argument(doc = "", shortName = "IRKnownVariants", fullName = "indelRealignmentKnownVariants", optional = true)
    protected List<String> indelRealignmentKnownVariants;

    @Argument(doc = "", shortName = "BQSRKnownVariants", fullName = "baseRecalibrationKnownVariants", optional = true)
    protected List<String> baseRecalibrationKnownVariants;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline( Pipeline pipeline ) {
        final ReadsSource readsSource = new ReadsSource(bam, pipeline);
        final String headerString = readsSource.getHeaderString();
        final SAMFileHeader readsHeader = ReadUtils.samHeaderFromString(headerString);
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary()):
                                                                                                 getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        // TODO: Using a header string for now as a stub. Need a custom coder for SAMFileHeader,
        // TODO: and SAMFileHeader needs to be Serializable
        final PCollectionView<String> headerSingleton = pipeline.apply(Create.of(headerString)).setCoder(StringUtf8Coder.of()).apply(View.<String>asSingleton());
        final PCollection<Read> initialReads = readsSource.getReadPCollection(intervals);

        final PCollection<Read> markedReads = initialReads.apply(new MarkDuplicatesStub().withHeader(headerSingleton));
    }

    private List<SimpleInterval> getAllIntervalsForReference(SAMSequenceDictionary sequenceDictionary) {
        return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
                .stream()
                .map(SimpleInterval::new)
                .collect(Collectors.toList());
    }

    // NOTE: need a way to ensure that certain tools are guaranteed to have a header -- one option is
    // an interface with a factory method
    private static class MarkDuplicatesStub extends PTransform<PCollection<Read>, PCollection<Read>> {

        private PCollectionView<String> header;

        @Override
        public PCollection<Read> apply( PCollection<Read> input ) {
            return input.apply(ParDo.named("MarkDuplicates").
                    of(new DoFn<Read, Read>() {
                        @Override
                        public void processElement( ProcessContext c ) throws Exception {
                            c.output(c.element());
                        }
                    }).withSideInputs(header));
        }

        public PTransform<PCollection<Read>, PCollection<Read>> withHeader( final PCollectionView<String> header ) {
            this.header = header;
            return this;
        }
    }
}
