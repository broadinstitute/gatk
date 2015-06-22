package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowReadFn;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Handles lifting reads from a bams into a PCollection and provides hooks to apply {@link ReadFilter} and
 * a {@link PTransformSAM}.
 *
 *Subclasses must override {@link #getTool()} and optionally override {@link #getReadFilters(SAMFileHeader)} and {@link #getReadTransformers()}
 */
public abstract class DataflowReadsPipeline extends DataflowCommandLineProgram {

    private static final long serialVersionUID = 1l;

    @Argument(doc="a prefix for the dataflow output files", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    /**
     * Returns a transform which performs the work of the given tool.
     * Subclasses must override this.
     * @return the Transform that should be performed in this pipeline.
     */
    protected abstract PTransformSAM<?> getTool();

    /**
     * @return {@link ReadFilter}s to apply before running the tool
     * @param header
     */
    protected ImmutableList<ReadFilter> getReadFilters( SAMFileHeader header ){
        return ImmutableList.of();
    }

    /**
     * @return {@link ReadTransformer}s to apply before running the tool
     */
    protected ImmutableList<ReadTransformer> getReadTransformers(){
        return ImmutableList.of();
    }

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        final ReadsDataflowSource readsDataflowSource = new ReadsDataflowSource(bam, pipeline);
        final SAMFileHeader header = readsDataflowSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary):
                getAllIntervalsForReference(sequenceDictionary);

        final PCollection<GATKRead> preads = readsDataflowSource.getReadPCollection(intervals, ValidationStringency.SILENT);

        final PCollection<?> presult = applyTransformsToPipeline(header, preads);

        final PCollection<String> pstrings = presult.apply(DataflowUtils.convertToString());
        pstrings.apply(TextIO.Write.to(outputFile));
    }

    private List<SimpleInterval> getAllIntervalsForReference(SAMSequenceDictionary sequenceDictionary) {
        return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
                .stream()
                .map(SimpleInterval::new)
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    protected PCollection<?> applyTransformsToPipeline(final SAMFileHeader header, final PCollection<GATKRead> preadsIn) {

        PCollection<GATKRead> preads = preadsIn;
        for (final ReadFilter filter : getReadFilters(header)) {
            preads = preads.apply(Filter.by(wrapReadFilter(filter)));
        }

        for (final ReadTransformer transformer : getReadTransformers()) {
            preads = preads.apply(wrapTransformer(transformer, header));
        }

        final PTransformSAM<?> tool = getTool();
        tool.setHeader(header);

        return preads.apply(tool);
    }

    private static SerializableFunction<GATKRead, Boolean> wrapReadFilter( final ReadFilter filter ) {
        return new SerializableFunction<GATKRead, Boolean>() {
            private static final long serialVersionUID = 1l;

            @Override
            public Boolean apply( GATKRead input ) {
                return filter.apply(input);  // Necessary because we can't use wildcard types for our SerializableFunction
            }
        };
    }

    private static PTransform<? super PCollection<GATKRead>,PCollection<GATKRead>> wrapTransformer(final ReadTransformer transformer, final SAMFileHeader header){
        return ParDo.of(new DataFlowReadFn<GATKRead>(header){
                    private static final long serialVersionUID = 1l;

                    @Override
                    protected void apply(final GATKRead read) {
                        output(transformer.apply(read));
                    }
                }
        );
    }
}
