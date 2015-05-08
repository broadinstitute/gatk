package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.dataflow.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Handles lifting reads from a bams into a PCollection and provides hooks to apply {@link ReadFilter}, {@link ReadConverter}, and
 * a {@link PTransformSAM}.
 *
 *Subclasses must override {@link #getTool()} and optionally override {@link #getReadFilters()} and {@link #getReadTransformers()}
 */
public abstract class DataflowReadsPipeline extends DataflowCommandLineProgram {

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
    abstract protected PTransformSAM<?> getTool();

    /**
     * @return {@link ReadFilter}s to apply before running the tool
     */
    protected ImmutableList<ReadFilter> getReadFilters(){
        return ImmutableList.of();
    }

    /**
     * @return {@link ReadTransformer}s to apply before running the tool
     */
    protected ImmutableList<ReadTransformer> getReadTransformers(){
        return ImmutableList.of();
    }

    private static SerializableFunction<Read,Boolean> wrapFilter(ReadFilter filter, String headerString){
        return new SAMSerializableFunction<>(headerString, filter);
    }

    @Override
    final protected void setupPipeline(Pipeline pipeline) {
        ReadsSource readsSource = new ReadsSource(bam, pipeline);
        String headerString = readsSource.getHeaderString();
        SAMSequenceDictionary sequenceDictionary = ReadUtils.samHeaderFromString(headerString).getSequenceDictionary();
        List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary):
                getAllIntervalsForReference(sequenceDictionary);

        PCollection<Read> preads = readsSource.getReadPCollection(intervals);

        PCollection<?> presult = applyTransformsToPipeline(headerString, preads);

        PCollection<String> pstrings = presult.apply(DataflowUtils.convertToString());
        pstrings.apply(TextIO.Write.to(outputFile));
    }

    private List<SimpleInterval> getAllIntervalsForReference(SAMSequenceDictionary sequenceDictionary) {
        return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
                .stream()
                .map(SimpleInterval::new)
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    protected PCollection<?> applyTransformsToPipeline(String headerString, PCollection<Read> preads) {
        for (ReadFilter filter : getReadFilters()) {
            preads = preads.apply(Filter.by(wrapFilter(filter, headerString)));
        }

        for (ReadTransformer transformer : getReadTransformers()){
            preads = preads.apply(wrapTransformer(transformer, headerString));
        }

        PTransformSAM<?> f = getTool();
        f.setHeaderString(headerString);

        return preads.apply(f);
    }

    private static PTransform<? super PCollection<Read>,PCollection<Read>> wrapTransformer(ReadTransformer transformer, String headerString){
        return ParDo.of(new DataFlowSAMFn<Read>(headerString){

                    @Override
                    protected void apply(SAMRecord read) {
                        output(ReadConverter.makeRead(transformer.apply(read)));
                    }
                }
        );
    }
}
