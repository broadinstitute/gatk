package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;


@CommandLineProgramProperties(
        summary ="Marks duplicates on dataflow",
        oneLineSummary ="Mark Duplicates",
        programGroup = DataFlowProgramGroup.class)
public final class MarkDuplicatesDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;

    @Argument(doc="output BAM file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline(final Pipeline pipeline) {
        final ReadsDataflowSource readsSource = new ReadsDataflowSource(bam, pipeline);
        final SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary):
                IntervalUtils.getAllIntervalsForReference(sequenceDictionary);

        final PCollectionView<SAMFileHeader> headerPcolView = pipeline.apply(Create.of(header)).apply(View.<SAMFileHeader>asSingleton());

        final PCollection<GATKRead> preads = readsSource.getReadPCollection(intervals);
        final PCollection<GATKRead> results = preads.apply(new MarkDuplicates(headerPcolView));

        // TODO: support writing large output files (need a sharded BAM writer)
        SmallBamWriter.writeToFile(pipeline, results, header, outputFile);
    }
}