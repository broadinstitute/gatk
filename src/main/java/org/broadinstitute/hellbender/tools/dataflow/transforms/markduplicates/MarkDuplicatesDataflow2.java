package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.values.KV;
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
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
        summary ="Marks duplicates on dataflow",
        oneLineSummary ="Mark Duplicates",
        programGroup = DataFlowProgramGroup.class)
public final class MarkDuplicatesDataflow2 extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(doc="output BAM file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected String jobName() {
        // surprisingly, Dataflow requires that the name be unique.
        int randomThreeDigits = new Random().nextInt(900)+100;
        return "MarkDuplicatesDataflow2-"+System.getProperty("user.name")+"-"+randomThreeDigits;
    }

    @Override
    protected void setupPipeline(final Pipeline pipeline) {
        // Load the reads.
        final ReadsDataflowSource readsDataflowSource = new ReadsDataflowSource(bam, pipeline);
        final SAMFileHeader readsHeader = readsDataflowSource.getHeader();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
            : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        final PCollectionView<SAMFileHeader> headerSingleton = ReadsDataflowSource.getHeaderView(pipeline, readsHeader);

        final PCollection<KV<String, Iterable<GATKRead>>> readsByShard = readsDataflowSource.getGroupedReadPCollection(intervals, 1_000_000, 200_000, pipeline);
        final PCollection<GATKRead> results = readsByShard.apply(new MarkDuplicatesFromShardsDataflowTransform(headerSingleton));

        // TODO: support writing large output files (need a sharded BAM writer)
        //SmallBamWriter.writeToFile(pipeline, results, readsHeader, outputFile);
    }


    /**
     * Returns the list of intervals from the given sequence dictionary.
     */
    private static List<SimpleInterval> getAllIntervalsForReference(final SAMSequenceDictionary sequenceDictionary) {
    return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
            .stream()
            .map(SimpleInterval::new)
            .collect(Collectors.toList());
    }


}