package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
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
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKGroupedReadNullCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.File;
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

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Argument(doc = "File to write duplication metrics to.", optional=true,
        shortName = "M", fullName = "METRICS_FILE")
    protected File metricsFile;

    // use this for large files now, since we would crash otherwise.
    @Argument(doc = "Skips saving the output.", optional=true,
        fullName = "dontSave")
    protected boolean dontSave;

    // the default values are fine, but this allows for fine-tuning if necessary.
    @Argument(doc = "How many bases do we query from disk at once", optional=true,
        fullName = "readShardSize")
    protected int readShardSize = 1_000_000;

    // the default values are fine, but this allows for fine-tuning if necessary.
    @Argument(doc = "How many bases do we process at a once", optional=true,
        fullName = "processShardSize")
    protected int processShardSize = 5_000;


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
        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
            new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
        final PCollectionView<OpticalDuplicateFinder> finderPcolView = pipeline.apply(Create.of(finder)).setName("OpticalDuplicateFinder").apply(View.<OpticalDuplicateFinder>asSingleton());

        // we read this many at a time
        int basesPerShard = readShardSize;
        // and then output groups this large
        int outputBasesPerWorkUnit = processShardSize;

        final PCollection<KV<String, Iterable<GATKRead>>> readsByShard = readsDataflowSource.getGroupedReadPCollection(intervals, basesPerShard, outputBasesPerWorkUnit, pipeline);
        final PCollection<GATKRead> results = readsByShard.apply(new MarkDuplicatesFromShardsDataflowTransform(headerSingleton, finderPcolView));

            // TODO: support writing large output files (need a sharded BAM writer)
        if (!dontSave) {
            SmallBamWriter.writeToFile(pipeline, results, readsHeader, outputFile);
        }

        if (metricsFile != null) {
            final PCollection<KV<String,DuplicationMetrics>> metrics = results.apply(new MarkDuplicatesUtils.GenerateMetricsTransform(headerSingleton));
            MarkDuplicatesUtils.writeMetricsToFile(pipeline, metrics, readsHeader, metricsFile);
        }

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