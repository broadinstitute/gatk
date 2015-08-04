package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.*;
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
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;

import java.io.File;
import java.util.List;


@CommandLineProgramProperties(
        summary ="Marks duplicates on dataflow",
        oneLineSummary ="Mark Duplicates",
        programGroup = DataFlowProgramGroup.class)
public final class MarkDuplicatesDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;

    // Used to set an attribute on the GATKRead marking this read as an optical duplicate.
    private static final String OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME = "OD";

    @Argument(doc="output BAM file", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @Argument(doc = "File to write duplication metrics to.", optional=true,
              shortName = "M", fullName = "METRICS_FILE")
    protected File metricsFile;

    @Argument(doc = "Regular expression that can be used to parse read names in the incoming SAM file. Read names are " +
             "parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used " +
             "to estimate the rate of optical duplication in order to give a more accurate estimated library size. " +
             "Set this option to null to disable optical duplicate detection. " +
             "The regular expression should contain three capture groups for the three variables, in order. " +
             "It must match the entire read name. " +
             "Note that if the default regex is specified, a regex match is not actually done, but instead the read name " +
             " is split on colon character. " +
             "For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. " +
             "For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.",
             optional = true)
    public String READ_NAME_REGEX = OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX;

    @Argument(doc = "The maximum offset between two duplicate clusters in order to consider them optical duplicates. This " +
             "should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the " +
             "Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.")
    public int OPTICAL_DUPLICATE_PIXEL_DISTANCE = OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE;

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

        final OpticalDuplicateFinder finder = READ_NAME_REGEX != null ?
            new OpticalDuplicateFinder(READ_NAME_REGEX, OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
        final PCollectionView<OpticalDuplicateFinder> finderPcolView = pipeline.apply(Create.of(finder)).apply(View.<OpticalDuplicateFinder>asSingleton());

        final PCollection<GATKRead> results = preads.apply(new MarkDuplicates(headerPcolView, finderPcolView));

        // TODO: support writing large output files (need a sharded BAM writer)
        SmallBamWriter.writeToFile(pipeline, results, header, outputFile);

        if (metricsFile != null) {
            final PCollection<KV<String,DuplicationMetrics>> metrics = results.apply(new MarkDuplicatesUtils.GenerateMetricsTransform(headerPcolView));
            MarkDuplicatesUtils.writeMetricsToFile(pipeline, metrics, header, metricsFile);
        }
    }
}
