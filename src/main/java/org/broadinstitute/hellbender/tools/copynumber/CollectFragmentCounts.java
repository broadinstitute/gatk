package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.CachedBinarySearchIntervalList;
import org.broadinstitute.hellbender.tools.exome.SampleCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableColumn;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Naive implementation of fragment-based coverage collection. The count for each interval is calculated by counting
 * how many different fragment centers intersect with this interval. The start and end positions of fragments are
 * inferred from read information. We only allow properly paired, first of pair reads - thus we do not double count
 * and we exclude reads whose fragment's position cannot be automatically inferred from its SAM record.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collect fragment counts, by counting how many fragment centers intersect with a" +
                " given interval. The fragments are inferred from SAM records of only properly paired intervals.",
        oneLineSummary = "Collect fragment counts",
        programGroup = CopyNumberProgramGroup.class
)
public class CollectFragmentCounts extends ReadWalker {

    private static Logger logger = LogManager.getLogger(CollectFragmentCounts.class);
    private SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

    public static final String COLUMN_SEPARATOR = "\t";
    public static final String LINE_SEPARATOR = "\n";

    @Argument(
            doc = "Output fragment-counts file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCountsFile = null;

    /**
     * Per interval counts, ordered according to {@link CollectFragmentCounts#cachedIntervalList}
     */
    private int[] intervalCounts;

    /**
     * List of non-overlapping searchable intervals
     */
    private CachedBinarySearchIntervalList<SimpleInterval> cachedIntervalList;

    /**
     * Name of the sample contained in the BAM file
     */
    private String sampleName;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>(super.getDefaultReadFilters());
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.FIRST_OF_PAIR); // this will make sure we don't double count
        filters.add(ReadFilterLibrary.PROPERLY_PAIRED);
        // this will only keep reads in pairs that are properly oriented and mapped on same chromosome
        // and lie within a few standard deviations from the mean of fragment size distributions
        return filters;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        if (readArguments.getReadFilesNames().size() != 1) {
            throw new UserException("This tool only accepts a single bam/sam/cram as input");
        }

        final SampleCollection sampleCollection = new SampleCollection(getHeaderForReads());
        if (sampleCollection.sampleCount() > 1) {
            throw new UserException.BadInput("We do not support bams with more than one sample.");
        }
        sampleName = sampleCollection.sampleIds().get(0);

        final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        cachedIntervalList = new CachedBinarySearchIntervalList<>(intervals);
        intervalCounts = new int[intervals.size()];

        logger.info("Starting to collect fragment counts...");
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        //getting a center of the fragment

        //TODO collect information on reads that do not have a properly paired mate
        int centerOfFragment = ReadOrientation.getCenterOfFragment(read);
        final Locatable centerFragmentLocation = new SimpleInterval(read.getContig(), centerOfFragment, centerOfFragment);

        final IndexRange intersectionRange = cachedIntervalList.findIntersectionRange(centerFragmentLocation);
        if (intersectionRange.size() > 1) {
            // should not reach here since intervals are checked for overlapping;
            // doing a check to protect against future code changes
            throw new GATKException.ShouldNeverReachHereException("At most one interval can intersect with a center of a fragment.");
        }
        intersectionRange.forEach(intervalIndex -> intervalCounts[intervalIndex]++);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Finished collecting fragment files");
        logger.info("Writing results to disk");
        writeResults(outputCountsFile, sampleName);

        logger.info("Fragment counts written to " + outputCountsFile.toString());
        return "SUCCESS";
    }

    private void writeResults(final File outputFile, final String sampleName) {
        //TODO replace this with the TableWriter
        try (PrintWriter countsFileWriter = new PrintWriter(outputFile)) {
            countsFileWriter.println(composeHeader(getCommandLine(), sampleName));
            for (int i = 0; i < intervalCounts.length; i++) {
                final SimpleInterval nextInterval = cachedIntervalList.getSortedIntervals().get(i);
                //creating the target just to get a properly formatted name that matches other current tools
                final Target targetFromInterval = new Target(nextInterval);
                String nextIntervalCountLine = String.join(COLUMN_SEPARATOR,
                        nextInterval.getContig(),
                        Integer.toString(nextInterval.getStart()),
                        Integer.toString(nextInterval.getEnd()),
                        targetFromInterval.getName(),
                        Integer.toString(intervalCounts[i]));
                countsFileWriter.println(nextIntervalCountLine);
            }
        } catch (IOException exc) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, exc);
        }
    }

    private static String composeHeader(final String commandLine, final String sampleName) {
        final List<String> tableColumnNames = new ArrayList<>();
        tableColumnNames.add(TargetTableColumn.CONTIG.toString());
        tableColumnNames.add(TargetTableColumn.START.toString());
        tableColumnNames.add(TargetTableColumn.END.toString());
        tableColumnNames.add(TargetTableColumn.NAME.toString());
        tableColumnNames.add(sampleName);
        return String.format(
                String.join(LINE_SEPARATOR,
                        "##fileFormat  = tsv",
                        "##commandLine = %s",
                        "##title       = Fragment read counts per interval",
                        String.join(COLUMN_SEPARATOR, tableColumnNames)),
                commandLine);
    }

    /**
     * Helper class to calculate fragment center of a properly paired read
     */
    @VisibleForTesting
    protected enum ReadOrientation {

        /**
         * Read was located on forward strand
         */
        FORWARD(read -> read.getUnclippedStart() + read.getFragmentLength() / 2),

        /**
         * Read was located on reverse strand
         */
        REVERSE(read -> read.getUnclippedStart() + (read.getLength() - 1)  + read.getFragmentLength() / 2);

        private final Function<GATKRead, Integer> readToFragmentCenterMapper;

        ReadOrientation(Function<GATKRead, Integer> readToCenterMapper) {
            this.readToFragmentCenterMapper = readToCenterMapper;
        }

        /**
         * Get a function that maps the read to the center of the fragment
         */
        protected Function<GATKRead, Integer> getReadToFragmentCenterMapper() {
            return readToFragmentCenterMapper;
        }

        /**
         * Get {@link ReadOrientation} instance corresponding to the orientation of the read
         */
        protected static ReadOrientation getReadOrientation(GATKRead read) {
            return read.getFragmentLength() > 0 ? FORWARD : REVERSE;
        }

        /**
         * Compute center of the fragment that read corresponds to
         */
        protected static int getCenterOfFragment(GATKRead read) {
            return getReadOrientation(read).getReadToFragmentCenterMapper().apply(read);
        }
    }
}
