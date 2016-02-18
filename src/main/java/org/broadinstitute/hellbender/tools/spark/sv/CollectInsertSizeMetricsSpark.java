package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;

@CommandLineProgramProperties(
        summary        = "Program to collect insert size distribution information in SAM/BAM file(s)",
        oneLineSummary = "Collect Insert Size Distribution on Spark",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class CollectInsertSizeMetricsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "A local path to file to write insert size metrics to, with extension.",
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              optional = true)
    public String OUTPUT = null;

    @Argument(doc = "A local path to PDF file where histogram plot will be saved in.",
              shortName = "HIST",
              fullName = "HistogramPlotPDF",
              optional = true)
    public String HISTOGRAM_PLOT_FILE = null;

    @Argument(doc = "Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS_TOL*MEDIAN_ABSOLUTE_DEVIATION. " +
                    "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
                    "artifacts to make the mean and sd grossly misleading regarding the real distribution.",
              shortName = "TOL",
              fullName = "HistogramPlotDeviationsTolerance",
              optional = true)
    public double DEVIATIONS_TOL = 10.0;

    // read filtering criteria
    @Argument(doc = "If set to true, use pairs of reads that are not properly oriented.",
              shortName = "nPP",
              fullName = "useNonProperlyPairedReads",
              optional = true)
    public boolean useNonProperlyPairedReads = false;

    @Argument(doc = "If set to true, include unmapped reads as well.",
              shortName = "U",
              fullName = "useUnmappedReads",
              optional = true)
    public boolean useUnmappedReads = false;

    @Argument(doc = "If set to true, include reads whose mate is unmapped as well.",
              shortName = "Um",
              fullName = "useMateUnmappedReads",
              optional = true)
    public boolean useMateUnmappedReads = false;

    @Argument(doc = "If set to true, include duplicated reads as well.",
              shortName = "Dup",
              fullName = "useDupReads",
              optional = true)
    public boolean useDupReads = false;

    @Argument(doc = "If set to true, include secondary alignments.",
              shortName = "S",
              fullName = "useSecondaryAlignments",
              optional = true)
    public boolean useSecondaryAlignments = false;

    @Argument(doc = "If set to true, include supplementary alignments.",
              shortName = "SS",
              fullName = "useSupplementaryAlignments",
              optional = true)
    public boolean useSupplementaryAlignments = false;

    @Argument(doc = "If set non-zero value, only include reads passing certain mapping quality threshold. " +
                    "If set to zero, reads with zero mapping quality will be included in calculating metrics.",
              shortName = "Q",
              fullName = "MQPassingThreshold",
              optional = true)
    public int MQPassingThreshold = 0;

//    TODO
//    @Argument(doc = "If set to non-zero value, only include reads with mate passing certain mapping quality threshold.",
//              shortName = "Qm",
//              fullName = "mateMQPassingThreshold",
//              optional = true)
//    public int mateMQPassingThreshold = 0;

//    @Argument(doc="The level(s) at which to accumulate metrics. Possible values are {SAMPLE, LIBRARY, READ GROUP}. " +
//                  "Default: All reads.",
//              shortName="LEVEL",
//              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
//              optional = false)
//    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    /**
     * Which end of a read pair to use for collecting insert size metrics.
     * For truncated genomic regions of putative SV breakpoints, not all reads have both ends land in the region,
     *     so third case is possible: will use either end when only one end is available in the region specified,
     *     and only first end if both are available
     */
    public enum UseEnd {
        FIRST(1), SECOND(2);//, EITHER(0);
        private final int value;
        UseEnd(int value){
            this.value = value;
        }
        public int getValue(){
            return value;
        }
    }

    @Argument(doc = "Which end of pairs to use for collecting information. " +
                    "Possible values:{FIRST, SECOND, EITHER}." +
                    "Option EITHER picks up information from 1st end when both ends are available, " +
                    "and pick either end when only one end is available. " +
                    "(Remember, in SV analysis, a truncated region may be investigated, so this is possible.)",
              shortName = "E",
              fullName = "whichEndOfPairToUse",
              optional = true)
    public UseEnd useEnd = UseEnd.FIRST;

    // path to Picard R script for producing histograms in PDF files.
    private static final String R_SCRIPT = "insertSizeHistogram.R";

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> filteredReads = getReads();

        // Class where real metric-collection work is delegated to.
        final InsertSizeMetricsCollectorSpark collector = new InsertSizeMetricsCollectorSpark(filteredReads,
                                                                                              // METRIC_ACCUMULATION_LEVEL,
                                                                                              getHeaderForReads().getReadGroups(),
                                                                                              DEVIATIONS_TOL);

        final SAMFileHeader readsHeader = getHeaderForReads();
        final String inputFileName = getReadSourceName();

        if(OUTPUT != null) {
            writeMetricsFile(collector, readsHeader, inputFileName);
        }
        if(HISTOGRAM_PLOT_FILE != null){
            writeHistogramPDF(readsHeader, inputFileName);
        }
    }

    /**
     *  Implicitly called in getReads(): behavior from base (GATKTools). Return type serializable.
     *  Not calling base(GATKTools) because custom filters is not a superset of filters defined there.
     */
    @Override
    public ReadFilter makeReadFilter() {

        final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(getHeaderForReads());

        final PrimaryFilter pfilter = new PrimaryFilter();

        final SVCustomReadFilter sfilter = new SVCustomReadFilter(useEnd,
                                                                  !useNonProperlyPairedReads,
                                                                  !useUnmappedReads,
                                                                  !useMateUnmappedReads,
                                                                  !useDupReads,
                                                                  !useSecondaryAlignments,
                                                                  !useSupplementaryAlignments,
                                                                  MQPassingThreshold);

        return alignmentAgreesWithHeader.and(pfilter).and(sfilter);
    }

    /**
     * Similar to WellformedReaFilter except doesn't check for header but checks if mapping quality is available.
     */
    private static final class PrimaryFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        private final ReadFilter primary = ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK
                                           .and(ReadFilterLibrary.VALID_ALIGNMENT_START)
                                           .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                                           .and(ReadFilterLibrary.HAS_READ_GROUP)
                                           .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                                           .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                                           .and(ReadFilterLibrary.SEQ_IS_STORED)
                                           .and(ReadFilterLibrary.CIGAR_IS_SUPPORTED)
                                           .and(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);

        @Override
        public boolean test(final GATKRead read){
            return primary.test(read);
        }
    }

    /**
     * Customized serializable reads filter, based on cmd line arguments provided
     */
    private static final class SVCustomReadFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        private final ReadFilter collapsedFilter;

        public SVCustomReadFilter(final UseEnd whichEnd,
                                  final boolean filterNonProperlyPairdReads,
                                  final boolean filterUnmappedReads,
                                  final boolean filterMateUnmappedReads,
                                  final boolean filterDupReads,
                                  final boolean filterSecondaryAlignments,
                                  final boolean filterSupplementaryAlignments,
                                  final int MQThreshold){

                final UseEnd endVal = whichEnd;

                ReadFilter tempFilter = read -> 0!=read.getFragmentLength();
                           tempFilter = tempFilter.and(read -> read.isPaired());
                           tempFilter = tempFilter.and(read -> endVal == (read.isFirstOfPair() ? UseEnd.FIRST : UseEnd.SECOND));

                if(filterNonProperlyPairdReads)   { tempFilter = tempFilter.and(read -> read.isProperlyPaired());}
                if(filterUnmappedReads)           { tempFilter = tempFilter.and(read -> !read.isUnmapped());}
                if(filterMateUnmappedReads)       { tempFilter = tempFilter.and(read -> !read.mateIsUnmapped());}
                if(filterDupReads)                { tempFilter = tempFilter.and(read -> !read.isDuplicate());}
                if(filterSecondaryAlignments)     { tempFilter = tempFilter.and(read -> !read.isSecondaryAlignment());}
                if(filterSupplementaryAlignments) { tempFilter = tempFilter.and(read -> !read.isSupplementaryAlignment());}

                if(0!=MQThreshold){ tempFilter = tempFilter.and(read -> read.getMappingQuality() >= MQThreshold);}

                collapsedFilter = tempFilter;

                // TODO: efficiently checking mate mapping quality (may require index information)
                // TODO: pick only isReverseStrand/mateIsReverseStrand (strand bias)
                // TODO: case 0 of "whichEnd", most complicated since behavior depends on if mate is available in interesting region
                // TODO: filter or not based on length==0
                // TODO: filter based on MATE_ON_SAME_CONTIG MATE_DIFFERENT_STRAND GOOD_CIGAR NON_ZERO_REFERENCE_LENGTH_ALIGNMENT
        }

        @Override
        public boolean test(final GATKRead read){
            return collapsedFilter.test(read);
        }
    }

    @VisibleForTesting
    void writeMetricsFile(final InsertSizeMetricsCollectorSpark collector, final SAMFileHeader readsHeader, final String inputFileName){

        final MetricsFile<InsertSizeMetrics, Integer> metrics = getMetricsFile();

        collector.produceMetricsFile(metrics);

        MetricsUtils.saveMetrics(metrics, OUTPUT, getAuthHolder());

        if (metrics.getAllHistograms().isEmpty()) {
            logger.warn("No valid reads found in input file.");
        }
    }

    @VisibleForTesting
    void writeHistogramPDF(final SAMFileHeader readsHeader, final String inputFileName){

        if(0.0 == DEVIATIONS_TOL){
            logger.warn("MAD tolerance for histogram set to 0, no plot to generate.");
        }

        final File histogramPlotPDF = new File(HISTOGRAM_PLOT_FILE);
        IOUtil.assertFileIsWritable(histogramPlotPDF);

            /* Temporarily commented out because Picard R script for this purpose doesn't produce subtitle based on library
            // If we're working with a single library, assign that library's name as a suffix to the plot title
            final List<SAMReadGroupRecord> readGroups = readsHeader.getReadGroups();

            //A subtitle for the plot, usually corresponding to a library.
            String plotSubtitle = "";
            if (readGroups.size() == 1) {
                plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
            }*/

        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, CollectInsertSizeMetricsSpark.class));
        executor.addArgs(OUTPUT,                                // text-based metrics file
                         histogramPlotPDF.getAbsolutePath(),    // PDF graphics file
                         inputFileName);                        // input bam file
        executor.exec();
    }
}