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
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.R.RScriptExecutorException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.Set;
import java.util.EnumSet;

// TODO: filter reads based on only isReverseStrand/mateIsReverseStrand (strand bias)
// TODO: filter reads based on {MATE_ON_SAME_CONTIG, MATE_DIFFERENT_STRAND, GOOD_CIGAR, NON_ZERO_REFERENCE_LENGTH_ALIGNMENT}
// TODO: filter reads based on length value (if too large), and/or minimum_pct like in Picard.
// TODO: case EITHER for enum EndToUser. For truncated genomic regions of putative SV breakpoints, not all reads have
//       both ends land in the region, so third case is possible: will use either end when only one end is available in
//       the region specified, and only first end if both are available.
// TODO: user argument validation (eg. maxMADTolerance)
@CommandLineProgramProperties(
        summary        = "Program to collect insert size distribution information in SAM/BAM file(s)",
        oneLineSummary = "Collect Insert Size Distribution on Spark",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class CollectInsertSizeMetricsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "A local path to file to write insert size metrics to, with extension.",
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              optional = false)
    public String output = null;

    @Argument(doc = "A local path to PDF file where histogram plot will be saved in.",
              shortName = "HIST",
              fullName = "HistogramPlotPDF",
              optional = false)
    public String histogramPlotFile = null;

    @Argument(doc = "Generate mean, sd and plots by trimming the data down to MEDIAN + maxMADTolerance*MEDIAN_ABSOLUTE_DEVIATION. " +
                    "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
                    "artifacts to make the mean and sd grossly misleading regarding the real distribution.",
              shortName = "TOL",
              fullName = "HistogramPlotDeviationsTolerance",
              optional = true)
    public double maxMADTolerance = 10.0;

    // read filtering criteria
    @Argument(doc = "If set to true, filter pairs of reads that are not properly--as judged by aligner--oriented.",
              shortName = "PP",
              fullName = "filterNonProperlyPairedReads",
              optional = true)
    public boolean filterNonProperlyPairedReads = false;

    @Argument(doc = "If set to true, include duplicated reads as well.",
              shortName = "Dup",
              fullName = "useDuplicateReads",
              optional = true)
    public boolean useDuplicateReads = false;

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
              shortName = "MAPQ",
              fullName = "MAPQThreshold",
              optional = true)
    public int MQPassingThreshold = 0;

    @Argument(doc="The level(s) at which to accumulate metrics. Possible values are {ALL_READS, SAMPLE, LIBRARY, READ GROUP}.",
              shortName="LEVEL",
              fullName = "MetricsAccumulationLevel",
              optional = false)
    public Set<MetricAccumulationLevel> metricAccumulationLevel = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    /**
     * Which end of a read pair to use for collecting insert size metrics.
     */
    public enum EndToUse {
        FIRST(1), SECOND(2);
        private final int value;
        EndToUse(int value){
            this.value = value;
        }
        public int getValue(){
            return value;
        }
    }

    @Argument(doc = "Which end of pairs to use for collecting information. " +
                    "Possible values:{FIRST, SECOND}.",
              shortName = "E",
              fullName = "whichEndOfPairToUse",
              optional = true)
    public EndToUse useEnd = EndToUse.FIRST;

    // path to Picard R script for producing histograms in PDF files.
    private static final String R_SCRIPT = "insertSizeHistogram.R";

    @Override
    public boolean requiresReads(){ return true; }

    @Override
    public void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> filteredReads = getReads();

        if(filteredReads.isEmpty()) {throw new GATKException("No valid reads found in input file."); }

        final SAMFileHeader readsHeader = getHeaderForReads();

        final MetricsFile<InsertSizeMetrics, Integer> metricsFile = getMetricsFile();

        final InsertSizeMetricsCollectorSpark collector = new InsertSizeMetricsCollectorSpark(filteredReads,
                                                                                              readsHeader,
                                                                                              metricAccumulationLevel,
                                                                                              maxMADTolerance,
                                                                                              metricsFile);

        MetricsUtils.saveMetrics(metricsFile, output, getAuthHolder());
        writeHistogramPDF();
    }

    /**
     *  Implicitly called in getReads(): behavior from base (GATKTools). Return type serializable.
     *  Not calling base(GATKTools) because custom filters is not a superset of filters defined there.
     */
    @Override
    public ReadFilter makeReadFilter() {

        final SVCustomReadFilter svFilter = new SVCustomReadFilter(useEnd,
                                                                   filterNonProperlyPairedReads,
                                                                   !useDuplicateReads,
                                                                   !useSecondaryAlignments,
                                                                   !useSupplementaryAlignments,
                                                                   MQPassingThreshold);

        return new WellformedReadFilter(getHeaderForReads()).and(svFilter);
    }

    /**
     * Customized serializable reads filter, based on cmd line arguments provided
     */
    private static final class SVCustomReadFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        private final ReadFilter combinedReadFilter;

        public SVCustomReadFilter(final EndToUse whichEnd,
                                  final boolean filterNonProperlyPairedReads,
                                  final boolean filterDuplicatedReads,
                                  final boolean filterSecondaryAlignments,
                                  final boolean filterSupplementaryAlignments,
                                  final int     MQThreshold){

            final EndToUse endVal = whichEnd;

            ReadFilter tempFilter = ReadFilterLibrary.MAPPED;
                       tempFilter = tempFilter.and(GATKRead::isPaired);
                       tempFilter = tempFilter.and(read -> 0!=read.getFragmentLength());
                       tempFilter = tempFilter.and(read -> endVal == (read.isFirstOfPair() ? EndToUse.FIRST : EndToUse.SECOND));

            if(filterNonProperlyPairedReads)  { tempFilter = tempFilter.and(GATKRead::isProperlyPaired); }
            if(filterDuplicatedReads)         { tempFilter = tempFilter.and(read -> !read.isDuplicate()); }
            if(filterSecondaryAlignments)     { tempFilter = tempFilter.and(read -> !read.isSecondaryAlignment()); }
            if(filterSupplementaryAlignments) { tempFilter = tempFilter.and(read -> !read.isSupplementaryAlignment()); }

            if(0!=MQThreshold)  { tempFilter = tempFilter.and(read -> read.getMappingQuality() >= MQThreshold);}

            combinedReadFilter = tempFilter;
        }

        @Override
        public boolean test(final GATKRead read){
            return combinedReadFilter.test(read);
        }
    }

    /**
     * Calls R script to plot histogram(s) in PDF.
     * @throws RScriptExecutorException
     */
    @VisibleForTesting
    void writeHistogramPDF() throws RScriptExecutorException {

        final File histogramPlotPDF = new File(histogramPlotFile);
        IOUtil.assertFileIsWritable(histogramPlotPDF);

        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, CollectInsertSizeMetricsSpark.class));
        executor.addArgs(output,                                // text-based metrics file
                         histogramPlotPDF.getAbsolutePath(),    // PDF graphics file
                         getReadSourceName());                  // input bam file
        executor.exec();
    }
}
