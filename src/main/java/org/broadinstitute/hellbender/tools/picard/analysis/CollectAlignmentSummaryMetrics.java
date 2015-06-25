package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.utils.illumina.IlluminaAdapterPair;

import java.io.File;
import java.util.*;

/**
 * A command line tool to read a BAM file and produce standard alignment metrics that would be applicable to any alignment.  
 * Metrics to include, but not limited to:
 * <ul>
 * <li>Total number of reads (total, period, no exclusions)</li>
 * <li>Total number of PF reads (PF == does not fail vendor check flag)</li>
 * <li>Number of PF noise reads (does not fail vendor check and has noise attr set)</li>
 * <li>Total aligned PF reads (any PF read that has a sequence and position)</li>
 * <li>High quality aligned PF reads (high quality == mapping quality >= 20)</li>
 * <li>High quality aligned PF bases (actual aligned bases, calculate off alignment blocks)</li>
 * <li>High quality aligned PF Q20 bases (subset of above where base quality >= 20)</li>
 * <li>Median mismatches in HQ aligned PF reads (how many aligned bases != ref on average)</li>
 * <li>Reads aligned in pairs (vs. reads aligned with mate unaligned/not present)</li>
 * <li>Read length (how to handle mixed lengths?)</li>
 * <li>Bad Cycles - how many machine cycles yielded combined no-call and mismatch rates of >= 80%</li>
 * <li>Strand balance - reads mapped to positive strand / total mapped reads</li>
 * </ul>
 * Metrics are written for the first read of a pair, the second read, and combined for the pair.
 * 
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
@CommandLineProgramProperties(
        usage = CollectAlignmentSummaryMetrics.USAGE,
        usageShort = CollectAlignmentSummaryMetrics.USAGE,
        programGroup = QCProgramGroup.class
)
public final class CollectAlignmentSummaryMetrics extends SinglePassSamProgram {
    static final String USAGE = "Produces from a SAM or BAM a file containing summary alignment metrics";
    
    private static final Log log = Log.getInstance(CollectAlignmentSummaryMetrics.class);

    // Usage and parameters

    @Argument(doc="Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.")
    public int MAX_INSERT_SIZE = 100000;

    @Argument(doc="List of adapter sequences to use when processing the alignment metrics")
	public List<String> ADAPTER_SEQUENCE = CollectionUtil.makeList(
            IlluminaAdapterPair.SINGLE_END.get5PrimeAdapter(),
            IlluminaAdapterPair.SINGLE_END.get3PrimeAdapter(),
            IlluminaAdapterPair.PAIRED_END.get5PrimeAdapter(),
            IlluminaAdapterPair.PAIRED_END.get3PrimeAdapter(),
            IlluminaAdapterPair.INDEXED.get5PrimeAdapter(),
            IlluminaAdapterPair.INDEXED.get3PrimeAdapter()
    );

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(shortName="BS", doc="Whether the SAM or BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    private AlignmentSummaryMetricsCollector collector;

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);

        if (header.getSequenceDictionary().isEmpty()) {
            log.warn(INPUT.getAbsoluteFile() + " has no sequence dictionary.  If any reads " +
                    "in the file are aligned then alignment summary metrics collection will fail.");
        }

        final boolean doRefMetrics = REFERENCE_SEQUENCE != null;
        collector = new AlignmentSummaryMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), doRefMetrics,
                ADAPTER_SEQUENCE, MAX_INSERT_SIZE,  IS_BISULFITE_SEQUENCED);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        collector.acceptRecord(rec, ref);
    }

    @Override
    protected void finish() {
        collector.finish();

        final MetricsFile<AlignmentSummaryMetrics, Long> file = getMetricsFile();
        collector.addAllLevelsToFile(file);

        file.write(OUTPUT);
    }
}
