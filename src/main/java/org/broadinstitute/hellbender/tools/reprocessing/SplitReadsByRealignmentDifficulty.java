package org.broadinstitute.hellbender.tools.reprocessing;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnels;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiplePassReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to two new SAM/BAM/CRAM files (unsorted).",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file (unsorted)",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public class SplitReadsByRealignmentDifficulty extends MultiplePassReadWalker {

    public static final String MAPPABILITY_INTERVAL_LIST_SHORT_NAME = "umil";
    public static final String MAPPABILITY_INTERVAL_LIST_LONG_NAME = "uniqueMapIntervals";
    public static final String READ_LOG_COUNT_SHORT_NAME = "RLG";
    public static final String READ_LOG_COUNT_LONG_NAME =  "reads_log_count";
    public static final String EXPECTED_UNCERTAIN_READS_SHORT_NAME = "EUR";
    public static final String EXPECTED_UNCERTAIN_READS_LONG_NAME = "expected_uncertain_reads";
    public static final String BLOOM_FILTER_FPR_SHORT_NAME = "FPR";
    public static final String BLOOM_FILTER_FPR_LONG_NAME = "bloom_filter_fpr";

    @Argument(fullName = MAPPABILITY_INTERVAL_LIST_LONG_NAME,
            shortName = MAPPABILITY_INTERVAL_LIST_SHORT_NAME,
            doc="Mappability interval list.  This must be a local path.")
    public GATKPath mappabilityIntervalPath;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write naive output to this file")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME+"uncertain",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME+"u",
            doc="Write uncertain output to this file")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath outputUncertainReads;

    @Argument(fullName = READ_LOG_COUNT_LONG_NAME,
            shortName = READ_LOG_COUNT_SHORT_NAME,
            doc="How many reads have been processed before logging an update")
    public int readsLogCount=1000000;

    @Argument(fullName = EXPECTED_UNCERTAIN_READS_LONG_NAME,
            shortName = EXPECTED_UNCERTAIN_READS_SHORT_NAME,
            doc="Expected number of uncertain read names (approximately 15% of total read pairs). " +
                "Used to size the Bloom filter. Overestimating is safe and uses slightly more memory; " +
                "underestimating increases false positive rate but remains correct.",
            optional = true)
    public long expectedUncertainReads = 250_000_000L;

    @Argument(fullName = BLOOM_FILTER_FPR_LONG_NAME,
            shortName = BLOOM_FILTER_FPR_SHORT_NAME,
            doc="False positive rate for the Bloom filter (0.0 to 1.0). Lower values use more memory " +
                "but result in fewer naive reads being incorrectly classified as uncertain. " +
                "Default of 0.001 (0.1%) uses ~314MB for 175M uncertain reads.",
            optional = true)
    public double bloomFilterFpr = 0.001;

    private SAMFileGATKReadWriter outputWriter;
    private SAMFileGATKReadWriter outputWriterUncertain;
    private BloomFilter<String> uncertainBloomFilter;
    private int readCountPass1 = 0;
    private int readCountPass2 = 0;
    private OverlapDetector<SimpleInterval> overlapDetector;
    private int naiveCount = 0;
    private int uncertainCount = 0;
    private long uncertainReadNamesAdded = 0;


    private void logHeapUsage(final String phase) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.info("Used memory (MB) after " + phase + ": " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
    }

    private int overlapCount(final Locatable i) {
        int result = 0;
        final Set<SimpleInterval> overlaps = overlapDetector.getOverlaps(i);
        for (SimpleInterval overlap: overlaps) {
            int tmp = overlap.size();
            if (tmp > result) {
                result = tmp;
            }
        }
        return result;
    }

    @Override
    public void traverseReads() {
        // Naive: We believe that we can naively reprocess.  Ie, alignment stays the same
        // Uncertain: We believe that we will need to reprocess this read.  Ie, alignment may change
        // Unknown:  We do not yet know whether the read is naive or uncertain.  Typically, this is when we are
        //  awaiting information on the mate read

        //  For the below table this will apply to both the read in question and the mate.
        // If a read is >100 overlapped and the mate is >100 overlapped: Naive
        // If a read is >100 overlapped and the mate is not on same contig:  Uncertain (takes precedence)
        // If a read is >100 overlapped and the mate is not overlapped: Naive
        // If a read is unmapped: uncertain
        // If a read does not overlap and mate does not overlap:  uncertain

        outputWriter = createSAMWriter(output, false);
        outputWriterUncertain = createSAMWriter(outputUncertainReads, false);

        // Initialize Bloom filter for tracking uncertain read names
        // We only track uncertain reads (minority ~15%) and default to naive for anything not in the filter
        uncertainBloomFilter = BloomFilter.create(
                Funnels.stringFunnel(StandardCharsets.UTF_8),
                expectedUncertainReads,
                bloomFilterFpr
        );
        logger.info("Initialized Bloom filter: expectedInsertions=" + expectedUncertainReads + 
                    ", fpr=" + bloomFilterFpr + 
                    ", estimated size (MB)=" + String.format("%.1f", expectedUncertainReads * (-Math.log(bloomFilterFpr) / (Math.log(2) * Math.log(2))) / 8.0 / 1024 / 1024));

        try {
            final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(mappabilityIntervalPath.toString(), new BEDCodec(), false);
            final List<BEDFeature> intervalsAsBed = bedReader.iterator().toList();
            final List<SimpleInterval> intervals = intervalsAsBed.stream().map(f -> new SimpleInterval(f.getContig(), f.getStart(), f.getEnd())).collect(Collectors.toList());
            overlapDetector = OverlapDetector.create(intervals);
        } catch (final IOException ioe) {
            throw new GATKException("Could not load bed file.", ioe);
        }

        // Step 1: Do a pass of all reads and determine which read names are uncertain
        // Only uncertain reads are added to the Bloom filter; naive reads require no storage
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                sortPrimaryReads(read)
        );

        // Step 2: Do a pass of all reads and determine where to write the read, keeping mates together
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                writeReads(read)
        );
    }

    private void sortPrimaryReads(final GATKRead read) {

        readCountPass1++;
        boolean isPrimary = (!read.isSecondaryAlignment() && !read.isSupplementaryAlignment() && !read.isUnmapped());
        boolean isReadOverlap = (overlapCount(read) > 100);
        
        // Only track uncertain reads in the Bloom filter
        // Naive reads (overlapping primary) require no storage - they're the default case
        if (!(isReadOverlap && isPrimary)) {
            // This read is uncertain - add to Bloom filter
            uncertainBloomFilter.put(read.getName());
            uncertainReadNamesAdded++;
        }
        // Naive reads (isReadOverlap && isPrimary) don't need to be tracked

        if ((readCountPass1 > 0) && ((readCountPass1 % readsLogCount) == 0)) {
            logReadCounts(readCountPass1, "Sorting primary reads (" + readCountPass1 + " reads)");
        }
    }

    private void writeReads(final GATKRead read) {
        readCountPass2++;
        // Write this read to the proper location. Keeping those with the same name in the same output file.
        // Use Bloom filter to check if read was classified as uncertain.
        // If in Bloom filter -> uncertain (may include some false positives from naive, which is safe)
        // If not in Bloom filter -> definitely naive (no false negatives)
        if (uncertainBloomFilter.mightContain(read.getName())) {
            outputWriterUncertain.addRead(read);
            uncertainCount++;
        } else {
            outputWriter.addRead(read);
            naiveCount++;
        }

        if (((readCountPass2 % readsLogCount) == 0) && (readCountPass2 > 0)) {
            logHeapUsage("Writing reads (" + readCountPass2 + " reads)");
        }
    }

    private void logReadCounts(int readCount, String phase) {
        logger.info("Reads processed: " + readCount);
        logger.info("Uncertain read names added to Bloom filter: " + uncertainReadNamesAdded);
        logHeapUsage(phase);
    }

    @Override
    public void closeTool() {
        logReadCounts(readCountPass2, "finish");
        logger.info("Naive reprocessing reads: " + naiveCount);
        logger.info("Uncertain reprocessing reads: " + uncertainCount);
        logger.info("Uncertain read names in Bloom filter: " + uncertainReadNamesAdded);
        
        // Estimate false positive impact
        long estimatedFalsePositives = (long) (naiveCount * bloomFilterFpr);
        logger.info("Estimated false positives (naive reads sent to uncertain due to Bloom filter FPR): ~" + estimatedFalsePositives);

        if ( outputWriter != null ) {
            outputWriter.close();
        }
        if ( outputWriterUncertain != null ) {
            outputWriterUncertain.close();
        }
    }
}
