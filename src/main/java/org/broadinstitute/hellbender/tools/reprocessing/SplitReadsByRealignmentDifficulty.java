package org.broadinstitute.hellbender.tools.reprocessing;

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
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashCode;

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

    private SAMFileGATKReadWriter outputWriter;
    private SAMFileGATKReadWriter outputWriterUncertain;
    private HashedStringSet128 overlapReadNameSet = new HashedStringSet128(4000L, 0.75f);
    private HashedStringSet128 noOverlapReadNameSet = new HashedStringSet128(4000L, 0.75f);
    private int readCountPass1 = 0;
    private int readCountPass2 = 0;
    private OverlapDetector<SimpleInterval> overlapDetector;
    private int naiveCount = 0;
    private int uncertainCount = 0;
    private int unknownReads = 0;
    private int unknownReadLogCount = 0;


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

        try {
            final FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(mappabilityIntervalPath.toString(), new BEDCodec(), false);
            final List<BEDFeature> intervalsAsBed = bedReader.iterator().toList();
            final List<SimpleInterval> intervals = intervalsAsBed.stream().map(f -> new SimpleInterval(f.getContig(), f.getStart(), f.getEnd())).collect(Collectors.toList());
            overlapDetector = OverlapDetector.create(intervals);
        } catch (final IOException ioe) {
            throw new GATKException("Could not load bed file.", ioe);
        }

        // Step 1:  Do a pass of all reads and determine which read names are primary and overlap (or not)
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                sortPrimaryReads(read)
        );

        //TODO: Check intersection of the two sets....

        // Step 2:  Do a pass of all reads and determine where to write the read, keeping mates together
        forEachRead((GATKRead read, ReferenceContext reference, FeatureContext features) ->
                writeReads(read)
        );
    }

    private void sortPrimaryReads(final GATKRead read) {

        readCountPass1++;
        boolean isPrimary = (!read.isSecondaryAlignment() && !read.isSupplementaryAlignment() && !read.isUnmapped());
        boolean isReadOverlap = (overlapCount(read) > 100);
        if (isReadOverlap && isPrimary) {
            overlapReadNameSet.add(read.getName());
        } else {
            // NOTE:  non-primary reads can be put in here and will have the same name as overlapping reads.  The
            //  reason this still works is that we check the overlapping set first in the writing phase.
            noOverlapReadNameSet.add(read.getName());
        }

        if ((readCountPass1 > 0) && ((readCountPass1 % readsLogCount) == 0)) {
            logReadCounts(readCountPass1, "Sorting primary reads (" + readCountPass1 + " reads)");
        }
    }

    private void writeReads(final GATKRead read) {
        readCountPass2++;
        // Write this read to the proper location.  Keeping those with the same name in the same output file.
        if (overlapReadNameSet.contains(read.getName())) {
            outputWriter.addRead(read);
            naiveCount++;
        } else if (noOverlapReadNameSet.contains(read.getName())) {
            outputWriterUncertain.addRead(read);
            uncertainCount++;
        } else {
            unknownReads++;
            unknownReadLogCount ++;
            if (unknownReadLogCount <= 10) {
                logger.warn("Found a read that could not be placed (only reporting the first 10): " + read);
            }
        }

        if (((readCountPass2 % readsLogCount) == 0) && (readCountPass2 > 0)) {
            logHeapUsage("Writing reads (" + readCountPass2 + " reads)");
        }
    }

    private void logReadCounts(int readCount, String phase) {
        logger.info("Reads processed: " + readCount);
        logger.info("overlap set: " + overlapReadNameSet.size());
        logger.info("no overlap set: " + noOverlapReadNameSet.size());
        logger.info("unknown (should be zero): " + unknownReads);
        logHeapUsage(phase);
    }

    @Override
    public void closeTool() {
        //TODO: Add percentages
        logReadCounts(readCountPass2, "finish");
        logger.info("Naive reprocessing reads: " + naiveCount);
        logger.info("Uncertain reprocessing reads: " + uncertainCount);

        if ( outputWriter != null ) {
            outputWriter.close();
        }
        if ( outputWriterUncertain != null ) {
            outputWriterUncertain.close();
        }
    }


    /** Exact 128-bit membership using two primitive longs per key.  Mostly AI generated*/
    // TODO: Really only checks the low 64 bits
    public final static class HashedStringSet128 {
        private final Long2LongOpenHashMap longOpenHashMap;

        public HashedStringSet128(long expectedSize, float loadFactor) {
            int initial = expectedSize > Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)expectedSize;
            this.longOpenHashMap = new Long2LongOpenHashMap(initial, loadFactor);
            this.longOpenHashMap.defaultReturnValue(Long.MIN_VALUE); // sentinel
        }

        public boolean add(String readName) {
            HashCode hc = Hashing.murmur3_128(0).hashString(readName, StandardCharsets.UTF_8);
            long low  = hc.asLong();            // low 64
            long high = toHigh64(hc.asBytes()); // high 64

//            long prev = longOpenHashMap.putIfAbsent(low, high);
            long existing = longOpenHashMap.get(low);
            if (existing == Long.MIN_VALUE) {
                // not present, safe to insert
                longOpenHashMap.put(low, high);
                // new entry
            } else if (existing == high) {
                // already present with same high bits
                // do nothing
            } else {
                // same low64 but different high64 → extremely rare collision
                throw new GATKException.ShouldNeverReachHereException("128-bit hash collision detected on read names.");
            }
            if (existing == Long.MIN_VALUE) return true;        // brand new
            if (existing == high)            return false;      // already present
            // Extremely rare: two different 128-bit values sharing low64.
            // Handle by upgrading to a tiny side structure if you truly need multiple per low64.
            throw new IllegalStateException("Low64 collision encountered; consider a tiny multimap fallback.");
        }

        public boolean contains(String readName) {
            HashCode hc = Hashing.murmur3_128(0).hashString(readName, StandardCharsets.UTF_8);
            long low  = hc.asLong();
            long high = toHigh64(hc.asBytes());
            return longOpenHashMap.get(low) == high;
        }

        public long size() { return longOpenHashMap.size(); }

        private static long toHigh64(byte[] b128) {
            // bytes[0..15] are the 128 bits; Guava orders low/high consistently.
            // Reconstruct the OTHER 64 bits:
            long x = 0;
            for (int i = 8; i < 16; i++) x = (x << 8) | (b128[i] & 0xffL);
            return x;
        }
    }
}
