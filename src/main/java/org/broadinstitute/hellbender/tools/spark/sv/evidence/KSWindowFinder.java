package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class KSWindowFinder {
    protected static final int BLOCK_SIZE = 250;
    protected static final float KS_SIGNIFICANCE = 3.1E-4f;
    private static final int EVIDENCE_WEIGHT = 5;

    // When multiplied by the library read frequency this number gives an upper bound on the number of observations
    // in a window.  we don't want to evaluate windows with too many observations because they're mostly due
    // to mapping defects causing huge pileups.
    // Here's how this number is actually derived:
    // WINDOW_SIZE = 2 * BLOCK_SIZE
    // EXPECTED_READS_PER_WINDOW = WINDOW_SIZE * library read frequency
    // Deflate according to a heuristic derived from the shape of the actual distribution -- due to mapping
    // error there's a shift down in the median from all the huge pileups sucking reads away.
    // This may need to be revisited whenever we change the mapper.
    // ACTUAL_MEDIAN_READS_PER_WINDOW = EXPECTED_READS_PER_WINDOW / 2
    // Finally, we put, as a maximum, a copy number of 5 (most things with copy numbers that high are mapping
    // artifacts).  that's 2.5x the expected number of reads because that's a diploid number.
    protected static final float MAX_APPARENT_BASES_PER_WINDOW = BLOCK_SIZE * 2.5f;

    protected final ReadMetadata readMetadata;
    private final SVReadFilter filter;

    // We're calculating a K-S statistic over 500-base windows, but the windows are tiled every 250 bases.
    // So what we do is to keep a pair of histograms for each library: one for first half of the window, and one
    // for the second half.
    // The reads are coordinate sorted, and whenever the next read presented to the test method crosses out of the
    // current window, it's time to do a K-S test.
    // We add the 2nd-half histogram to the 1st-half histogram, so that the 1st-half histogram now has the complete
    // data for the entire window.  We do the K-S test, and then zero the 1st-half histogram.
    // At this point the 2nd-half histogram becomes the 1st-half histogram for the next window, and the newly
    // cleared 1st-half histogram becomes the 2nd-half histogram for the next window.
    // We swap roles simply by toggling fillIdx between 0 and 1 each time we cross a boundary.
    // So in the code that follows, fillIdx points to the member of the pair that is currently accumulating data --
    // i.e., the 2nd-half histogram.  The temporary variable oldIdx in the checkHistograms method below points to
    // the other member of the pair -- i.e., the 1st-half histogram.
    protected final Map<String, IntHistogram[]> libraryToHistoPairMap;
    protected int fillIdx;
    // curContig+curEnd mark the genomic location of the end of the current window
    protected String curContig;
    protected int curEnd;

    public KSWindowFinder( final ReadMetadata readMetadata, final SVReadFilter filter ) {
        this.readMetadata = readMetadata;
        this.filter = filter;
        final Map<String, LibraryStatistics> libraryStatisticsMap = readMetadata.getAllLibraryStatistics();
        libraryToHistoPairMap = new HashMap<>(SVUtils.hashMapCapacity(libraryStatisticsMap.size()));
        libraryStatisticsMap.forEach(( libName, stats ) -> {
            final IntHistogram[] pair = new IntHistogram[2];
            pair[0] = stats.createEmptyHistogram();
            pair[1] = stats.createEmptyHistogram();
            libraryToHistoPairMap.put(libName, pair);
        });
        curContig = null;
        curEnd = 0;
        fillIdx = 0;
    }

    /**
     * Used in production to gather evidence derived from K-S calculation on windows.
     */
    public void testReadAndGatherEvidence( final GATKRead read, final List<BreakpointEvidence> evidenceList ) {
        if ( !filter.isTemplateLenTestable(read) ) return;

        if ( !isSameBlock(read) ) {
            checkHistograms(evidenceList);
            advanceBlock(read);
        }
        addObservation(read);
    }

    /**
     * Used in production to gather evidence derived from K-S calculation on windows.
     */
    public void checkHistograms( final List<BreakpointEvidence> evidenceList ) {
        if ( curContig == null ) return;

        final int oldIdx = fillIdx ^ 1; // bit magic sets oldIdx to 1 if fillIdx is 0, and vice versa
        final int start = Math.max(1, curEnd - 2 * BLOCK_SIZE);
        final SVInterval curInterval = new SVInterval(readMetadata.getContigID(curContig), start, curEnd);
        for ( final Map.Entry<String, IntHistogram[]> entry : libraryToHistoPairMap.entrySet() ) {
            final IntHistogram[] histoPair = entry.getValue();
            final IntHistogram oldHisto = histoPair[oldIdx];
            oldHisto.addObservations(histoPair[fillIdx]);
            final long readCount = oldHisto.getTotalObservations();
            if ( readCount > 0 &&
                    readMetadata.getLibraryStatistics(entry.getKey()).getCDF()
                            .isDifferentByKSStatistic(oldHisto, KS_SIGNIFICANCE) ) {
                evidenceList.add(
                        new BreakpointEvidence.TemplateSizeAnomaly(curInterval, EVIDENCE_WEIGHT, (int)readCount));
            }
            oldHisto.clear();
        }
    }

    @VisibleForTesting
    Map<String, IntHistogram[]> getLibraryToHistoPairMap() { return libraryToHistoPairMap; }

    protected boolean isTestable( final GATKRead read ) { return filter.isTemplateLenTestable(read); }

    protected boolean isSameBlock( final GATKRead read ) {
        return read.getContig().equals(curContig) && read.getStart() < curEnd;
    }

    protected void addObservation( final GATKRead read ) {
        final IntHistogram[] histoPair = libraryToHistoPairMap.get(readMetadata.getLibraryName(read.getReadGroup()));
        histoPair[fillIdx].addObservation(Math.abs(read.getFragmentLength()));
    }

    protected void advanceBlock( final GATKRead read ) {
        final int oldIdx = fillIdx;
        fillIdx ^= 1; // bit magic to switch halves
        curEnd += BLOCK_SIZE;
        if ( !read.getContig().equals(curContig) || read.getStart() >= curEnd ) {
            curContig = read.getContig();
            curEnd = read.getStart() + BLOCK_SIZE;
            // clear 1st-half window when we switch contigs or blast forward due to a gap in coverage
            libraryToHistoPairMap.values().forEach(histPair -> histPair[oldIdx].clear());
        }
    }
}