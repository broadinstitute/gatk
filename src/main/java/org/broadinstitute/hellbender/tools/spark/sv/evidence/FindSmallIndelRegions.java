package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.netflix.servo.util.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

@CommandLineProgramProperties(summary="Find regions likely to contain small indels due to fragment length anomalies.",
        oneLineSummary="Find regions containing small indels.",
        programGroup=StructuralVariationSparkProgramGroup.class)
public final class FindSmallIndelRegions extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private static final FindBreakpointEvidenceSparkArgumentCollection params =
                new FindBreakpointEvidenceSparkArgumentCollection();

    private static final int EVIDENCE_SIZE_GUESS = 1000;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final JavaRDD<GATKRead> allReads = getUnfilteredReads();
        final Set<Integer> crossContigIgnoreSet = Collections.emptySet();
        final SAMFileHeader header = getHeaderForReads();
        final int maxTrackedFragmentLength = params.maxTrackedFragmentLength;
        final SVReadFilter filter = new SVReadFilter(params);
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigIgnoreSet, header, maxTrackedFragmentLength, allReads, filter);
        final Broadcast<ReadMetadata> broadcastReadMetadata = ctx.broadcast(readMetadata);

        final List<BreakpointEvidence> evidenceList =
            allReads
                .mapPartitions(readItr -> {
                    final List<BreakpointEvidence> evList = new ArrayList<>(EVIDENCE_SIZE_GUESS);
                    final Finder finder = new Finder(broadcastReadMetadata.value(), filter);
                    while ( readItr.hasNext() ) {
                        finder.test(readItr.next(), evList);
                    }
                    finder.checkHistograms(evList);
                    return evList.iterator();
                })
                .collect();

        System.out.println();
        final List<SAMSequenceRecord> seqs = header.getSequenceDictionary().getSequences();
        for ( final BreakpointEvidence evidence : evidenceList ) {
            final SVInterval interval = evidence.getLocation();
            final String contigName = seqs.get(interval.getContig()).getSequenceName();
            System.out.println(contigName + "\t" + (interval.getStart() + 1) + "\t" + interval.getEnd());
        }
    }

    public static final class Finder {
        private static final int BLOCK_SIZE = 250;
        private static final int EVIDENCE_WEIGHT = 5;
        private static final float KS_SIGNIFICANCE = 1.E-5f;

        private final ReadMetadata readMetadata;
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
        private final Map<String, IntHistogram[]> libraryToHistoPairMap;
        private int fillIdx;
        // curContig+curEnd mark the genomic location of the end of the current window
        private String curContig;
        private int curEnd;

        public Finder( final ReadMetadata readMetadata, final SVReadFilter filter ) {
            this.readMetadata = readMetadata;
            this.filter = filter;
            final Map<String, FragmentLengthStatistics> libraryStatisticsMap = readMetadata.getAllLibraryStatistics();
            libraryToHistoPairMap = new HashMap<>(SVUtils.hashMapCapacity(libraryStatisticsMap.size()));
            libraryStatisticsMap.forEach( (libName, stats) -> {
                final IntHistogram[] pair = new IntHistogram[2];
                pair[0] = stats.createEmptyHistogram();
                pair[1] = stats.createEmptyHistogram();
                libraryToHistoPairMap.put(libName, pair);
            });
            curContig = null;
            curEnd = 0;
            fillIdx = 0;
        }

        public void test( final GATKRead read, final List<BreakpointEvidence> evidenceList ) {
            if ( !filter.isTemplateLenTestable(read) ) return;

            final boolean sameContig;
            if ( !(sameContig = read.getContig().equals(curContig)) || read.getStart() >= curEnd ) {
                checkHistograms(evidenceList);
                if ( sameContig ) {
                    curEnd += BLOCK_SIZE;
                } else {
                    curContig = read.getContig();
                    curEnd = read.getStart() + BLOCK_SIZE;
                    // clear 1st-half window when we switch contigs
                    final int oldIdx = fillIdx ^ 1;
                    libraryToHistoPairMap.values().forEach(histPair -> histPair[oldIdx].clear());
                }
            }
            final IntHistogram[] histoPair =
                    libraryToHistoPairMap.get(readMetadata.getLibraryName(read.getReadGroup()));
            histoPair[fillIdx].addObservation(Math.abs(read.getFragmentLength()));
        }

        public void checkHistograms( final List<BreakpointEvidence> evidenceList ) {
            final int oldIdx = fillIdx ^ 1; // bit magic sets oldIdx to 1 if fillIdx is 0, and vice versa
            SVInterval curInterval = null;
            for ( final Map.Entry<String,IntHistogram[]> entry : libraryToHistoPairMap.entrySet() ) {
                final IntHistogram[] histoPair = entry.getValue();
                final IntHistogram oldHisto = histoPair[oldIdx];
                oldHisto.addObservations(histoPair[fillIdx]);
                final long readCount = oldHisto.getTotalObservations();
                if ( readCount > 0 &&
                        readMetadata.getLibraryStatistics(entry.getKey())
                                .isDifferentByKSStatistic(oldHisto,KS_SIGNIFICANCE) ) {
                    if ( curInterval == null ) {
                        int start = curEnd - 2 * BLOCK_SIZE;
                        if ( start < 1 ) start = 1;
                        curInterval = new SVInterval(readMetadata.getContigID(curContig), start, curEnd);
                    }
                    evidenceList.add(
                            new BreakpointEvidence.TemplateSizeAnomaly(curInterval,EVIDENCE_WEIGHT,(int)readCount));
                }
                oldHisto.clear();
            }
            fillIdx = oldIdx;
        }

        @VisibleForTesting
        Map<String, IntHistogram[]> getLibraryToHistoPairMap() { return libraryToHistoPairMap; }
    }
}
