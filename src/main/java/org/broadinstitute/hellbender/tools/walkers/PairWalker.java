package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.PrintDistantMates;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.HopscotchSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

@CommandLineProgramProperties(
        summary = "A walker that presents read pairs as pairs.  Run PrintDistantMates, and use "+
                "its output as input to the PairWalker along with the original inputs to ensure "+
                "seeing all pairs as pairs.  This is unnecessary if you only care about pairs "+
                "where both reads lie within your intervals.",
        oneLineSummary = "A walker that presents read pairs as pairs.",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public abstract class PairWalker extends ReadWalker {

    @Argument(fullName = "pair-padding", shortName = "pp",
                doc = "Amount of extra padding (in bp) to add to pick up mates near your intervals.")
    protected int intervalPadding = 1000;

    private RegionChecker regionChecker = null;
    private final HopscotchSet<PairBufferEntry> pairBufferSet = new HopscotchSet<>(1000000);
    private final HopscotchSet<String> distantPairsProcessed = new HopscotchSet<>(1000000);

    @Override public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
        readFilters.add(ReadFilterLibrary.PRIMARY_LINE);
        readFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        return readFilters;
    }

    @Override
    public boolean requiresReads() { return true; }

    @Override
    protected List<SimpleInterval> transformTraversalIntervals( final List<SimpleInterval> intervals,
                                                                final SAMSequenceDictionary dictionary ) {
        regionChecker = new RegionChecker(intervals, dictionary);

        final List<SimpleInterval> paddedIntervals = new ArrayList<>(intervals.size());
        SimpleInterval prevInterval = null;
        for ( final SimpleInterval interval : intervals ) {
            final SimpleInterval curInterval = interval.expandWithinContig(intervalPadding, dictionary);
            if ( prevInterval == null ) {
                prevInterval = curInterval;
            } else if ( prevInterval.getContig().equals(curInterval.getContig()) &&
                        prevInterval.getStart() <= curInterval.getEnd() + 1 &&
                        prevInterval.getStart() <= curInterval.getEnd() + 1 ) {
                prevInterval = prevInterval.mergeWithContiguous(curInterval);
            } else {
                paddedIntervals.add(prevInterval);
                prevInterval = curInterval;
            }
        }
        if ( prevInterval != null ) {
            paddedIntervals.add(prevInterval);
        }
        return Collections.unmodifiableList(paddedIntervals);
    }

    @Override
    public void apply( final GATKRead read,
                       final ReferenceContext referenceContext,
                       final FeatureContext featureContext ) {
        if ( !read.isPaired() || read.isSecondaryAlignment() || read.isSupplementaryAlignment() ) {
            applyUnpaired(read);
            return;
        }
        final boolean inInterval = regionChecker == null || regionChecker.isInInterval(read);
        final PairBufferEntry pb = new PairBufferEntry(read, inInterval);
        final PairBufferEntry pbMate = pairBufferSet.find(pb);
        if ( pbMate == null ) {
            pairBufferSet.add(pb);
        } else {
            if ( pb.isInInterval() || pbMate.isInInterval() ) {
                if ( !pb.isDistantMate() && !pbMate.isDistantMate() ) {
                    apply(pbMate.getRead(), pb.getRead());
                } else {
                    final String readName = pb.getRead().getName();
                    if ( distantPairsProcessed.contains(readName) ) {
                        distantPairsProcessed.remove(readName);
                    } else {
                        distantPairsProcessed.add(readName);
                        apply(pbMate.getRead(), pb.getRead());
                    }
                }
            }
            pairBufferSet.remove(pbMate);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        int nUnpaired = 0;
        for ( final PairBufferEntry pb : pairBufferSet ) {
            if ( pb.isInInterval() ) {
                applyUnpaired(pb.getRead());
                nUnpaired += 1;
            }
        }
        if ( nUnpaired > 0 ) {
            logger.info("There were " + nUnpaired + " incomplete pairs.");
        }
        return null;
    }

    /**
     * This method is called once for each pair of reads.
     * The pairs will NOT be in strict coordinate order.
     * The reads supplied will be the primary lines for each pair.
     **/
    public abstract void apply( final GATKRead read, final GATKRead mate );

    /**
     * Unpaired reads, secondary and supplemental alignments, and other detritus comes through here.
     * Also, if you haven't used PrintDistantMates, you'll get all of pairs for which we didn't
     * find the mate at the end of the traversal via this method.
     */
    public abstract void applyUnpaired( final GATKRead read );

    // Maintains an iterator over a set of coordinate-sorted, disjoint intervals.
    // It allows you to ask whether a read overlaps any of the intervals by calling isInInterval.
    // With each call, it advances the iterator as appropriate, and checks for overlap.
    // The sequence of reads passed to isInInterval must also be coordinate-sorted.
    @VisibleForTesting
    static final class RegionChecker {
        private final List<SimpleInterval> intervals;
        private final SAMSequenceDictionary dictionary;
        private Iterator<SimpleInterval> intervalIterator;
        private SimpleInterval currentInterval;

        public RegionChecker( final List<SimpleInterval> intervals, final SAMSequenceDictionary dictionary ) {
            this.intervals = intervals;
            this.dictionary = dictionary;
            reset();
        }

        public boolean isInInterval( final GATKRead read ) {
            if ( currentInterval == null || read.isUnmapped() ) return false;
            final String readContig = read.getContig();
            if ( !currentInterval.getContig().equals(readContig) ) {
                final int readContigId = dictionary.getSequenceIndex(readContig);
                int currentContigId;
                while ( (currentContigId =
                        dictionary.getSequenceIndex(currentInterval.getContig())) < readContigId ) {
                    if ( !intervalIterator.hasNext() ) {
                        currentInterval = null;
                        return false;
                    }
                    currentInterval = intervalIterator.next();
                }
                if ( currentContigId > readContigId ) {
                    return false;
                }
            }

            // we've advanced the iterator so that the current contig is the same as the read's contig
            final int readStart = read.getStart();
            while ( currentInterval.getEnd() < readStart ) {
                if ( !intervalIterator.hasNext() ) {
                    currentInterval = null;
                    return false;
                }
                currentInterval = intervalIterator.next();
                if ( !currentInterval.getContig().equals(readContig) ) {
                    return false;
                }
            }

            // we've advanced the iterator so that the current end is no smaller than the read start
            return currentInterval.overlaps(read);
        }

        public void reset() {
            intervalIterator = intervals.iterator();
            currentInterval = intervalIterator.hasNext() ? intervalIterator.next() : null;
        }
    }

    // Tracks the first read of a pair we observe, while we're waiting for the mate to show up.
    // Note that it's "keyed" on read name:  Two PairBufferEntries are equal if they point to reads
    // with the same name.
    private static final class PairBufferEntry {
        private final GATKRead read;
        private final boolean inInterval;
        private final boolean distantMate;

        public PairBufferEntry( final GATKRead read, final boolean inInterval ) {
            this.distantMate = PrintDistantMates.isDistantMate(read);
            this.read = this.distantMate ? PrintDistantMates.undoDistantMateAlterations(read) : read;
            this.inInterval = inInterval;
        }

        public GATKRead getRead() { return read; }
        public boolean isInInterval() { return inInterval; }
        public boolean isDistantMate() { return distantMate; }

        @Override
        public boolean equals( final Object obj ) {
            return this == obj ||
                    (obj instanceof PairBufferEntry &&
                            ((PairBufferEntry)obj).read.getName().equals(read.getName()));
        }

        @Override
        public int hashCode() {
            return read.getName().hashCode();
        }
    }
}
