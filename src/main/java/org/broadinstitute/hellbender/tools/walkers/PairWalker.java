package org.broadinstitute.hellbender.tools.walkers;

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
        summary = "Prints read pairs within a region, along with their mates.",
        oneLineSummary = "Print read pairs within a region and their mates.",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public abstract class PairWalker extends ReadWalker {

    @Argument(fullName = "pair-padding", shortName = "pp",
                doc = "Amount of extra padding (in bp) to add to pick up mates near your intervals.")
    protected int intervalPadding = 1000;

    private RegionChecker regionChecker = null;
    private final HopscotchSet<PairBuffer> pairBufferSet = new HopscotchSet<>(1000000);

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
        final boolean inInterval = regionChecker == null || regionChecker.isInInterval(read);
        if ( read.isPaired() ) {
            final PairBuffer pb = pairBufferSet.findOrAdd(new PairBuffer(read.getName()), k -> (PairBuffer)k);
            pb.add(read, inInterval);
            if ( pb.isComplete() ) {
                if ( pb.isApplicable(regionChecker) ) {
                    apply(pb.getRead1(), pb.getRead2());
                }
                pairBufferSet.remove(pb);
            }
        } else if ( inInterval ) {
            applyUnpaired(read);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        int nUnpaired = 0;
        for ( final PairBuffer pb : pairBufferSet ) {
            if ( pb.isInInterval() ) {
                final GATKRead read = pb.getRead1();
                if ( read != null ) applyUnpaired(read);
                else applyUnpaired(pb.getRead2());
                nUnpaired += 1;
            }
        }
        if ( nUnpaired > 0 ) {
            logger.info("There were " + nUnpaired + " incomplete pairs.");
        }
        return null;
    }

    public abstract void apply( final GATKRead read, final GATKRead mate );

    public abstract void applyUnpaired( final GATKRead read );

    private static final class RegionChecker {
        private final List<SimpleInterval> intervals;
        private final SAMSequenceDictionary dictionary;
        private Iterator<SimpleInterval> intervalIterator;
        private SimpleInterval currentInterval;

        public RegionChecker( final List<SimpleInterval> intervals, final SAMSequenceDictionary dictionary ) {
            this.intervals = intervals;
            this.dictionary = dictionary;
            reset();
        }

        public boolean isInCurrentInterval( final GATKRead read ) {
            return currentInterval.overlaps(read);
        }

        public boolean isInInterval( final GATKRead read ) {
            if ( currentInterval == null || read.isUnmapped() ) return false;
            final String readContig = read.getContig();
            if ( !currentInterval.getContig().equals(readContig) ) {
                final int readContigId = dictionary.getSequenceIndex(readContig);
                int currentContigId;
                while ( (currentContigId = dictionary.getSequenceIndex(currentInterval.getContig())) < readContigId ) {
                    if ( intervalIterator.hasNext() ) {
                        currentInterval = intervalIterator.next();
                    } else {
                        currentInterval = null;
                        return false;
                    }
                }
                if ( currentContigId > readContigId ) {
                    return false;
                }
            }
            final int readStart = read.getStart();
            do {
                if ( currentInterval.getEnd() < readStart ) {
                    if ( intervalIterator.hasNext() ) {
                        currentInterval = intervalIterator.next();
                    } else {
                        currentInterval = null;
                        return false;
                    }
                } else {
                    return currentInterval.overlaps(read);
                }
            } while ( currentInterval.getContig().equals(readContig) );
            return false;
        }

        public void reset() {
            intervalIterator = intervals.iterator();
            currentInterval = intervalIterator.hasNext() ? intervalIterator.next() : null;
        }
    }

    private static final class PairBuffer {
        private final String qName;
        private GATKRead read1;
        private GATKRead read2;
        private boolean inInterval;
        private boolean distantMate1;
        private boolean distantMate2;

        public PairBuffer( final String qName ) {
            this.qName = qName;
            read1 = null;
            read2 = null;
            inInterval = false;
            distantMate1 = false;
            distantMate2 = false;
        }

        public boolean isComplete() { return read1 != null && read2 != null; }
        public boolean isInInterval() { return inInterval; }
        public GATKRead getRead1() { return read1; }
        public GATKRead getRead2() { return read2; }

        public void add( final GATKRead readArg, final boolean inInterval ) {
            final GATKRead read = PrintDistantMates.untangleRead(readArg);
            if ( read.isFirstOfPair() ) {
                read1 = read;
                if ( read != readArg ) distantMate1 = true;
            }
            else {
                read2 = read;
                if ( read != readArg ) distantMate2 = true;
            }
            this.inInterval |= inInterval;
        }

        public boolean isApplicable( final RegionChecker regionChecker ) {
            return inInterval && (!(distantMate1 || distantMate2) || !isDuplicative(regionChecker));
        }

        private boolean isDuplicative( final RegionChecker regionChecker ) {
            if ( regionChecker == null ||
                    (regionChecker.isInCurrentInterval(read1) && regionChecker.isInCurrentInterval(read2)) ) {
                return read1.getStart() < read2.getStart() ? distantMate1 : distantMate2;
            }
            return false;
        }

        @Override public boolean equals( final Object obj ) {
            return this == obj || obj instanceof PairBuffer && qName.equals(((PairBuffer)obj).qName);
        }

        @Override public int hashCode() { return qName.hashCode(); }
    }
}
