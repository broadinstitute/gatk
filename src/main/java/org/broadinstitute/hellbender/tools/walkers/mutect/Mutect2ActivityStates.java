package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.engine.AlignmentAndReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.AssemblyRegionIterator;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;

/**
 * Mutect2's activity states of the loci of a shard's intervals, decided without a pileup at each locus.  Each read is walked
 * once, when it arrives, into a {@link LocusEvidence}, and each locus is decided by {@link Mutect2Engine#evaluateLocus} once
 * every read starting at or before it has arrived, as a pileup would be.  {@link AssemblyRegionIterator} forms the regions.
 *
 * A read's evidence is recorded when it arrives, before any region containing it is processed, so qualities that
 * processing a region changes in place, such as those of overlapping mates, don't affect later loci.
 */
final class Mutect2ActivityStates implements AssemblyRegionIterator.ActivityStates {
    private final PeekableIterator<GATKRead> reads;
    private final Iterator<SimpleInterval> intervals;
    private final SAMFileHeader header;
    private final ReferenceDataSource reference;
    private final FeatureManager features;
    private final Mutect2Engine engine;
    private final LocusEvidence evidence;

    // the interval whose loci we are deciding, a reference context on it and the next locus to decide
    private SimpleInterval interval = null;
    private ReferenceContext intervalReference = null;
    private int nextLocus = 0;

    // for pileup detection, which needs the pileups of each region's loci: the activity states of the loci not before the
    // region that was filled last, or null without pileup detection
    private final Queue<ActivityProfileState> pendingStates;

    /**
     * @param reads the reads of the shard, filtered, transformed and downsampled
     * @param trackPileups whether to attach each region's pileups to it, as pileup detection needs
     */
    Mutect2ActivityStates(final Iterator<GATKRead> reads, final List<SimpleInterval> intervals, final SAMFileHeader header,
                          final ReferenceDataSource reference, final FeatureManager features, final Mutect2Engine engine,
                          final boolean trackPileups) {
        Utils.validateArg(!intervals.isEmpty(), "the shard has no intervals");
        this.reads = new PeekableIterator<>(reads);
        this.intervals = intervals.iterator();
        this.header = Utils.nonNull(header);
        this.reference = Utils.nonNull(reference);
        this.features = features;
        this.engine = Utils.nonNull(engine);
        evidence = engine.makeLocusEvidence(intervals.get(0).getContig(), intervals.get(0).getStart(), reference, features);
        pendingStates = trackPileups ? new ArrayDeque<>() : null;
    }

    @Override
    public boolean hasNext() {
        return (interval != null && nextLocus <= interval.getEnd()) || intervals.hasNext();
    }

    @Override
    public ActivityProfileState next() {
        if (interval == null || nextLocus > interval.getEnd()) {
            interval = intervals.next();
            intervalReference = new ReferenceContext(reference, interval);
            nextLocus = interval.getStart();
        }
        final int locus = nextLocus++;

        evidence.discardLociBefore(locus);
        while (reads.hasNext() && reads.peek().getStart() <= locus) {
            final GATKRead read = reads.next();
            if (pendingStates != null) {
                PileupBasedAlleles.addMismatchPercentageToRead(read, header, new ReferenceContext(reference, new SimpleInterval(read)));
            }
            evidence.add(read);
        }

        final ActivityProfileState state = engine.evaluateLocus(evidence, new SimpleInterval(interval.getContig(), locus, locus),
                features, intervalReference);
        if (pendingStates != null) {
            pendingStates.add(state);
        }
        return state;
    }

    /**
     * For pileup detection, attach the pileups of the region's loci to it, as {@link AssemblyRegionIterator} does, from
     * the region's reads, which include every read covering them
     */
    @Override
    public void addPileupData(final AssemblyRegion region) {
        if (pendingStates == null) {
            return;
        }
        while (!pendingStates.isEmpty() && pendingStates.peek().getLoc().getStart() < region.getStart()) {
            pendingStates.poll();
        }

        final LocusIteratorByState libs = new LocusIteratorByState(region.getReads().iterator(), DownsamplingMethod.NONE,
                ReadUtils.getSamplesFromHeader(header), header, true);
        final IntervalAlignmentContextIterator pileups = new IntervalAlignmentContextIterator(libs,
                new IntervalLocusIterator(List.of(region.getSpan()).iterator()), header.getSequenceDictionary());
        final List<AlignmentAndReferenceContext> alignmentData = new ArrayList<>(region.getSpan().size());
        for (final ActivityProfileState state : pendingStates) {
            if (state.getLoc().getStart() > region.getEnd()) {
                break;
            }
            final AlignmentContext pileup = pileups.next();
            Utils.validate(pileup.getStart() == state.getLoc().getStart(), "the pileups and activity states of a region are out of step");
            final AlignmentAndReferenceContext data = new AlignmentAndReferenceContext(pileup, new ReferenceContext(reference, new SimpleInterval(pileup)));
            data.setActivityScore(state.getOriginalActiveProb());
            alignmentData.add(data);
        }
        region.addAllAlignmentData(alignmentData);
    }
}
