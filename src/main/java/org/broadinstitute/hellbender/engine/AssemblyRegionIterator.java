package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.ReadCachingIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

public class AssemblyRegionIterator implements Iterator<AssemblyRegion> {
    private static final Logger logger = LogManager.getLogger(AssemblyRegionIterator.class);

    private final Shard<GATKRead> readShard;
    private final SAMFileHeader readHeader;
    private final ReferenceDataSource reference;
    private final FeatureManager features;
    private final AssemblyRegionEvaluator evaluator;
    private final int minRegionSize;
    private final int maxRegionSize;
    private final int assemblyRegionPadding;
    private final double activeProbThreshold;
    private final int maxProbPropagationDistance;

    private final LocusIteratorByState locusIterator;
    private final ActivityProfile activityProfile;

    private Queue<AssemblyRegion> nextRegions;
    private AssemblyRegion previousRegion;
    private final ReadCachingIterator readCachingIterator;
    private Queue<GATKRead> readCache;
    private final ReadCoordinateComparator readComparator;

    public AssemblyRegionIterator(final Shard<GATKRead> readShard,
                                  final SAMFileHeader readHeader,
                                  final ReferenceDataSource reference,
                                  final FeatureManager features,
                                  final AssemblyRegionEvaluator evaluator,
                                  final int minRegionSize,
                                  final int maxRegionSize,
                                  final int assemblyRegionPadding,
                                  final double activeProbThreshold,
                                  final int maxProbPropagationDistance) {

        Utils.nonNull(readShard);
        Utils.nonNull(readHeader);
        //Utils.nonNull(reference);
        //Utils.nonNull(features);
        Utils.nonNull(evaluator);
        Utils.validateArg(minRegionSize >= 1, "minRegionSize must be >= 1");
        Utils.validateArg(maxRegionSize >= 1, "maxRegionSize must be >= 1");
        Utils.validateArg(minRegionSize <= maxRegionSize, "minRegionSize must be <= maxRegionSize");
        Utils.validateArg(assemblyRegionPadding >= 0, "assemblyRegionPadding must be >= 0");
        Utils.validateArg(activeProbThreshold >= 0.0, "activeProbThreshold must be >= 0.0");
        Utils.validateArg(maxProbPropagationDistance >= 0, "maxProbPropagationDistance must be >= 0");

        this.readShard = readShard;
        this.readHeader = readHeader;
        this.reference = reference;
        this.features = features;
        this.evaluator = evaluator;
        this.minRegionSize = minRegionSize;
        this.maxRegionSize = maxRegionSize;
        this.assemblyRegionPadding = assemblyRegionPadding;
        this.activeProbThreshold = activeProbThreshold;
        this.maxProbPropagationDistance = maxProbPropagationDistance;

        this.readCachingIterator = new ReadCachingIterator(readShard.iterator());
        this.readCache = new ArrayDeque<>();
        this.readComparator = new ReadCoordinateComparator(readHeader);
        this.locusIterator = new LocusIteratorByState(readCachingIterator, DownsamplingMethod.NONE, false, ReadUtils.getSamplesFromHeader(readHeader), readHeader, false);
        this.activityProfile = new BandPassActivityProfile(null, maxProbPropagationDistance, activeProbThreshold, BandPassActivityProfile.MAX_FILTER_SIZE, BandPassActivityProfile.DEFAULT_SIGMA, readHeader);
        this.nextRegions = new ArrayDeque<>();

        loadNextAssemblyRegions();
    }

    @Override
    public boolean hasNext() {
        return nextRegions.size() > 0;
    }

    @Override
    public AssemblyRegion next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when there were no more elements");
        }

        final AssemblyRegion next = nextRegions.poll();
        
        if ( nextRegions.isEmpty() ) {
            loadNextAssemblyRegions();
        }

        return next;
    }

    private void loadNextAssemblyRegions() {
        final List<AssemblyRegion> readyRegions = new ArrayList<>();

        while ( locusIterator.hasNext() && readyRegions.isEmpty() ) {
            final AlignmentContext pileup = locusIterator.next();

            if ( ! activityProfile.isEmpty() ) {
                final boolean forceConversion = pileup.getLocation().getStart() != activityProfile.getEnd() + 1;
                readyRegions.addAll(activityProfile.popReadyAssemblyRegions(assemblyRegionPadding, minRegionSize, maxRegionSize, forceConversion));
            }

            if ( readShard.getPaddedInterval().contains(pileup.getLocation()) ) {
                final SimpleInterval pileupInterval = new SimpleInterval(pileup);
                final ReferenceContext pileupRefContext = new ReferenceContext(reference, pileupInterval);
                final FeatureContext pileupFeatureContext = new FeatureContext(features, pileupInterval);

                final ActivityProfileState profile = evaluator.isActive(pileup, pileupRefContext, pileupFeatureContext);
                activityProfile.add(profile);
            }
        }

        if ( ! locusIterator.hasNext() ) {
            readyRegions.addAll(activityProfile.popReadyAssemblyRegions(assemblyRegionPadding, minRegionSize, maxRegionSize, true));
        }

        if ( ! readyRegions.isEmpty() ) {
            fillAssemblyRegionsWithReads(readyRegions);
        }

        for ( final AssemblyRegion region : readyRegions ) {
            nextRegions.offer(region);
        }
    }

    private void fillAssemblyRegionsWithReads( final List<AssemblyRegion> assemblyRegions ) {
        if ( previousRegion != null ) {
            for ( final GATKRead previousRegionRead : previousRegion.getReads() ) {
                for ( final AssemblyRegion region : assemblyRegions ){
                    if ( region.getExtendedSpan().overlaps(previousRegionRead) ) {
                        region.add(previousRegionRead);
                    }
                }
            }
        }

        readCache.addAll(readCachingIterator.consumeCachedReads());

        final SimpleInterval firstRegionExtendedSpan = assemblyRegions.get(0).getExtendedSpan();
        final SimpleInterval lastRegionExtendedSpan = assemblyRegions.get(assemblyRegions.size() - 1).getExtendedSpan();
        while ( ! readCache.isEmpty() ) {
            final GATKRead nextRead = readCache.peek();

            if ( isBefore(nextRead, firstRegionExtendedSpan) ) {
                // Discard reads that end before the start of the first region's span
                readCache.poll();
                continue;
            }
            if ( isAfter(nextRead, lastRegionExtendedSpan) ) {
                // Stop once we encounter a read that starts after the end of the last region's span
                break;
            }

            readCache.poll();

            for ( final AssemblyRegion region : assemblyRegions ) {
                if ( region.getExtendedSpan().overlaps(nextRead) ) {
                    region.add(nextRead);
                }
            }
        }

        previousRegion = assemblyRegions.get(assemblyRegions.size() - 1);
    }

    private boolean isBefore( final Locatable first, final Locatable second ) {
        final int contigComparison = compareContigs(first, second);
        return contigComparison == -1 || (contigComparison == 0 && first.getEnd() < second.getStart());
    }

    private boolean isAfter( final Locatable first, final Locatable second ) {
        final int contigComparison = compareContigs(first, second);
        return contigComparison == 1 || (contigComparison == 0 && first.getStart() > second.getEnd());
    }

    private int compareContigs( final Locatable first, final Locatable second ) {
        final int firstContigIndex = readHeader.getSequenceIndex(first.getContig());
        final int secondContigIndex = readHeader.getSequenceIndex(second.getContig());

        if (firstContigIndex == secondContigIndex) {
            return 0;
        } else if (firstContigIndex > secondContigIndex) {
            return 1;
        }
        return -1;
    }
    
    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported by AssemblyRegionIterator");
    }


}
