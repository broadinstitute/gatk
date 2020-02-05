package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Given a {@link MultiIntervalShard} of {@link GATKRead}, iterates over each locus within
 * that shard, and calculates the {@link ActivityProfileState} there, using the provided {@link AssemblyRegionEvaluator}
 * to determine if each site is active.
 *
 * Loads the reads from the shard as lazily as possible to minimize memory usage.
 *
 * NOTE: the provided shard must have appropriate read filters set on it for this traversal type (ie., unmapped
 * and malformed reads must be filtered out).
 */
public class ActivityProfileStateIterator implements Iterator<ActivityProfileState> {
    private static final Logger logger = LogManager.getLogger(ActivityProfileStateIterator.class);

    private final ReferenceDataSource reference;
    private final FeatureManager features;
    private final AssemblyRegionEvaluator evaluator;

    private final Iterator<AlignmentContext> locusIterator;

    /**
     * Constructs an AssemblyRegionIterator over a provided read shard
     *
     * @param readShard MultiIntervalShard containing the reads that will go into the assembly regions.
     *                  Must have a MAPPED filter set on it.
     * @param readHeader header for the reads
     * @param reference source of reference bases (may be null)
     * @param features source of arbitrary features (may be null)
     * @param evaluator evaluator used to determine whether a locus is active
     */
    public ActivityProfileStateIterator(final MultiIntervalShard<GATKRead> readShard,
                                        final SAMFileHeader readHeader,
                                        final ReferenceDataSource reference,
                                        final FeatureManager features,
                                        final AssemblyRegionEvaluator evaluator,
                                        final boolean includeReadsWithDeletionsInIsActivePileups) {

        Utils.nonNull(readShard);
        Utils.nonNull(readHeader);
        Utils.nonNull(evaluator);

        this.reference = reference;
        this.features = features;
        this.evaluator = evaluator;

        // We wrap our LocusIteratorByState inside an IntervalAlignmentContextIterator so that we get empty loci
        // for uncovered locations. This is critical for reproducing GATK 3.x behavior!
        LocusIteratorByState libs = new LocusIteratorByState(readShard.iterator(), DownsamplingMethod.NONE, ReadUtils.getSamplesFromHeader(readHeader), readHeader, includeReadsWithDeletionsInIsActivePileups);
        final IntervalLocusIterator intervalLocusIterator = new IntervalLocusIterator(readShard.getIntervals().iterator());
        this.locusIterator = new IntervalAlignmentContextIterator(libs, intervalLocusIterator, readHeader.getSequenceDictionary());
    }

    @Override
    public boolean hasNext() {
        return locusIterator.hasNext();
    }

    @Override
    public ActivityProfileState next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when there were no more elements");
        }

        final AlignmentContext pileup = locusIterator.next();
        final SimpleInterval pileupInterval = new SimpleInterval(pileup);
        final ReferenceContext pileupRefContext = new ReferenceContext(reference, pileupInterval);
        final FeatureContext pileupFeatureContext = new FeatureContext(features, pileupInterval);
        return evaluator.isActive(pileup, pileupRefContext, pileupFeatureContext);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported by AssemblyRegionIterator");
    }
}
