package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.SAMFileHeader;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Pileup tracker where the elements comes from a single sample.
 *
 * <p>This tracker is used when the all the pileup elements are known to come from a single sample, providing more
 * efficient implementation for:
 *
 * <ul>
 *     <li>{@link #splitBySample(SAMFileHeader)}</li>
 *     <li>{@link #getSamples(SAMFileHeader)}</li>
 *     <li>{@link #getTrackerForSample(String, SAMFileHeader)}</li>
 * </ul>
 *
 * <p>Note: this tracker inherits from {@link UnifiedPileupElementTracker} for re-use the implementation for iterator methods.
 *
 * <p>Note: this classes is fairly low-level, developers should probably confirm that their changes do not belong in
 * a higher-level class such as {@link org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState}
 * or {@link ReadPileup}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class SingleSamplePileupElementTracker extends UnifiedPileupElementTracker {

    private final String sampleName;

    /**
     * Instantiates an empty element tracker for the sample.
     *
     * @param sampleName the name of this sample.
     */
    SingleSamplePileupElementTracker(final String sampleName) {
        this(sampleName, Collections.emptyList(), true);
    }

    /**
     * Instantiates an unsorted element tracker for the sample.
     *
     * <p>Note: if a {@link #sortedIterator()} is requested, sorting will be performed and cached into the tracker.
     *
     * @param sampleName the name of this sample.
     * @param sampleElements list of pileup elements.
     */
    SingleSamplePileupElementTracker(final String sampleName, final List<PileupElement> sampleElements) {
        this(sampleName, sampleElements, false);
    }

    /**
     * Instantiates a sorted/unsorted element tracker for them sample. The sample may be {@code null}.
     *
     * <p>Note: if {@code preSorted=true} and a {@link #sortedIterator()} is requested, sorting will be performed
     * and cached into the tracker.
     *
     * @param sampleName the name of this sample.
     * @param sampleElements list of pileup elements.
     * @param preSorted {@code true} if the elements are already sorted; {@code false} otherwise.
     */
    SingleSamplePileupElementTracker(final String sampleName, final List<PileupElement> sampleElements, final boolean preSorted) {
        super(sampleElements, preSorted);
        this.sampleName = sampleName;
    }

    @Override
    public PileupElementTracker splitBySample(final SAMFileHeader header) {
        return this;
    }

    @Override
    public PileupElementTracker makeFilteredTracker(final Predicate<PileupElement> filter) {
        return new SingleSamplePileupElementTracker(sampleName, getElementStream().filter(filter).collect(Collectors.toList()), sorted);
    }

    @Override
    public Set<String> getSamples(final SAMFileHeader header) {
        return Collections.singleton(sampleName);
    }

    @Override
    public PileupElementTracker getTrackerForSample(final String sample,
            final SAMFileHeader header) {
        return (Objects.equals(sampleName, sample)) ? this : new SingleSamplePileupElementTracker(sample, null);
    }
}
