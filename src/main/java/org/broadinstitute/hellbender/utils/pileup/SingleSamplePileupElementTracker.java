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
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class SingleSamplePileupElementTracker extends UnifiedPileupElementTracker {

    private final String sampleName;

    /** Instantiates an empty element tracker for the sample. */
    SingleSamplePileupElementTracker(final String sampleName) {
        this(sampleName, Collections.emptyList(), true);
    }

    /** Instantiates an unsorted element tracker for the sample. */
    SingleSamplePileupElementTracker(final String sampleName, final List<PileupElement> sampleElements) {
        this(sampleName, sampleElements, false);
    }

    /** Instantiates a sorted/unsorted element tracker for them sample. The sample may be {@code null}. */
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
