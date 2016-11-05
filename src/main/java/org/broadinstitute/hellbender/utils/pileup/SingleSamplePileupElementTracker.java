package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.SAMFileHeader;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class SingleSamplePileupElementTracker extends UnifiedPileupElementTracker {

    private final String sampleName;

    SingleSamplePileupElementTracker(final String sampleName) {
        this(sampleName, Collections.emptyList(), true);
    }

    SingleSamplePileupElementTracker(final String sampleName, final List<PileupElement> sampleElements) {
        this(sampleName, sampleElements, false);
    }

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
