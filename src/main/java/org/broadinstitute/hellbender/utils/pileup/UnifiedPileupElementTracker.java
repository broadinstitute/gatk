package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class UnifiedPileupElementTracker extends PileupElementTracker {

    /** If {@code true}, the tracker is already sorted by {@link #READ_START_COMPARATOR}. */
    protected boolean sorted;

    private List<PileupElement> elements;

    UnifiedPileupElementTracker() {
        this(Collections.emptyList(), true);
    }

    UnifiedPileupElementTracker(final List<PileupElement> pileup) {
        this(pileup, false);
    }

    UnifiedPileupElementTracker(final List<PileupElement> pileup, final boolean preSorted) {
        this.elements = pileup;
        this.sorted = false;
    }

    @Override
    public Stream<PileupElement> getElementStream() {
        return elements.stream();
    }

    @Override
    public Iterator<PileupElement> sortedIterator() {
        sortTracker();
        return iterator();
    }

    protected void sortTracker() {
        // only sort if it is not sorted yet
        if (!sorted) {
            elements = getElementStream()
                    .sorted(READ_START_COMPARATOR)
                    .collect(Collectors.toList());
            sorted = true;
        }
    }

    @Override
    public PileupElementTracker splitBySample(final SAMFileHeader header) {
        final Set<String> samples = getSamples(header);
        switch (samples.size()) {
            case 0: // no samples means empty pileup
                return this;
            case 1: // only one sample
                return new SingleSamplePileupElementTracker(samples.iterator().next(), elements, sorted);
            default: // more than one sample
                return new PerSamplePileupElementTracker(
                        samples.stream().collect(Collectors.toMap(s -> s, s -> getTrackerForSample(s, header))));
        }
    }

    @Override
    public PileupElementTracker makeFilteredTracker(final Predicate<PileupElement> filter) {
        return new UnifiedPileupElementTracker(getElementStream().filter(filter).collect(Collectors.toList()), sorted);
    }

    @Override
    public PileupElementTracker getTrackerForSample(final String sample,
            final SAMFileHeader header) {
        return new SingleSamplePileupElementTracker(sample,
                getElementStream().filter(pe -> Objects.equals(ReadUtils.getSampleName(pe.getRead(), header), sample))
                .collect(Collectors.toList()));
    }

    @Override
    public Set<String> getSamples(final SAMFileHeader header) {
        return getElementStream().map(PileupElement::getRead)
                .map(r -> ReadUtils.getSampleName(r, header))
                .collect(Collectors.toSet());
    }
}
