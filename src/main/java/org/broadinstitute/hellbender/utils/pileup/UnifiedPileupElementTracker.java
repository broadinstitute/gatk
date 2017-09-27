package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
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
 * Pileup tracker where the elements are unified, that is, not already split by sample or belonging to an unique one.
 *
 * <p>This tracker is used when the origin of the pileup elements provided are a mixture coming from different samples.
 * It has the following advantages:
 *
 * <ul>
 *     <li>It caches the sorted list of elements after calling {@link #sortedIterator()}</li>
 *     <li>{@link #makeFilteredTracker(Predicate)} keeps the sorted status for do not re-sort the list.</li>
 *     <li>{@link #getTrackerForSample(String, SAMFileHeader)} returns a {@link SingleSamplePileupElementTracker} (more efficient implementation of some methods).</li>
 *     <li>{@link #splitBySample(SAMFileHeader)} returns a {@link SingleSamplePileupElementTracker} if a single sample is present (more efficient implementation of some methods)</li>
 *     <li>{@link #splitBySample(SAMFileHeader)} returns a {@link PerSamplePileupElementTracker} if a multiple samples are present (more efficient implementation of some methods)</li>
 * </ul>
 *
 * <p>Note: this classes is fairly low-level, developers should probably confirm that their changes do not belong in
 * a higher-level class such as {@link org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState}
 * or {@link ReadPileup}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class UnifiedPileupElementTracker extends PileupElementTracker {

    /** If {@code true}, the tracker is already sorted by {@link #READ_START_COMPARATOR}. */
    protected boolean sorted;

    private List<PileupElement> elements;

    /** Instantiates an empty element tracker. */
    UnifiedPileupElementTracker() {
        this(Collections.emptyList(), true);
    }

    /** Instantiates an unsorted element tracker.
     *
     * <p>Note: if {@link #sortedIterator()} is requested, sorting will be performed and cached into the tracker.
     *
     * @param pileup list of pileup elements.
     */
    UnifiedPileupElementTracker(final List<PileupElement> pileup) {
        this(pileup, false);
    }

    /**
     * Instantiates a sorted/unsorted element tracker.
     *
     * <p>Note: if {@code preSorted=true} and a {@link #sortedIterator()} is requested, sorting will be performed
     * and cached into the tracker.
     *
     * @param pileup list of pileup elements.
     * @param preSorted {@code true} if the elements are already sorted; {@code false} otherwise.
     */
    UnifiedPileupElementTracker(final List<PileupElement> pileup, final boolean preSorted) {
        this.elements = Utils.nonNull(pileup);
        this.sorted = preSorted;
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

    /** Sorts the tracker. */
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
