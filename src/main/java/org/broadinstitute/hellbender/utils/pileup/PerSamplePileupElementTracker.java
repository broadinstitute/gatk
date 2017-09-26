package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.MergingIterator;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.AbstractMap;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Pileup tracker where the elements comes from multiple samples
 *
 * <p>This tracker is used when the all the pileup elements are already split by sample, providing more
 * efficient implementation for:
 *
 * <ul>
 *     <li>{@link #sortedIterator()}</li>
 *     <li>{@link #splitBySample(SAMFileHeader)}</li>
 *     <li>{@link #getSamples(SAMFileHeader)}</li>
 *     <li>{@link #fixOverlaps()}</li>
 *     <li>{@link #getTrackerForSample(String, SAMFileHeader)}</li>
 * </ul>
 *
 * <p>Note: this classes is fairly low-level, developers should probably confirm that their changes do not belong in
 * a higher-level class such as {@link org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState}
 * or {@link ReadPileup}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
class PerSamplePileupElementTracker extends PileupElementTracker {

    private final Map<String, PileupElementTracker> stratified;

    /** Instantiates a per-sample element tracker. */
    PerSamplePileupElementTracker(final Map<String, PileupElementTracker> stratified) {
        Utils.nonNull(stratified, "null stratified PileupElementTracker");
        Utils.nonEmpty(stratified.keySet(), "empty stratified PileupElementTracker");
        this.stratified = stratified;
    }

    @Override
    public Stream<PileupElement> getElementStream() {
        return stratified.values().stream().flatMap(PileupElementTracker::getElementStream);
    }

    @Override
    public Iterator<PileupElement> sortedIterator() {
        return stratified.values().stream()
                // it is pre-sorted
                .flatMap(v -> StreamSupport.stream(Spliterators.spliteratorUnknownSize(v.sortedIterator(), Spliterator.SORTED), false))
                .sorted(READ_START_COMPARATOR).iterator();
    }

    @Override
    public Set<String> getSamples(final SAMFileHeader header) {
        return stratified.keySet();
    }

    @Override
    public PileupElementTracker splitBySample(final SAMFileHeader header) {
        return this;
    }

    @Override
    public void fixOverlaps() {
        stratified.values().forEach(PileupElementTracker::fixOverlaps);
    }

    @Override
    public PileupElementTracker makeFilteredTracker(final Predicate<PileupElement> filter) {
        return new PerSamplePileupElementTracker(stratified.entrySet().stream()
                .map(entry -> new AbstractMap.SimpleEntry<>(entry.getKey(), entry.getValue().makeFilteredTracker(filter)))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));
    }

    public PileupElementTracker getTrackerForSample(final String sample, final SAMFileHeader header) {
        return stratified.getOrDefault(sample, new SingleSamplePileupElementTracker(sample));
    }
}
