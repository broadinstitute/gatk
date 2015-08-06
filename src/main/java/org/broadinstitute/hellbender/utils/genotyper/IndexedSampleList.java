package org.broadinstitute.hellbender.utils.genotyper;

import org.broadinstitute.hellbender.utils.collections.IndexedSet;

import java.util.Collection;

/**
 * Simple implementation of a sample-list using an indexed-set.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IndexedSampleList implements SampleList {

    private final IndexedSet<String> samples;

    /**
     * Constructs an empty sample-list.
     */
    public IndexedSampleList() {
        samples = new IndexedSet<>(0);
    }

    /**
     * Constructs a sample-list from a collection of samples.
     *
     * <p>
     *     Repeats in the input collection are ignored (just the first occurrence is kept).
     *     Sample names will be sorted based on the traversal order
     *     of the original collection.
     * </p>
     *
     * @param samples input sample collection.
     *
     * @throws IllegalArgumentException if {@code samples} is {@code null} or it contains {@code nulls}.
     */
    public IndexedSampleList(final Collection<String> samples) {
        //note: no checking here - IndexedSet constructor does it
        this.samples = new IndexedSet<>(samples);
    }

    /**
     * Constructs a sample-list from an array of samples.
     *
     * <p>
     *     Repeats in the input array are ignored (just the first occurrence is kept).
     *     Sample names will be sorted based on the traversal order
     *     of the original array.
     * </p>
     *
     * @param samples input sample array.
     *
     * @throws IllegalArgumentException if {@code samples} is {@code null} or it contains {@code nulls}.
     */
    public IndexedSampleList(final String... samples) {
        //note: no checking here - IndexedSet constructor does it
        this.samples = new IndexedSet<>(samples);
    }

    @Override
    public int numberOfSamples() {
        return samples.size();
    }

    @Override
    public int indexOfSample(final String sample) {
        return samples.indexOf(sample);
    }

    @Override
    public String getSample(final int sampleIndex) {
        return samples.get(sampleIndex);
    }
}