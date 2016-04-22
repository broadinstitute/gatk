package org.broadinstitute.hellbender.utils.genotyper;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * An immutable, indexed set of samples.
 *
 * <p>
 *     Implementing classes must guarantee that the sample list will remain <b>constant</b> through the life of the object.
 * </p>
 */
//Note: Names in this interface are unusual because of name clash in a subclass.
// For example the name of SampleList.numberOfSamples() cannot be simply size(), as would be usual,
// because {@link ReadLikelihoods} implements AlleleList and SampleList and then size() would be ambiguous.
public interface SampleList  {

    static final SampleList EMPTY_LIST = new SampleList() {
        @Override
        public int numberOfSamples() {
            return 0;
        }

        @Override
        public int indexOfSample(final String sample) {
            Utils.nonNull(sample);
            return -1;
        }

        @Override
        public String getSample(final int sampleIndex) {
            throw new IllegalArgumentException("index is out of valid range");  //we know this without checking because it's an empty list
        }
    };

    /**
     * Empty list.
     *
     * @return never {@code null}
     */
    public static SampleList emptySampleList() {
        return EMPTY_LIST;
    }

    /**
     * Returns number of elements in the list.
     */
    public int numberOfSamples();

    /**
     * Returns the index of an object.
     * @param sample the sample of interest.
     *
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     *
     * @return {@code -1} if such a sample is not an element of this set, otherwise is index in the set thus a
     * values within [0,{@link #numberOfSamples()}).
     */
    public int indexOfSample(final String sample);

    /**
     * Returns the element given its index within the set.
     * @param sampleIndex the target samples's index.
     *
     * @throws IllegalArgumentException if {@code index} is not valid; in [0,{@link #numberOfSamples()}).
     *
     * @return never {@code null}; as null is not a valid element.
     */
    public String getSample(final int sampleIndex);

    /**
     * Checks whether two sample lists are in fact the same.
     * @param first one list to compare.
     * @param second another list to compare.
     *
     * @throws IllegalArgumentException if if either list is {@code null}.
     *
     * @return {@code true} iff both list are equal.
     */
    public static boolean equals(final SampleList first, final SampleList second) {
        Utils.nonNull(first, "first list is null");
        Utils.nonNull(second, "second list is null");
        final int sampleCount = first.numberOfSamples();
        if (sampleCount != second.numberOfSamples()) {
            return false;
        }

        for (int i = 0; i < sampleCount; i++) {
            final String firstSample = first.getSample(i);
            final String secondSample = second.getSample(i);
            if (!firstSample.equals(secondSample)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns a {@link List} unmodifiable view of a sample-list
     *
     * @throws IllegalArgumentException if {@code list} is {@code null}.
     *
     * @return Unmodifiable view of the sample list. Never {@code null}.
     */
    default public List<String> asListOfSamples() {
        return new AbstractList<String>() {
                @Override
                public String get(final int index) {
                    return SampleList.this.getSample(index);
                }

                @Override
                public int size() {
                    return SampleList.this.numberOfSamples();
                }
        };
    }

    /**
     * Returns a {@link Set} unmodifiable view of the sample-list
     *
     * @return Unmodifiable view of the sample set. Never null.
     */
    default public Set<String> asSetOfSamples() {
        return new AbstractSet<String>() {
            @Override
            public Iterator<String> iterator() {
                return new Iterator<String>() {
                    private int index = 0;

                    @Override
                    public boolean hasNext() {
                        return index < SampleList.this.numberOfSamples();
                    }

                    @Override
                    public String next() {
                        if (index >= SampleList.this.numberOfSamples()) {
                            throw new NoSuchElementException("iterating beyond sample list end");
                        }
                        return SampleList.this.getSample(index++);
                    }

                    @Override
                    public void remove() {
                        throw new UnsupportedOperationException("unsupported operation exception");
                    }
                };
            }

            @Override
            public int size() {
                return SampleList.this.numberOfSamples();
            }

            @Override
            public boolean contains(final Object obj) {
                //note: this handles null too (instanceof returns false)
                return (obj instanceof String) && SampleList.this.indexOfSample(((String) obj)) >= 0;
            }
        };
    }

    /**
     * Creates a list with a single sample.
     *
     * @param sampleName the sample name.
     * @return never {@code sampleName}
     */
    public static SampleList singletonSampleList(final String sampleName) {
        Utils.nonNull(sampleName, "the sample name cannot be null");
        return new SampleList() {
            @Override
            public int numberOfSamples() {
                return 1;
            }

            @Override
            public int indexOfSample(final String sample) {
                return sampleName.equals(sample) ? 0 : -1;
            }

            @Override
            public String getSample(final int sampleIndex) {
                if (sampleIndex == 0) {
                    return sampleName;
                }
                throw new IllegalArgumentException("index is out of bounds");
            }
        };
    }
}
