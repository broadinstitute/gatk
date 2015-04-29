package org.broadinstitute.hellbender.tools.exome;

import java.util.*;

/**
 * Class that represent a sample with Id and read-groups.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class Sample {

    /**
     * Holds the sample id.
     */
    private final String id;

    /**
     * The sample read-groups.
     */
    private final Collection<String> readGroups;

    /**
     * Creates a sample with an arbitrary number of read-groups.
     *
     * @param id the sample id.
     * @param readGroups the read-groups that make reference to this sample.
     *
     * @throws IllegalArgumentException if {@code id} is {@code null}, {@code readGroups} is {@code null},
     *    it contains any {@code null} or some repeated read-group ids.
     */
    public Sample(final String id, final Collection<String> readGroups) {

        if (id == null) {
            throw new IllegalArgumentException("the id cannot be null");
        } else if (readGroups == null) {
            throw new IllegalArgumentException("the readGroups collection cannot be null");
        }
        this.id = id;
        final Set<String> readGroupSet = new LinkedHashSet<>(readGroups.size());
        for (final String readGroupId : readGroups) {
            if (readGroupId == null) {
                throw new IllegalArgumentException("the input read-groups collection contains null values");
            } else if (!readGroupSet.add(readGroupId)) {
                throw new IllegalArgumentException(
                        String.format("the same read-groups '%s' is found multiple times in the read-groups collection",
                                readGroupId));

            }
        }
        this.readGroups = Collections.unmodifiableCollection(readGroupSet);
    }

    /**
     * Creates a sample with no read-groups.
     * @param id the sample id.
     * @throws IllegalArgumentException if {@code id} is {@code null}.
     */
    public Sample(final String id) {
        this(id,Collections.emptySet());
    }

    /**
     * Creates a sample with a single read-group.
     * @param id the sample id.
     * @param readGroup the only read-group for the new sample.
     * @throws IllegalArgumentException if any, {@code id} or {@code readGroup}, is {@code null}.
     */
    public Sample(final String id, final String readGroup) {
        this(id,Collections.singleton(readGroup));
    }

    /**
     * Returns the sample id.
     * @return never {@code null}.
     */
    public String getId() {
        return id;
    }

    /**
     * Returns the list of read-groups of the sample.
     * @return never {@code null} but an empty list if there is no read-groups.
     *  This list is immutable.
     */
    public Collection<String> readGroups() {
        return readGroups;
    }
}
