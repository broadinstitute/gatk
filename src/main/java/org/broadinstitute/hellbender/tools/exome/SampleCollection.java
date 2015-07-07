package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Collection of samples and corresponding read-groups.
 *
 * <p>
 *     This collection provides an interface to make some simple queries such as what samples and read-group
 *     it contains, which sample a read-group belongs to, etc.
 * </p>
 * <p>
 *     Each sample has a unique and fixed numerical 0-based index in the collection running from 0 to <code>{@link #sampleCount()}-1</code>.
 *     Samples are assigned their index based on the lexicographical order of their ids.
 * </p>
 *
 * <p>
 *     Each read-group has also a numerical 0-based index from 0 to <code>{@link #readGroupCount()-1}</code>.
 *     The order is based on the enclosing sample's index and then lexicographical between read-groups of the
 *     same sample.
 * </p>
 * <p>
 *     Accordingly a read-groups returned by {@link Sample#readGroups()} on a sample object obtained from this collection
 *     will follow the same lexicographical order.
 * </p>
 * <p>
 *     "<i>Orphan</i>" read-groups, those that do not belong to any sample, are last, as if the "{@code null}"
 *     sample's id is greater than any-other id in lexicographical order.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCollection {

    /**
     * Map from sample-id to its index <code>[0 .. {@link #sampleCount()}).</code>
     */
    private final Map<String,Integer> sampleIndexById;

    /**
     * Map from sample-id to its index <code>[0 .. {@link #sampleCount()}).</code>
     */
    private final Map<String,Integer> readGroupIndexById;

    /**
     * Map from group id to sample index.
     */
    private final Map<String,Integer> sampleIndexByGroupId;

    /**
     * List of samples sorted by index.
     */
    private final List<Sample> samples;

    /**
     * List of sample ids sorted by index.
     */
    private final List<String> sampleIds;

    /**
     * List of read group ids.
     */
    private final List<String> readGroupIds;

    /**
     * Creates a new sample-db from the read file header.
     *
     * @param header the input read source header.
     *
     * @throws IllegalArgumentException if {@code header} is {@code null}.
     *  is reserved for internal use.
     */
    public SampleCollection(final SAMFileHeader header) {
        if (header == null) {
            throw new IllegalArgumentException("the header provided cannot be null.");
        }

        // Build map Sample-id -> List(ReadGroup)
        final Map<String,List<SAMReadGroupRecord>> readGroupsBySampleId = header.getReadGroups().stream()
                .filter(rg -> rg.getSample() != null)
                .collect(Collectors.groupingBy(SAMReadGroupRecord::getSample));

        // Sorted list of sample ids.
        final List<String> sampleIds = readGroupsBySampleId.keySet().stream()
                .sorted()
                .collect(Collectors.toList());

        // Sample object sorted by sample id.
        // Their read groups are sorted lexicographically.
        final List<Sample> samples = sampleIds.stream()
                .map(id -> new Sample(id,
                        readGroupsBySampleId.get(id).stream()
                            .map(SAMReadGroupRecord::getId)
                            .sorted()
                            .collect(Collectors.toList())))
                .collect(Collectors.toList());

        // Sorted orphan (without sample-id) read-groups ids.
        final List<String> orphanReadGroupIds = header.getReadGroups().stream()
                .filter(rg -> rg.getSample() == null)
                .map(SAMReadGroupRecord::getId)
                .sorted()
                .collect(Collectors.toList());

        // Read-group ids sorted by sample-id and then lexicographically.
        // Orphan (sample-less) are last.
        final List<String> readGroupIds =
                Stream.concat(
                        samples.stream().flatMap(s -> s.readGroups().stream()),
                        orphanReadGroupIds.stream())
                .collect(Collectors.toList());

        // Create map of sample ids to their index:
        final Map<String,Integer> sampleIndexById = IntStream.range(0, samples.size()).boxed()
                .collect(Collectors.toMap(sampleIds::get, i -> i));

        // Create map of read-group ids to their indices:
        final Map<String,Integer> readGroupIndexById = IntStream.range(0, readGroupIds.size()).boxed()
                .collect(Collectors.toMap(readGroupIds::get,i -> i));

        // Read-group id -> sample-index.
        final Map<String,Integer> sampleIndexByGroupId = new HashMap<>(readGroupIds.size());
        readGroupIds.forEach(id -> {
            final String sampleId = header.getReadGroup(id).getSample();
            sampleIndexByGroupId.put(id, sampleId == null ? -1 : sampleIndexById.get(sampleId));
        });

        this.sampleIds = Collections.unmodifiableList(sampleIds);
        this.readGroupIds = Collections.unmodifiableList(readGroupIds);
        this.samples = Collections.unmodifiableList(samples);
        this.sampleIndexById = Collections.unmodifiableMap(sampleIndexById);
        this.sampleIndexByGroupId = Collections.unmodifiableMap(sampleIndexByGroupId);
        this.readGroupIndexById = Collections.unmodifiableMap(readGroupIndexById);
    }

    /**
     * Returns number of samples in this collection.
     * @return 0 or greater.
     */
    public int sampleCount() {
        return sampleIds.size();
    }

    /**
     * Returns the number of read groups in this collection.
     * @return 0 or greater.
     */
    public int readGroupCount() {
        return readGroupIds.size();
    }

    /**
     * Returns the index of the sample given a read group id.
     * @param readGroupId query read-group id.
     * @throws IllegalArgumentException if {@code readGroupId} is {@code null}.
     * @return {@code -1} if there is not such a read group.
     */
    public int sampleIndexByGroupId(final String readGroupId) {
        if (readGroupId == null) {
            throw new IllegalArgumentException("the read-group id cannot be null");
        }
        return sampleIndexByGroupId.getOrDefault(readGroupId,-1);
    }

    /**
     * Returns the read group index given the read-group id.
     * @param readGroupId the query read-group id.
     * @throws IllegalArgumentException if {@code readGroupId} is {@code null}.
     * @return {@code -1} if there is not such a read-group.
     */
    public int readGroupIndexById(final String readGroupId) {
        if (readGroupId == null) {
            throw new IllegalArgumentException("the input read-group id cannot be null");
        }
        return readGroupIndexById.getOrDefault(readGroupId,-1);
    }

    /**
     * Returns list of sample ids sorted by index.
     * @return never {@code null} but potentially a zero length list. The returned object is
     * read-only list an cannot not (or should not be) modified the the caller.
     */
    public List<String> sampleIds() {
        return sampleIds;
    }

    /**
     * Returns list of samples sorted by index.
     * @return never {@code null} but potentially a zero length list. The returned object is
     * read-only list an cannot not (or should not be) modified the the caller.
     */
    public List<Sample> samples() {
        return samples;
    }

    /**
     * Returns list of read-group ids sorted by index.
     * @return never {@code null} but potentially a zero length list. The returned object is
     * read-only list an cannot not (or should not be) modified the the caller.
     */
    public List<String> readGroups() {
        return readGroupIds;
    }

    /**
     * Returns a sample given its id.
     * @param id the requested sample's id.
     * @return {@code null} if there is not such a sample.
     */
    protected Sample sample(final String id) {
        final Integer index = sampleIndexById.get(id);
        if (index == null) {
            return null;
        } else {
            return samples.get(index);
        }
    }

    /**
     * Returns the sample index for a read.
     *
     * <p>The input read is assumed to come for a source with a header compatible with
     * the one used to compose this collection; that is, the correspondence between read-group ids and
     * samples in the read record's associated header should match the one represented in the collection.</p>
     *
     * <p>When this is not the case, this collection's read-group to sample correspondence takes preference, thus the
     *    sample-id in the read record's read-group is ignored to this effect.
     * </p>
     *
     * <p>This method returns {@code -1} if there is no sample-id that corresponds to the read's read-group. or
     *    if the read does not have a read-group.
     * </p>
     *
     * @param read the query read.
     * @throws IllegalArgumentException if {@code read} is {@code null} or it makes reference to an unknown
     *    read-group.
     * @return {@code -1} if there the read does not have a corresponding sample in this collection,
     * otherwise the corresponding sample's index; always within <code>[0..{@link #sampleCount()})</code>.
     */
    public int sampleIndexByRead(final GATKRead read) {
        final String groupId = read.getReadGroup();
        if (groupId == null) {
            return -1;
        } else {
            final Integer result = sampleIndexByGroupId.get(groupId);
            if (result != null) {
                return result;
            } else {
                throw new IllegalArgumentException(
                        String.format("the input read has a read-group '%s' unknown to this collection: ", groupId));
            }
        }
    }

    /**
     * Returns the read-group index for a read.
     *
     * <p>This method returns {@code -1} if the read has no read-group.</p>
     *
     * @param read the query read.
     * @throws IllegalArgumentException if {@code read} is {@code null} or if it refers to an unknown
     *  read-group.
     * @return {@code -1} if there read does not have a read-group, otherwise the corresponding read-group
     * index; always in <code>[0..{@link #readGroupCount()})</code>.
     */
    public int readGroupIndexByRead(final GATKRead read) {
        final String groupId = read.getReadGroup();
        if (groupId == null) {
            return -1;
        } else {
            final Integer result = readGroupIndexById.get(groupId);
            if (result != null) {
                return result;
            } else {
                throw new IllegalArgumentException(
                        String.format("the input read has a read-group '%s' unknown to this collection: ", groupId));
            }
        }
    }
}
