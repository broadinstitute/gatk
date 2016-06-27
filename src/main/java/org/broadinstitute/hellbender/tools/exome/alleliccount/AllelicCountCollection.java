package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple data structure to pass and read/write a List of {@link AllelicCount} objects.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountCollection {
    private final List<AllelicCount> counts;

    public AllelicCountCollection() {
        counts = new ArrayList<>();
    }

    /**
     * Constructor from from file. Checks whether the input file has just the basic columns, or full columns.
     * @param inputFile file to read from
     */
    public AllelicCountCollection(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);

        try (final AllelicCountReader reader = new AllelicCountReader(inputFile)) {
            counts = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    /**
     * Adds a new {@link AllelicCount} to counts.
     */
    public void add(final AllelicCount allelicCount) {
        counts.add(Utils.nonNull(allelicCount));
    }

    /** Returns an unmodifiable view of the list of AllelicCounts.   */
    public List<AllelicCount> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * @return a map from SimpleIntervals in counts to their integer index
     */
    public Map<SimpleInterval, Integer> getSimpleIntervalToIndexMap() {
        return IntStream.range(0, counts.size()).boxed()
                .collect(Collectors.toMap(i -> counts.get(i).getInterval(), i -> i));
    }

    /**
     * Writes out pulldown data to specified file at different verbosity levels.
     * @param verbosity  verbosity level
     * @param outputFile  file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile, final AllelicCountTableColumn.AllelicCountTableVerbosity verbosity) {
        try (final AllelicCountWriter writer = new AllelicCountWriter(outputFile, verbosity)) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCountCollection)) {
            return false;
        }

        final AllelicCountCollection allelicCountCollection = (AllelicCountCollection) o;
        return counts.equals(allelicCountCollection.counts);
    }

    @Override
    public int hashCode() {
        return counts.hashCode();
    }
}
