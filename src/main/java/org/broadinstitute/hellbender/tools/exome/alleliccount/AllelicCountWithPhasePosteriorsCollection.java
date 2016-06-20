package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Simple data structure to pass and read/write a List of {@link AllelicCountWithPhasePosteriors} objects.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountWithPhasePosteriorsCollection {
    private final List<AllelicCountWithPhasePosteriors> counts;

    public AllelicCountWithPhasePosteriorsCollection() {
        counts = new ArrayList<>();
    }

    /**
     * Constructor from from file. Checks whether the input file has just the basic columns, or full columns.
     * @param inputFile file to read from
     */
    public AllelicCountWithPhasePosteriorsCollection(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);

        try (final AllelicCountWithPhasePosteriorsReader reader = new AllelicCountWithPhasePosteriorsReader(inputFile)) {
            counts = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    /**
     * Adds a new {@link AllelicCountWithPhasePosteriors} to counts.
     */
    public void add(final AllelicCountWithPhasePosteriors allelicCount) {
        counts.add(Utils.nonNull(allelicCount));
    }

    /** Returns an unmodifiable view of the list of AllelicCountWithPhasePosteriors.   */
    public List<AllelicCountWithPhasePosteriors> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * Writes to specified file at different verbosity levels.
     * @param verbosity  verbosity level
     * @param outputFile  file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile, final AllelicCountTableColumn.AllelicCountTableVerbosity verbosity) {
        try (final AllelicCountWithPhasePosteriorsWriter writer = new AllelicCountWithPhasePosteriorsWriter(outputFile, verbosity)) {
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
        if (!(o instanceof AllelicCountWithPhasePosteriorsCollection)) {
            return false;
        }

        final AllelicCountWithPhasePosteriorsCollection allelicCountCollection = (AllelicCountWithPhasePosteriorsCollection) o;
        return counts.equals(allelicCountCollection.counts);
    }

    @Override
    public int hashCode() {
        return counts.hashCode();
    }
}
