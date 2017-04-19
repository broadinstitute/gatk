package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

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
 * Simple data structure to pass and read/write a List of {@link AllelicCount} objects.
 * All {@link AllelicCount} fields (including ref/alt nucleotide) must be specified if reading/writing from/to file.
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
     * Constructor from file. All {@link AllelicCount} fields must be specified (including ref/alt nucleotide).
     * @param inputFile file to read from
     * @throws UserException.BadInput if not all {@link AllelicCount} fields are specified
     * @throws UserException.CouldNotReadInputFile if input file cannot be read
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
     * Adds a new {@link AllelicCount}.
     */
    public void add(final AllelicCount allelicCount) {
        counts.add(Utils.nonNull(allelicCount));
    }

    /**
     * Returns an unmodifiable view of the list of {@link AllelicCount}s.
     */
    public List<AllelicCount> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * Writes counts to specified file.  All {@link AllelicCount} fields must be specified (including ref/alt nucleotide).
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     * @throws IllegalArgumentException if not all {@link AllelicCount} fields are specified
     * @throws UserException.CouldNotCreateOutputFile if output file cannot be created
     */
    public void write(final File outputFile) {
        try (final AllelicCountWriter writer = new AllelicCountWriter(outputFile)) {
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
