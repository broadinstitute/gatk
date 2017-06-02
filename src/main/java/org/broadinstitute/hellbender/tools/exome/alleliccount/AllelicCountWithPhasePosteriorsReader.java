package org.broadinstitute.hellbender.tools.exome.alleliccount;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Reads {@link AllelicCountWithPhasePosteriors} instances from a tab-separated table file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountWithPhasePosteriorsReader extends TableReader<AllelicCountWithPhasePosteriors> {

    private final AllelicCountTableColumn.AllelicCountTableVerbosity verbosity;

    /**
     * Opens a reader on an a pre-existing allelic counts with phase posteriors tab-separated file.
     *
     * @param file the source file where to read from.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     * @throws UserException.BadInput if not all mandatory columns can be found.
     */
    public AllelicCountWithPhasePosteriorsReader(final File file) throws IOException {
        super(file); /* the constructor of TableReader parses the header */

        /* detect verbosity level */
        if (columns().containsAll(PhasePosteriorsTableColumn.appendPhasePosteriorColumns(AllelicCountTableColumn.FULL_COLUMNS).names())) {
            verbosity = AllelicCountTableColumn.AllelicCountTableVerbosity.FULL;
        } else if (columns().containsAll(PhasePosteriorsTableColumn.appendPhasePosteriorColumns(AllelicCountTableColumn.INTERMEDIATE_COLUMNS).names())) {
            verbosity = AllelicCountTableColumn.AllelicCountTableVerbosity.INTERMEDIATE;
        } else if (columns().containsAll(PhasePosteriorsTableColumn.appendPhasePosteriorColumns(AllelicCountTableColumn.BASIC_COLUMNS).names())) {
            verbosity = AllelicCountTableColumn.AllelicCountTableVerbosity.BASIC;
        } else {
            final Set<String> missingColumns = Sets.difference(
                    new HashSet<>(PhasePosteriorsTableColumn.appendPhasePosteriorColumns(AllelicCountTableColumn.BASIC_COLUMNS).names()), new HashSet<>(columns().names()));
            throw new UserException.BadInput("Bad header in AllelicCountWithPhasePosteriors file. Not all mandatory columns are present." +
                    " Missing: " + StringUtils.join(missingColumns, ", "));
        }
    }

    @Override
    protected AllelicCountWithPhasePosteriors createRecord(final DataLine dataLine) {
        final AllelicCount count = AllelicCountReader.createAllelicCountWithVerbosity(dataLine, verbosity);
        //if probabilities in file are 0.0, use Double.MIN_VALUE instead so we can take log
        final double refMinorProb = Math.max(dataLine.getDouble(PhasePosteriorsTableColumn.REF_MINOR_PROB), Double.MIN_VALUE);
        final double altMinorProb = Math.max(dataLine.getDouble(PhasePosteriorsTableColumn.ALT_MINOR_PROB), Double.MIN_VALUE);
        final double outlierProb = Math.max(dataLine.getDouble(PhasePosteriorsTableColumn.OUTLIER_PROB), Double.MIN_VALUE);
        return new AllelicCountWithPhasePosteriors(count, Math.log(refMinorProb), Math.log(altMinorProb), Math.log(outlierProb));
    }

    public AllelicCountTableColumn.AllelicCountTableVerbosity getVerbosity() {
        return verbosity;
    }
}
