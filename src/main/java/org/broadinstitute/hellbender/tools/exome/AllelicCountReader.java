package org.broadinstitute.hellbender.tools.exome;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.AllelicCountTableColumns.AllelicCountTableVerbosity;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Reads {@link AllelicCount} instances from a tab-separated table file.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountReader extends TableReader<AllelicCount> {

    private final AllelicCountTableVerbosity verbosity;

    /**
     * Opens a reader on an a pre-existing allelic counts tab-separated file.
     *
     * @param file the source file where to read from.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     * @throws UserException.BadInput if not all mandatory columns can be found.
     */
    public AllelicCountReader(final File file) throws IOException {
        super(file); /* the constructor of TableReader parses the header */

        /* detect verbosity level */
        if (columns().containsAll(AllelicCountTableColumns.FULL_COLUMN_NAME_ARRAY)) {
            verbosity = AllelicCountTableVerbosity.FULL;
        } else if (columns().containsAll(AllelicCountTableColumns.INTERMEDIATE_COLUMN_NAME_ARRAY)) {
            verbosity = AllelicCountTableVerbosity.INTERMEDIATE;
        } else if (columns().containsAll(AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY)) {
            verbosity = AllelicCountTableVerbosity.BASIC;
        } else {
            final Set<String> missingColumns = Sets.difference(
                    new HashSet<>(AllelicCountTableColumns.BASIC_COLUMN_NAME_ARRAY), new HashSet<>(columns().names()));
            throw new UserException.BadInput("Bad header in AllelicCount file. Not all mandatory columns are present." +
                    " Missing: " + StringUtils.join(missingColumns, ", "));
        }
    }

    @Override
    protected AllelicCount createRecord(final DataLine dataLine) {
        /* mandatory (basic) fields */
        final int position = dataLine.getInt(AllelicCountTableColumns.POSITION.name());
        final SimpleInterval interval = new SimpleInterval(
                dataLine.get(AllelicCountTableColumns.CONTIG.name()), position, position);
        final int refReadCount = dataLine.getInt(AllelicCountTableColumns.REF_COUNT.name());
        final int altReadCount = dataLine.getInt(AllelicCountTableColumns.ALT_COUNT.name());

        if (verbosity == AllelicCountTableVerbosity.BASIC) {
            return new AllelicCount(interval, refReadCount, altReadCount);
        } else {
            final Nucleotide refNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumns.REF_NUCLEOTIDE.name()).getBytes()[0]);
            final Nucleotide altNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumns.ALT_NUCLEOTIDE.name()).getBytes()[0]);
            final int readDepth = dataLine.getInt(AllelicCountTableColumns.READ_DEPTH.name());
            if (verbosity == AllelicCountTableVerbosity.INTERMEDIATE) {
                return new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide, readDepth);
            } else { /* verbosity == AllelicCountTableVerbosity.FULL */
                final double hetLogOdds = dataLine.getDouble(AllelicCountTableColumns.HET_LOG_ODDS.name());
                return new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide, readDepth,
                        hetLogOdds);
            }
        }
    }

}
