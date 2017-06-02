package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Writes sex genotype data to file and abstract writers.
 *
 * Extended genotyping inference data will be written if the corresponding optional columns
 * are passed to the constructor (in which case, all to-be-written {@link SexGenotypeData}
 * entries must have extended genotyping inference data).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SexGenotypeTableWriter extends TableWriter<SexGenotypeData> {

    private static final String PROB_FORMAT = "%.4e";
    private final List<String> sexGenotypesList;

    /**
     * Public constructor for writing to a writer.
     *
     * @param writer an instance of {@link Writer}
     * @param columns a collection of table columns
     * @throws IOException if a write error occurs
     */
    public SexGenotypeTableWriter(@Nonnull final Writer writer,
                                  @Nonnull final TableColumnCollection columns) throws IOException {
        super(writer, columns);
        if (!new HashSet<>(columns.names()).containsAll(SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS_SET)) {
            throw new UserException.BadInput("Some mandatory columns are not passed; mandatory columns: " +
                SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS_SET.stream().collect(Collectors.joining(", ")));
        }
        final Set<String> sexGenotypesSet = Sets.difference(new HashSet<>(columns.names()),
                new HashSet<>(SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS.names()));
        sexGenotypesList = new ArrayList<>();
        sexGenotypesList.addAll(sexGenotypesSet);
    }

    /**
     * Public constructor for writing to a file.
     *
     * @param file output file
     * @param columns a collection of table columns
     * @throws IOException if a write error occurs
     */
    public SexGenotypeTableWriter(@Nonnull final File file,
                                  @Nonnull final TableColumnCollection columns) throws IOException {
        super(file, columns);
        final Set<String> sexGenotypesSet = Sets.difference(new HashSet<>(columns.names()),
                new HashSet<>(SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS.names()));
        sexGenotypesList = new ArrayList<>();
        sexGenotypesList.addAll(sexGenotypesSet);
    }

    /**
     * Compose a {@link DataLine} from {@link SexGenotypeData}
     *
     * @param record an instance of {@link SexGenotypeData}
     * @param dataLine an instance of {@link DataLine} to populate
     */
    @Override
    protected void composeLine(@Nonnull final SexGenotypeData record, @Nonnull final DataLine dataLine) {
        dataLine.append(record.getSampleName())
                .append(record.getSexGenotype());
        if (sexGenotypesList.size() > 0) {
            if (!record.hasExtendedGenotypingInfo() || !record.getSexGenotypesSet().containsAll(sexGenotypesList)) {
                throw new UserException.BadInput("The sex genotype record for \"" + record.getSampleName() + "\" does not" +
                        " contain the required extended genotyping information");
            }
            final double[] orderedLogLikelihoods = sexGenotypesList.stream()
                    .map(record::getLogLikelihoodPerGenotype).mapToDouble(d -> d).toArray();
            Arrays.stream(orderedLogLikelihoods).forEach(d -> dataLine.append(formatProb(d)));
        }
    }

    /**
     * Format double numbers for output.
     *
     * @param value a double value
     * @return formatted string
     */
    static String formatProb(final double value) {
        return String.format(PROB_FORMAT, value);
    }

}
