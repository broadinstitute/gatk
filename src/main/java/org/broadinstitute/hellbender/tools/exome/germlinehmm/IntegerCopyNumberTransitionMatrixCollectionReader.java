package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import javax.annotation.Nonnull;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * This class reads a collection of copy number transition matrices. The collection is a tab-separated table
 * with a mandatory column "CONTIG_NAME" and at least an additional column named after a sex genotype string
 * identifier. The rows are different contigs. For example, for two sex genotypes "SEX_XX" and "SEX_XY", the
 * table is formatted as follows:
 *
 *      CONTIG_NAME SEX_XX  SEX_XY
 *      1   t_matrix_table_XX_1.tsv  t_matrix_table_XY_1.tsv
 *      2   t_matrix_table_XX_2.tsv  t_matrix_table_XY_2.tsv
 *      3   t_matrix_table_XX_3.tsv  t_matrix_table_XY_3.tsv
 *                           . . .
 *      X   t_matrix_table_XX_X.tsv  t_matrix_table_XY_X.tsv
 *      Y   t_matrix_table_XX_Y.tsv  t_matrix_table_XY_Y.tsv
 *
 * The T-matrix table files ("t_matrix_table_[...].tsv") are expected to reside in the same directory as the
 * collection table file. These files do not need to be unique, e.g. one can use the same copy number prior
 * table for all autosomal contigs.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionMatrixCollectionReader {

    private static final String CONTIG_NAME_TABLE_COLUMN = "CONTIG_NAME";
    private static final TableColumnCollection MANDATORY_COLUMNS = new TableColumnCollection(CONTIG_NAME_TABLE_COLUMN);

    /**
     * Reads the copy number transition prior table file and returns an instance of
     * {@link IntegerCopyNumberTransitionMatrixCollection}.
     *
     * @param tableFile input table
     * @return an instance of {@link IntegerCopyNumberTransitionMatrixCollection}
     */
    public static IntegerCopyNumberTransitionMatrixCollection read(@Nonnull final File tableFile) {
        try {
            return read(new FileReader(tableFile), tableFile.getParent());
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the copy number transition matrix collection" +
                    " table file: " + tableFile.getAbsolutePath());
        }
    }

    /**
     * Reads the copy number transition matrix collection from a reader and returns an instance of
     * {@link IntegerCopyNumberTransitionMatrixCollection}.
     *
     * @param reader a reader instance for the collection
     * @param parentPath the parent path for the transition matrices in the collection
     * @return an instance of {@link IntegerCopyNumberTransitionMatrixCollection}
     */
    public static IntegerCopyNumberTransitionMatrixCollection read(@Nonnull final Reader reader,
                                                                   @Nonnull final String parentPath) {
        Utils.nonNull(reader, "The integer copy number transition matrix collection table reader must be non-null");
        Utils.nonNull(parentPath, "The parent path of the integer copy number transition matrix collection table file" +
                " must be non-null");
        final IntegerCopyNumberTransitionMatrixCollection collection = new IntegerCopyNumberTransitionMatrixCollection();
        try (final TableReader<TransitionMatrixCollectionRow> tableReader = new TableReader<TransitionMatrixCollectionRow>(reader) {
            @Override
            protected TransitionMatrixCollectionRow createRecord(DataLine dataLine) {
                return new TransitionMatrixCollectionRow(dataLine);
            }
        }) {
            TableUtils.checkMandatoryColumns(tableReader.columns(), MANDATORY_COLUMNS, UserException.BadInput::new);
            final Set<String> sexGenotypes = Sets.difference(new HashSet<>(tableReader.columns().names()),
                    Collections.singleton(CONTIG_NAME_TABLE_COLUMN));
            if (sexGenotypes.isEmpty()) {
                throw new UserException.BadInput("The copy number transition matrix collection must contain annotations" +
                        " for at least one sex genotype");
            }
            tableReader.iterator().forEachRemaining(row ->
                sexGenotypes.forEach(sexGenotype -> collection.add(sexGenotype, row.getContig(),
                        IntegerCopyNumberTransitionMatrixData.read(new File(parentPath,
                                row.getTransitionMatrixTablePathForSexGenotype(sexGenotype)), 0)))
            );
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the copy number transition matrix collection table", ex);
        }
        return collection;
    }

    /**
     * A wrapper around {@link DataLine} that represents a row of the copy number transition matrix collection
     */
    private static final class TransitionMatrixCollectionRow {
        private final DataLine dataLine;

        TransitionMatrixCollectionRow(final DataLine dataLine) {
            this.dataLine = Utils.nonNull(dataLine, "Data line must be non-null");
        }

        String getContig() {
            return dataLine.get(CONTIG_NAME_TABLE_COLUMN);
        }

        String getTransitionMatrixTablePathForSexGenotype(final String sexGenotype) {
            return dataLine.get(sexGenotype);
        }
    }

}
