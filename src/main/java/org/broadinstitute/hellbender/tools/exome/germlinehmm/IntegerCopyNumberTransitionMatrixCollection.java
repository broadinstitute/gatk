package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import com.google.common.collect.Sets;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.ContigGermlinePloidyAnnotationTableColumn;
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
import java.util.*;
import java.util.stream.Collectors;

/**
 * This class holds a collection of {@link IntegerCopyNumberTransitionMatrix} for different sex genotypes
 * and contigs.
 *
 * A new entry (sex genotype, contig) => {@link IntegerCopyNumberTransitionMatrix} can be added to the collection
 * using {@link #add(String, String, IntegerCopyNumberTransitionMatrix)}. The collection can be in an <i>incomplete</i>
 * state
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class IntegerCopyNumberTransitionMatrixCollection {

    /**
     * The map {@link SexGenotypeContigPairKey} => integer copy number transition matrix
     */
    private final Map<SexGenotypeContigPairKey, IntegerCopyNumberTransitionMatrix>
            sexGenotypeAndContigToTransitionMatrixMap;

    private int maxCopyNumberInCollection;

    /**
     * An arbitrary negative value
     */
    private static final int MAX_COPY_NUMBER_IN_EMPTY_COLLECTION = -1;

    /**
     * Public constructor (creates an empty collection)
     */
    public IntegerCopyNumberTransitionMatrixCollection() {
        sexGenotypeAndContigToTransitionMatrixMap = new HashMap<>();
        maxCopyNumberInCollection = MAX_COPY_NUMBER_IN_EMPTY_COLLECTION;
    }

    /**
     * Adds a new entry to the collection.
     *
     * @param sexGenotype sex genotype string identifier
     * @param contig contig string identifier
     * @param transitionMatrix an instance of {@link IntegerCopyNumberTransitionMatrix}
     */
    public void add(@Nonnull final String sexGenotype,
                    @Nonnull final String contig,
                    @Nonnull final IntegerCopyNumberTransitionMatrix transitionMatrix) {
        final SexGenotypeContigPairKey key = SexGenotypeContigPairKey.of(sexGenotype, contig);
        if (sexGenotypeAndContigToTransitionMatrixMap.keySet().contains(key)) {
            throw new IllegalArgumentException(String.format("The integer copy number transition matrix for sex" +
                    " genotype \"%s\" and contig \"%s\" is already defined. Check for duplicate entries in the" +
                    " table file.", sexGenotype, contig));
        }
        sexGenotypeAndContigToTransitionMatrixMap.put(key, transitionMatrix);
        maxCopyNumberInCollection = FastMath.max(maxCopyNumberInCollection, transitionMatrix.getMaxCopyNumber());
    }

    /**
     * Fetches the appropriate {@link IntegerCopyNumberTransitionMatrix} from the collection.
     *
     * @param sexGenotype sex genotype string identifier
     * @param contig contig string identifier
     * @return a reference to the appropriate {@link IntegerCopyNumberTransitionMatrix} in the collection
     */
    public IntegerCopyNumberTransitionMatrix get(@Nonnull final String sexGenotype,
                                                 @Nonnull final String contig) {
        return sexGenotypeAndContigToTransitionMatrixMap.get(checkKey(sexGenotype, contig));
    }

    /**
     * Returns a set of sex genotypes that are annotated at least once
     *
     * @return a set sex genotype identifier strings
     */
    public Set<String> getSexGenotypes() {
        return sexGenotypeAndContigToTransitionMatrixMap.keySet().stream()
                .map(SexGenotypeContigPairKey::getSexGenotype)
                .collect(Collectors.toSet());
    }

    /**
     * Returns a set of contigs that are annotated at least once
     *
     * @return a set of contig identifier strings
     */
    public Set<String> getContigs() {
        return sexGenotypeAndContigToTransitionMatrixMap.keySet().stream()
                .map(SexGenotypeContigPairKey::getContig)
                .collect(Collectors.toSet());
    }

    /**
     * Reads the collection from a tab-separated file. Refer to the documentation of
     * {@link IntegerCopyNumberTransitionMatrixCollection} for an example of the input formatting.
     *
     * @param tableFile a tab-separated collection table
     *
     * @return an instance of {@link IntegerCopyNumberTransitionMatrixCollection}
     */
    public static IntegerCopyNumberTransitionMatrixCollection read(final File tableFile) {
        return IntegerCopyNumberTransitionMatrixCollectionReader.read(tableFile);
    }

    /**
     * Reads the collection from a reader. Refer to the documentation of
     * {@link IntegerCopyNumberTransitionMatrixCollection} for an example of the input formatting.
     *
     * @param reader a reader instance
     * @param parentPath the parent path in relation to which individual matrix files reside
     *
     * @return an instance of {@link IntegerCopyNumberTransitionMatrixCollection}
     */
    public static IntegerCopyNumberTransitionMatrixCollection read(final Reader reader,
                                                                   final String parentPath) {
        return IntegerCopyNumberTransitionMatrixCollectionReader.read(reader, parentPath);
    }

    /**
     * Asserts that all of contained sex genotypes and contigs have assigned copy number transition
     * matrices.
     *
     * @param sexGenotypes a set of sex genotype identifier strings
     * @param contigs a set of contig strings
     * @throws IllegalArgumentException if a required data entry is not available
     */
    private void assertCompleteness(@Nonnull final Set<String> sexGenotypes,
                                    @Nonnull final Set<String> contigs) {
        for (final String sexGenotype : sexGenotypes) {
            for (final String contig : contigs) {
                checkKey(sexGenotype, contig);
            }
        }
    }

    /**
     * Asserts that all pairwise choices of (complete sex genotype, contig) have assigned copy number
     * transition matrices.
     */
    public void assertCompleteness() {
        assertCompleteness(getSexGenotypes(), getContigs());
    }

    /**
     * Returns the maximum copy number state in the collection.
     *
     * @throws IllegalStateException if the collection is empty
     * @return a non-negative integer
     */
    public int getMaxCopyNumberInCollection() {
        if (maxCopyNumberInCollection == MAX_COPY_NUMBER_IN_EMPTY_COLLECTION) {
            throw new IllegalStateException("The collection is empty and maximum copy number is undefined");
        } else {
            return maxCopyNumberInCollection;
        }
    }

    /**
     * Checks where the pair (sexGenotype, contig) points to an entry in the collection, in which case,
     * an instance of {@link SexGenotypeContigPairKey} is generated and returned.
     *
     * @param sexGenotype sex genotype string identifier
     * @param contig contig string identifier
     * @throws IllegalArgumentException if the pair (sexGenotype, contig) has no entry
     * @return an instance of {@link SexGenotypeContigPairKey}
     */
    private SexGenotypeContigPairKey checkKey(@Nonnull final String sexGenotype,
                                              @Nonnull final String contig) {
        final SexGenotypeContigPairKey key = SexGenotypeContigPairKey.of(sexGenotype, contig);
        Utils.validateArg(sexGenotypeAndContigToTransitionMatrixMap.containsKey(key),
                String.format("The collection does not include the integer copy number transition matrix" +
                                " for sex genotype \"%s\" and contig \"%s\"", sexGenotype, contig));
        return key;
    }

    /**
     * Finds the highest copy number for all sex genotypes and contigs, returns a clone with all transition
     * matrices padded to the highest copy number.
     *
     * @return a new instance of {@link IntegerCopyNumberTransitionMatrixCollection}
     */
    public IntegerCopyNumberTransitionMatrixCollection cloneWithPaddedStates() {
        final IntegerCopyNumberTransitionMatrixCollection newCollection =
                new IntegerCopyNumberTransitionMatrixCollection();
        if (!sexGenotypeAndContigToTransitionMatrixMap.isEmpty()) {
            final int highestCopyNumber = sexGenotypeAndContigToTransitionMatrixMap.values().stream()
                    .mapToInt(IntegerCopyNumberTransitionMatrix::getMaxCopyNumber)
                    .max().getAsInt();
            sexGenotypeAndContigToTransitionMatrixMap.keySet()
                    .forEach(key -> {
                        final IntegerCopyNumberTransitionMatrix transitionMatrix = sexGenotypeAndContigToTransitionMatrixMap.get(key);
                        newCollection.add(key.getSexGenotype(), key.getContig(),
                                new IntegerCopyNumberTransitionMatrix(transitionMatrix.getTransitionMatrix(),
                                        highestCopyNumber - transitionMatrix.getMaxCopyNumber()));
                    });
        }
        return newCollection;
    }

    /**
     * This class reads a collection of copy number transition matrices. The collection is a tab-separated table
     * with a mandatory column "CONTIG" and at least an additional column named after a sex genotype string
     * identifier. The rows are different contigs. For example, for two sex genotypes "SEX_XX" and "SEX_XY", the
     * table is formatted as follows:
     *
     *      CONTIG    SEX_XX                      SEX_XY
     *      1         t_matrix_table_XX_1.tsv     t_matrix_table_XY_1.tsv
     *      2         t_matrix_table_XX_2.tsv     t_matrix_table_XY_2.tsv
     *      3         t_matrix_table_XX_3.tsv     t_matrix_table_XY_3.tsv
     *      ...       ...                         ...
     *      X         t_matrix_table_XX_X.tsv     t_matrix_table_XY_X.tsv
     *      Y         t_matrix_table_XX_Y.tsv     t_matrix_table_XY_Y.tsv
     *
     * The T-matrix table files ("t_matrix_table_[...].tsv") are expected to be locatable in the same directory as
     * the collection table file. These files do not need to be unique, e.g. one can use the same copy number prior
     * table for all autosomal contigs.
     *
     * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
     */
    private static final class IntegerCopyNumberTransitionMatrixCollectionReader {
        private static final String CONTIG_TABLE_COLUMN = ContigGermlinePloidyAnnotationTableColumn.CONTIG.name();
        private static final TableColumnCollection MANDATORY_COLUMNS = new TableColumnCollection(CONTIG_TABLE_COLUMN);

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
         * @param parentPath the parent path in relation to which individual matrix files reside
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
                TableUtils.checkMandatoryColumns(tableReader.columns(), MANDATORY_COLUMNS, extraMessage -> new UserException.BadInput(
                        String.format("Error parsing \"%s\": %s", parentPath, extraMessage)));
                final Set<String> sexGenotypes = Sets.difference(new HashSet<>(tableReader.columns().names()),
                        Collections.singleton(CONTIG_TABLE_COLUMN));
                if (sexGenotypes.isEmpty()) {
                    throw new UserException.BadInput("The copy number transition matrix collection must contain annotations" +
                            " for at least one sex genotype");
                }
                tableReader.iterator().forEachRemaining(row ->
                        sexGenotypes.forEach(sexGenotype -> collection.add(sexGenotype, row.getContig(),
                                IntegerCopyNumberTransitionMatrix.read(new File(parentPath,
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
                return dataLine.get(CONTIG_TABLE_COLUMN);
            }

            String getTransitionMatrixTablePathForSexGenotype(final String sexGenotype) {
                return dataLine.get(sexGenotype);
            }
        }
    }
}