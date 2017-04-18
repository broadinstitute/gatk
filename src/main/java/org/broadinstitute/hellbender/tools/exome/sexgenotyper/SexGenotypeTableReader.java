package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Reads sex genotype data from tab-separated files and abstract readers. For example, a basic sex genotype
 * table for homo sapiens must be formatted as follows:
 *
 *     <pre>
 *         SAMPLE_NAME          SEX_GENOTYPE
 *         arbitrary_name_0     SEX_XX
 *         arbitrary_name_1     SEX_XY
 *         arbitrary_name_2     SEX_XY
 *         arbitrary_name_3     SEX_XX
 *         ...                  ...
 *     </pre>
 *
 * Here, SEX_XX and SEX_XY are arbitrary string identifiers for XX and XY genotypes. The string identifier
 * the user uses in this table must match the identifiers used in the tab-separated contig ploidy
 * annotation table. See {@link ContigGermlinePloidyAnnotation}.
 *
 * Extended genotyping inference data will be loaded if the input tab-separated source contains such optional
 * columns. Example:
 *
 *     <pre>
 *         SAMPLE_NAME          SEX_GENOTYPE      SEX_XX                                SEX_XY
 *         arbitrary_name_0     SEX_XX            arbitrary_sex_xx_log_likelihood_0     arbitrary_sex_xy_log_likelihood_0
 *         arbitrary_name_1     SEX_XY            arbitrary_sex_xx_log_likelihood_1     arbitrary_sex_xy_log_likelihood_1
 *         arbitrary_name_2     SEX_XY            arbitrary_sex_xx_log_likelihood_2     arbitrary_sex_xy_log_likelihood_2
 *         arbitrary_name_3     SEX_XX            arbitrary_sex_xx_log_likelihood_3     arbitrary_sex_xy_log_likelihood_3
 *         ...                  ...               ...                               ...
 *     </pre>
 *
 * For each sex genotyping appearing in the set of all SEX_GENOTYPE entries, the extended table must include
 * a column with the same name that lists the genotyping log likelihood values (SEX_XX and SEX_XY in the above
 * example).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SexGenotypeTableReader extends TableReader<SexGenotypeData> {
    private final List<String> sexGenotypeNames;

    /**
     * Public constructor from a file.
     *
     * @param file input file
     * @throws IOException if a read error occurs
     */
    public SexGenotypeTableReader(@Nonnull final File file) throws IOException {
        this(new FileReader(file));
    }

    /**
     * Public constructor from a reader.
     *
     * @param sourceReader input reader
     * @throws IOException if a read error occurs
     */
    public SexGenotypeTableReader(@Nonnull final Reader sourceReader) throws IOException {
        this(null, sourceReader);
    }

    /**
     * Public constructor from a reader along with its string identifier.
     *
     * @param sourceName string identifier of the input reader
     * @param sourceReader input reader
     * @throws IOException if a read error occurs
     */
    public SexGenotypeTableReader(@Nullable final String sourceName, @Nonnull Reader sourceReader) throws IOException {
        super(sourceName, sourceReader);
        final TableColumnCollection columns = columns();
        TableUtils.checkMandatoryColumns(columns, SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS,
                UserException.BadInput::new);
        sexGenotypeNames = new ArrayList<>();
        sexGenotypeNames.addAll(Sets.difference(new HashSet<>(columns.names()),
                new HashSet<>(SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS.names())));
    }

    /**
     * Create an instance of {@link SexGenotypeData} from a {@link DataLine}. Automatically determines
     * if extended genotyping inference data is available.
     *
     * @param dataLine the input data line
     * @return an instance of {@link SexGenotypeData}
     */
    @Override
    protected SexGenotypeData createRecord(@Nonnull final DataLine dataLine) {
        final String sampleName = dataLine.get(SexGenotypeTableColumn.SAMPLE_NAME.name());
        final String sexGenotype = dataLine.get(SexGenotypeTableColumn.SEX_GENOTYPE.name());
        if (sexGenotypeNames.size() > 0) {
            final double[] logLikelihoods = new double[sexGenotypeNames.size()];
            IntStream.range(0, sexGenotypeNames.size())
                    .forEach(i -> logLikelihoods[i] = dataLine.getDouble(sexGenotypeNames.get(i)));
            return new SexGenotypeData(sampleName, sexGenotype, sexGenotypeNames, logLikelihoods);
        } else {
            return new SexGenotypeData(sampleName, sexGenotype, null, null);
        }
    }
}
