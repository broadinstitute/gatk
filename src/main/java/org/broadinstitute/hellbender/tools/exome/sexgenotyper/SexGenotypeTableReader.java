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
 * Reads sex genotype data from tab-separated files and abstract readers.
 *
 * Extended genotyping inference data will be loaded if the input tab-separated source contains
 * such optional columns.
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
        final TableColumnCollection columns = this.columns();
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
