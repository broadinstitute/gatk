package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.Sets;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import javax.annotation.Nonnull;
import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A collection of {@link SexGenotypeData} along with I/O methods.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SexGenotypeDataCollection {
    private final List<SexGenotypeData> sexGenotypeDataList;

    /**
     * Public constructor (empty collection).
     */
    public SexGenotypeDataCollection() {
        sexGenotypeDataList = new ArrayList<>();
    }

    /**
     * Public constructor from a file.
     *
     * @param sexGenotypeDataFile sex genotype data file
     * @throws IOException if a read error occurs
     */
    public SexGenotypeDataCollection(@Nonnull final File sexGenotypeDataFile) throws IOException {
        this(sexGenotypeDataFile.getAbsolutePath(), new FileReader(sexGenotypeDataFile));
    }

    /**
     * Public constructor from a reader.
     *
     * @param sexGenotypeDataSourceName a string identifer for the reader
     * @param sexGenotypeDataReader an instance of {@link Reader}
     * @throws IOException if a read error occurs
     */
    public SexGenotypeDataCollection(@Nonnull final String sexGenotypeDataSourceName,
                                     @Nonnull final Reader sexGenotypeDataReader) throws IOException {
        try (final SexGenotypeTableReader reader = new SexGenotypeTableReader(sexGenotypeDataSourceName, sexGenotypeDataReader)) {
            sexGenotypeDataList = reader.stream().collect(Collectors.toList());
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read sex genotype data from " + sexGenotypeDataSourceName);
        }
    }

    /**
     * Add an instance of {@link SexGenotypeData} to the collection
     *
     * @param sexGenotypeData an instance of {@link SexGenotypeData}
     */
    public void add(@Nonnull final SexGenotypeData sexGenotypeData) {
        sexGenotypeDataList.add(sexGenotypeData);
    }

    /**
     * Write the collection to a writer. If extended genotyping inference data is available, they will
     * be also written as optional columns.
     *
     * @param dataWriter an instance of {@link Writer}
     */
    public void write(@Nonnull final Writer dataWriter) {
        if (sexGenotypeDataList.isEmpty()) {
            throw new IllegalStateException("The sex genotype data collection is empty");
        }
        /* check if extended genotyping information is available; first, check for nulls */
        boolean extended = true;
        if (sexGenotypeDataList.stream().filter(dat -> !dat.hasExtendedGenotypingInfo()).count() > 0) {
            extended = false;
        }
        Set<String> commonGenotypes = new HashSet<>();
        /* check if there is a non-empty intersection */
        if (extended) {
            /* pool all genotypes */
            commonGenotypes = new HashSet<>(sexGenotypeDataList.get(0).getSexGenotypesSet());
            for (final SexGenotypeData dat : sexGenotypeDataList) {
                commonGenotypes = Sets.intersection(commonGenotypes, dat.getSexGenotypesSet());
            }
            if (commonGenotypes.isEmpty()) {
                extended = false;
            }
        }

        final TableColumnCollection columns;
        if (extended) {
            final List<String> columnNames = new ArrayList<>();
            columnNames.addAll(SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS.names());
            columnNames.addAll(commonGenotypes);
            columns = new TableColumnCollection(columnNames);
        } else {
            columns = new TableColumnCollection(SexGenotypeTableColumn.MANDATORY_SEX_GENOTYPE_COLUMNS.names());
        }
        try (final SexGenotypeTableWriter writer = new SexGenotypeTableWriter(dataWriter, columns)) {
            writer.writeAllRecords(sexGenotypeDataList);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not write sex genotype data", e);
        }
    }

    /**
     * Return the list of contained {@link SexGenotypeData} instances.
     *
     * @return a list of {@link SexGenotypeData}.
     */
    public List<SexGenotypeData> getSexGenotypeDataList() {
        return sexGenotypeDataList;
    }

    /**
     * Equality comparison is based on the pairwise equality of the contained {@link SexGenotypeData}
     * entries.
     *
     * @param obj object to compare with
     * @return boolean
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof SexGenotypeDataCollection)) {
            return false;
        }

        final SexGenotypeDataCollection col = (SexGenotypeDataCollection) obj;
        return sexGenotypeDataList.equals(col.getSexGenotypeDataList());
    }

    @Override
    public int hashCode() {
        return sexGenotypeDataList.hashCode();
    }
}
