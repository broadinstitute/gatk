package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn.AllelicCountTableVerbosity;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Writes {@link AllelicCount} instances to a tab-separated table file.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCountWriter extends TableWriter<AllelicCount> {

    private static final String PROB_FORMAT = "%.4f";

    private final AllelicCountTableVerbosity verbosity;

    AllelicCountWriter(final File file, final AllelicCountTableVerbosity verbosity, final TableColumnCollection columns) throws IOException {
        super(file, columns);
        this.verbosity = verbosity;
    }

    public AllelicCountWriter(final File file, final AllelicCountTableVerbosity verbosity)
            throws IOException {
        this(file, verbosity, AllelicCountTableColumn.getColumns(verbosity));
    }

    @Override
    protected void composeLine(final AllelicCount record, final DataLine dataLine) {
        composeLine(record, dataLine, verbosity);
    }

    /**
     * Compose the record with given verbosity.
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     * @param verbosity the desired {@link AllelicCountTableVerbosity} of the record
     */
    static void composeLine(final AllelicCount record, final DataLine dataLine, final AllelicCountTableVerbosity verbosity) {
        switch (verbosity) {
            case BASIC:
                composeLineBasic(record, dataLine);
                break;
            case INTERMEDIATE:
                composeLineIntermediate(record, dataLine);
                break;
            case FULL:
                composeLineFull(record, dataLine);
                break;
        }
    }

    /**
     * Compose the basic-verbosity record: (contig, position, ref read count, alt read count)
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLineBasic(final AllelicCount record, final DataLine dataLine) {
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount());
    }

    /**
     * Compose the intermediate-verbosity record: (contig, position, ref read count, alt read count,
     * ref nucleotide, alt nucleotide, read depth)
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLineIntermediate(final AllelicCount record, final DataLine dataLine) {
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount())
                .append(record.getRefNucleotide().name())
                .append(record.getAltNucleotide().name())
                .append(record.getReadDepth());
    }

    /**
     * Compose the full-verbosity record: (contig, position, ref read count, alt read count,
     * ref nucleotide, alt nucleotide, read depth, heterozygosity log odds)
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLineFull(final AllelicCount record, final DataLine dataLine) {
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount())
                .append(record.getRefNucleotide().name())
                .append(record.getAltNucleotide().name())
                .append(record.getReadDepth())
                .append(formatProb(record.getHetLogOdds()));
    }

    static String formatProb(final double value) {
        return String.format(PROB_FORMAT, value);
    }
}
