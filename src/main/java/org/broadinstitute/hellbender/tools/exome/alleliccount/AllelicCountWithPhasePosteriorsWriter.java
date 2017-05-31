package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn.AllelicCountTableVerbosity;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Writes {@link AllelicCountWithPhasePosteriors} instances to a tab-separated table file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountWithPhasePosteriorsWriter extends TableWriter<AllelicCountWithPhasePosteriors> {

    private final AllelicCountTableVerbosity verbosity;

    public AllelicCountWithPhasePosteriorsWriter(final File file, final AllelicCountTableVerbosity verbosity) throws IOException {
        super(file, PhasePosteriorsTableColumn.appendPhasePosteriorColumns(AllelicCountTableColumn.getColumns(verbosity)));
        this.verbosity = verbosity;
    }

    @Override
    protected void composeLine(final AllelicCountWithPhasePosteriors record, final DataLine dataLine) {
        AllelicCountWriter.composeLine(record, dataLine, verbosity);
        composeLinePhasePosteriors(record, dataLine);
    }

    /**
     * Compose the record for the phase posteriors.
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLinePhasePosteriors(final AllelicCountWithPhasePosteriors record, final DataLine dataLine) {
        dataLine.append(AllelicCountWriter.formatProb(record.getRefMinorProb()))
                .append(AllelicCountWriter.formatProb(record.getAltMinorProb()))
                .append(AllelicCountWriter.formatProb(record.getOutlierProb()));
    }
}
