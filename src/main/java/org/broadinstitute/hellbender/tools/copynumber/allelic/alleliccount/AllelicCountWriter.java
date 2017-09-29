package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Writes {@link AllelicCount} instances to a tab-separated file.
 * All {@link AllelicCount} fields must be specified (including ref/alt nucleotide).
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
final class AllelicCountWriter extends TableWriter<AllelicCount> {

    AllelicCountWriter(final File file) throws IOException {
        super(file, AllelicCountTableColumn.COLUMNS);
    }

    @Override
    protected void composeLine(final AllelicCount record, final DataLine dataLine) {
        Utils.validateArg(record.getAltNucleotide() != null && record.getRefNucleotide() != null,
                "AllelicCount must have all fields specified to be written to file.");
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount())
                .append(record.getRefNucleotide().name())
                .append(record.getAltNucleotide().name());
    }
}
