package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Writes {@link AllelicCount} instances to a tab-separated table file.
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
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefReadCount())
                .append(record.getAltReadCount())
                .append(record.getRefNucleotide().name())
                .append(record.getAltNucleotide().name());
    }
}
