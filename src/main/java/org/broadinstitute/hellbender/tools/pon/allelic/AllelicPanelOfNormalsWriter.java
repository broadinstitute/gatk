package org.broadinstitute.hellbender.tools.pon.allelic;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * Writes an {@link AllelicPanelOfNormals} to a tab-separated table file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AllelicPanelOfNormalsWriter extends TableWriter<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> {
    private static final String DOUBLE_FORMAT = "%6.8f";

    AllelicPanelOfNormalsWriter(final File file) throws IOException {
        super(file, AllelicPanelOfNormalsTableColumn.COLUMNS);
    }

    @Override
    protected void composeLine(final Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues> record,
                               final DataLine dataLine) {
        dataLine.append(record.getKey().getContig())
                .append(record.getKey().getEnd())
                .append(formatDouble(record.getValue().getAlpha()))
                .append(formatDouble(record.getValue().getBeta()));
    }

    static String formatDouble(final double value) {
        return String.format(DOUBLE_FORMAT, value);
    }
}
