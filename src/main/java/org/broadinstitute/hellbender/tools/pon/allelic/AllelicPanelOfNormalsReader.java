package org.broadinstitute.hellbender.tools.pon.allelic;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.Map;

/**
 * Reads an {@link AllelicPanelOfNormals} from a tab-separated table file.
 *
 * @author Samuel Lee &lt;mehrtash@broadinstitute.org&gt;
 */
final class AllelicPanelOfNormalsReader extends TableReader<Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues>> {
    /**
     * Opens a reader on an a pre-existing allelic panel of normals file.
     *
     * @param file the source file where to read from.
     *
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there is an issue trying to read the contents of the file.
     * @throws RuntimeException if there is a formatting issue within the file.
     */
    AllelicPanelOfNormalsReader(final File file) throws IOException {
        super(file); /* the constructor of TableReader parses the header */
    }

    @Override
    protected Map.Entry<SimpleInterval, AllelicPanelOfNormals.HyperparameterValues> createRecord(final DataLine dataLine) {
        final int position = dataLine.getInt(AllelicPanelOfNormalsTableColumn.POSITION);
        final SimpleInterval interval = new SimpleInterval(dataLine.get(AllelicPanelOfNormalsTableColumn.CONTIG), position, position);
        final double alpha = dataLine.getDouble(AllelicPanelOfNormalsTableColumn.ALPHA);
        final double beta = dataLine.getDouble(AllelicPanelOfNormalsTableColumn.BETA);
        final AllelicPanelOfNormals.HyperparameterValues hyperparameterValues = new AllelicPanelOfNormals.HyperparameterValues(alpha, beta);

        return new AbstractMap.SimpleEntry<>(interval, hyperparameterValues);
    }
}
