package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

/**
 * Reads {@link AllelicCount} instances from a tab-separated table file.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AllelicCountReader extends TableReader<AllelicCount> {

    /**
     * Opens a reader on an a pre-existing allelic counts tab-separated file.
     */
    AllelicCountReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected AllelicCount createRecord(final DataLine dataLine) {
        final String contig = dataLine.get(AllelicCountTableColumn.CONTIG);
        final int position = dataLine.getInt(AllelicCountTableColumn.POSITION);
        final int refReadCount = dataLine.getInt(AllelicCountTableColumn.REF_COUNT);
        final int altReadCount = dataLine.getInt(AllelicCountTableColumn.ALT_COUNT);
        final Nucleotide refNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumn.REF_NUCLEOTIDE.name()).getBytes()[0]);
        final Nucleotide altNucleotide = Nucleotide.valueOf(dataLine.get(AllelicCountTableColumn.ALT_NUCLEOTIDE.name()).getBytes()[0]);
        final SimpleInterval interval = new SimpleInterval(contig, position, position);
        return new AllelicCount(interval, refReadCount, altReadCount, refNucleotide, altNucleotide);
    }
}
