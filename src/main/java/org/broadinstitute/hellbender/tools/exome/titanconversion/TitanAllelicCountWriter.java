package org.broadinstitute.hellbender.tools.exome.titanconversion;

import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

public class TitanAllelicCountWriter extends TableWriter<AllelicCount> {

    public TitanAllelicCountWriter(final File file) throws IOException {
        super(file, TitanAllelicCountTableColumn.COLUMNS);
    }

    @Override
    protected void composeLine(final AllelicCount record, final DataLine dataLine) {

        // Chr	Position	Ref	RefCount	Nref	NrefCount	NormQuality
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefNucleotide().name())
                .append(record.getRefReadCount())
                .append(record.getAltNucleotide().name())
                .append(record.getAltReadCount())
                .append("");
    }
}
