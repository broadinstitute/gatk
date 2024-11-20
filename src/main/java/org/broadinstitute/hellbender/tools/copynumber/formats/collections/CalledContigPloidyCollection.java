package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledContigPloidy;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.GermlineCNVNamingConstants;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.function.BiConsumer;
import java.util.function.Function;

public class CalledContigPloidyCollection extends AbstractSampleRecordCollection<CalledContigPloidy> {

    public CalledContigPloidyCollection(final File inputFile) {
        super(inputFile,
                new TableColumnCollection(GermlineCNVNamingConstants.CONTIG_COLUMN,
                        GermlineCNVNamingConstants.PLOIDY_COLUMN,
                        GermlineCNVNamingConstants.PLOIDY_GQ_COLUMN),
                getContigPloidyRecordFromDataLineDecoder(inputFile),
                getContigPloidyRecordFromDataLineEncoder());
    }

    /**
     * Generates an instance of {@link CalledContigPloidy} from a {@link DataLine} entry read from
     * a ploidy calls file.
     */
    private static Function<DataLine, CalledContigPloidy> getContigPloidyRecordFromDataLineDecoder(final File inputFile) {
        return dataLine -> {
            try {
                return new CalledContigPloidy(dataLine.get(GermlineCNVNamingConstants.CONTIG_COLUMN),
                        dataLine.getInt(GermlineCNVNamingConstants.PLOIDY_COLUMN),
                        dataLine.getDouble(GermlineCNVNamingConstants.PLOIDY_GQ_COLUMN));
            } catch (final IllegalArgumentException ex) {
                throw new UserException.BadInput(
                        String.format("Error parsing baseline copy-number file (%s) at line %d.",
                                inputFile.getAbsolutePath(), dataLine.getLineNumber()));
            }
        };
    }

    /**
     * Generates an instance of {@link DataLine} from {@link CalledContigPloidy} for writing a ploidy calls collection to a file.
     */
    private static BiConsumer<CalledContigPloidy, DataLine> getContigPloidyRecordFromDataLineEncoder() {
        return (ploidyCall, dataLine) -> dataLine.append(ploidyCall.getContig())
                .append(ploidyCall.getPloidy())
                .append(ploidyCall.getPloidyQuality());
    }
}
