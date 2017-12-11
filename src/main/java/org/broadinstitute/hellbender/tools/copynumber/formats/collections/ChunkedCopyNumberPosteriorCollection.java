package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.CopyNumberPosteriorRecord;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberStateCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Collection of copy number posteriors for an individual chunk containing a subset of intervals considered in the analysis
 */
public class ChunkedCopyNumberPosteriorCollection extends AbstractSampleRecordCollection<CopyNumberPosteriorRecord> {

    public ChunkedCopyNumberPosteriorCollection(final File inputFile,
                                                final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        super(inputFile,
                Utils.nonNull(integerCopyNumberStateCollection).getTableColumnCollection(),
                getPosteriorRecordFromDataLineDecoder(integerCopyNumberStateCollection),
                getPosteriorRecordToDataLineEncoder(integerCopyNumberStateCollection));
    }

    /**
     * We assume that the posteriors are stored in a log space
     */
    private static Function<DataLine, CopyNumberPosteriorRecord> getPosteriorRecordFromDataLineDecoder(final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        return dataLine -> {
            final Map<IntegerCopyNumberState, Double> copyNumberStateDoubleMap = new HashMap<>();
            for(int i = 0; i < integerCopyNumberStateCollection.size(); i++) {
                copyNumberStateDoubleMap.putIfAbsent(integerCopyNumberStateCollection.get(i), dataLine.getDouble(i));
            }
            try {
                final CopyNumberPosteriorRecord record = new CopyNumberPosteriorRecord(copyNumberStateDoubleMap);
                return record;
            } catch (IllegalArgumentException ex) {
                throw new UserException.BadInput("Validation error occured on line %d of the posterior file: " + String.format(ex.getMessage(), dataLine.getLineNumber()));
            }
        };
    }

    /**
     * Store posteriors in a log space when writing to a file
     */
    private static BiConsumer<CopyNumberPosteriorRecord, DataLine> getPosteriorRecordToDataLineEncoder(final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        return (copyNumberPosteriorRecord, dataLine) -> {
            integerCopyNumberStateCollection.getCopyNumberStates().stream().forEach(state -> dataLine.append(copyNumberPosteriorRecord.getCopyNumberPosterior(state)));
        };
    }

}
