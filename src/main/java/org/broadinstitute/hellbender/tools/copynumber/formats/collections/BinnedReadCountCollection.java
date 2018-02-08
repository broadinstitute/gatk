package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountData;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountDataFactory;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountTableColumn;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.BinningSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;

/**
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class BinnedReadCountCollection extends AbstractSampleLocatableCollection<BinningSampleLocatableMetadata, ReadCountData> {

    public BinnedReadCountCollection(final File binnedCountFile) {
        super(binnedCountFile, getReadCountDataFromDataLineDecoder(),
                getReadCountDataToDataLineEncoder());
    }

    public BinnedReadCountCollection(final BinningSampleLocatableMetadata metadata,
                                     final List<ReadCountData> readCountData,
                                     final TableColumnCollection tableColumnCollection) {
        super(metadata, readCountData, tableColumnCollection,
                getReadCountDataFromDataLineDecoder(),
                getReadCountDataToDataLineEncoder());
    }

    private static BiFunction<DataLine, BinningSampleLocatableMetadata, ReadCountData> getReadCountDataFromDataLineDecoder() {
        return (dataLine, binningSampleLocatableMetadata) ->  {
            final TableColumnCollection columns = dataLine.columns();
            final ReadCountCovariateBinCollection covariateBinCollection =
                    new ReadCountCovariateBinCollection(binningSampleLocatableMetadata.getCovariateBinningConfigurations());

            final String contig = dataLine.get(ReadCountTableColumn.CONTIG);
            final int start = dataLine.getInt(ReadCountTableColumn.START);
            final int end = dataLine.getInt(ReadCountTableColumn.END);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);

            final Map<String, Integer> columnValues = new HashMap<>();
            columns.names().stream().forEach(column -> columnValues.put(column, dataLine.getInt(columns.indexOf(column))));
            return ReadCountDataFactory.getReadCountDataObject(
                    binningSampleLocatableMetadata.getReadCountType(), interval, columnValues, covariateBinCollection);
        };
    }

    private static BiConsumer<ReadCountData, DataLine> getReadCountDataToDataLineEncoder() {
        return (readCountData, dataLine) -> {
            dataLine.append(readCountData.getContig());
            dataLine.append(readCountData.getStart());
            dataLine.append(readCountData.getEnd());
            readCountData.appendCountsTo(dataLine);
        };
    }

    @Override
    protected Metadata.Type getMetadataType() {
        return Metadata.Type.BINNING_SAMPLE_LOCATABLE;
    }
}
