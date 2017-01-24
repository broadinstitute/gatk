package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.*;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinningConfiguration;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * Reader for a read counts table.
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class BinnedReadCountReader extends TableReader<ReadCountData> {

    private final Function<DataLine, ReadCountData> recordExtractor;

    private ReadCountType readCountType = null;

    private String sampleName = null;

    private ReadCountCovariateBinCollection covariateBinCollection = null;

    public BinnedReadCountReader(final String sourceName, final Reader sourceReader) throws IOException {
        super(sourceName, sourceReader);
        validateArguments();
        final Function<DataLine, SimpleInterval> intervalExtractor = intervalExtractor(columns(),
                (message) -> formatException(message));
        recordExtractor = composeRecordExtractor(intervalExtractor);
    }

    public BinnedReadCountReader(final File file) throws IOException {
        this(Utils.nonNull(file).getPath(), new FileReader(file));
    }

    public BinnedReadCountReader(final Reader sourceReader) throws IOException {
        this(null, sourceReader);
    }

    /**
     * Get the sample name
     *
     * @return sample name
     */
    public String getSampleName() {
        return sampleName;
    }

    /**
     * Constructs a per line interval extractor given the header column names.
     *
     * @param columns               the header column names.
     * @param errorExceptionFactory the error handler to be called when there is any problem resolving the interval.
     * @return never {@code null} if there is enough columns to extract the coordinate information, {@code null} otherwise.
     */
    private static Function<DataLine, SimpleInterval> intervalExtractor(
            final TableColumnCollection columns, final Function<String, RuntimeException> errorExceptionFactory) {

        final int contigColumnNumber = columns.indexOf(ReadCountTableColumn.CONTIG.toString());
        final int startColumnNumber = columns.indexOf(ReadCountTableColumn.START.toString());
        final int endColumnNumber = columns.indexOf(ReadCountTableColumn.END.toString());

        if (contigColumnNumber == -1 || startColumnNumber == -1 || endColumnNumber == -1) {
            return null;
        }

        return (v) -> {
            final String contig = v.get(contigColumnNumber);
            final int start = v.getInt(startColumnNumber);
            final int end = v.getInt(endColumnNumber);
            if (start <= 0) {
                throw errorExceptionFactory.apply(String.format("start position must be greater than 0: %d", start));
            } else if (start > end) {
                throw errorExceptionFactory.apply(String.format("end position '%d' must equal or greater than the start position '%d'", end, start));
            } else {
                return new SimpleInterval(contig, start, end);
            }
        };
    }

    private Function<DataLine, ReadCountData> composeRecordExtractor(final Function<DataLine, SimpleInterval> intervalExtractor) {
        TableColumnCollection columns = this.columns();
        final int[] countColumnIndexes = IntStream.range(0, columns.columnCount())
                .filter(i -> !ReadCountTableColumn.isStandardTargetColumnName(columns.nameAt(i))).toArray();
        return (v) -> {
            final SimpleInterval interval = intervalExtractor.apply(v);
            final Map<String, Integer> columnValues = new HashMap<>();
            IntStream.range(0, columns.columnCount()).forEach(index -> columnValues.put(columns.nameAt(index), v.getInt(countColumnIndexes[index])));
            return ReadCountDataFactory.getReadCountDataObject(readCountType, interval, columnValues, covariateBinCollection);
        };
    }

    @Override
    protected void processCommentLine(final String commentText, final long lineNumber) {
        Pair<ReadCountFileHeaderKey, String> headerKeyValuePair =  ReadCountFileHeaderKey.getHeaderValueForKey(commentText);
        if (headerKeyValuePair != null) {
            switch (headerKeyValuePair.getLeft()) {
                case SAMPLE_NAME:
                    this.sampleName = headerKeyValuePair.getRight();
                    break;
                case READ_COUNT_TYPE:
                    String key = headerKeyValuePair.getRight();
                    ReadCountType type = ReadCountType.getReadCountTypeByName(key);
                    if (type == null) {
                        throw formatException(String.format("%s is not a recognized read count type", key));
                    } else {
                        this.readCountType = type;
                    }
                    break;
                case BINNING_CONFIGURATION:
                    List<ReadCountCovariateBinningConfiguration> binningConfigurationsList =
                            ReadCountCovariateBinningConfiguration.parseParameters(headerKeyValuePair.getRight());
                    covariateBinCollection = new ReadCountCovariateBinCollection(binningConfigurationsList);
                    break;
                default:
                    throw new GATKException.ShouldNeverReachHereException("Should not be able to reach here");
            }
        }
    }

    @Override
    protected ReadCountData createRecord(final DataLine dataLine) {
        return recordExtractor.apply(dataLine);
    }

    private void validateArguments() {
        Utils.validate(readCountType != null, String.format("Header of the read counts file must contain key '%s'",
                ReadCountFileHeaderKey.READ_COUNT_TYPE.getHeaderKeyName()));
        Utils.validate(sampleName != null, String.format("Header of the read counts file must contain key '%s'",
                ReadCountFileHeaderKey.READ_COUNT_TYPE.getHeaderKeyName()));
        //TODO format ReadCountDataFactory.getColumnsOfReadCountType(readCountType).toString() properly
        Utils.validate(columns().containsAll(ReadCountDataFactory.getColumnsOfReadCountType(readCountType, covariateBinCollection).names()),
                "Read count table must contain the following columns: " + ReadCountDataFactory.getColumnsOfReadCountType(readCountType, covariateBinCollection).toString());
        //TODO format TargetTableColumn.MANDATORY_COLUMNS.toString() properly
        Utils.validate(columns().containsAll(ReadCountTableColumn.COLUMNS.names()),
                "Read count table must contain the following columns: " + ReadCountTableColumn.COLUMNS.toString());
    }
}

