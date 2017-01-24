package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBin;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import javax.annotation.Nullable;
import java.util.*;

/**
 * Utility class to construct {@link ReadCountData} objects
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public final class ReadCountDataFactory {

    // private constructor to prevent instantiation
    private ReadCountDataFactory() {}

    /**
     * Construct an empty instance of {@link ReadCountData} object given its type
     *
     * @param type type of the read count data
     * @param interval simple interval
     * @param covariateBinCollection covariate bin collection, {@code null} if read count data type is not {@link ReadCountType#BINNED}
     * @return empty instance of read count data object
     */
    public static ReadCountData getReadCountDataObject(final ReadCountType type,
                                                       final SimpleInterval interval,
                                                       @Nullable final ReadCountCovariateBinCollection covariateBinCollection) {
        switch (type) {
            case SIMPLE_COUNT:
                return new SimpleReadCountData(interval);
            case BINNED:
                Utils.nonNull(covariateBinCollection);
                return new BinnedReadCountData(interval, covariateBinCollection);
            default:
                throw new UserException.BadInput(String.format(" %s is not a recognized read count type.", type));
        }
    }

    /**
     * Construct a populated instance of {@link ReadCountData} object given its type
     *
     * @param type type of the read count data
     * @param interval simple interval
     * @param columnValues column values to populate read count data object with
     * @param covariateBinCollection covariate bin collection, {@code null} if read count data type is not {@link ReadCountType#BINNED}
     * @return populated instance of read count data object
     */
    public static ReadCountData getReadCountDataObject(final ReadCountType type,
                                                       final SimpleInterval interval,
                                                       final Map<String, Integer> columnValues,
                                                       @Nullable final ReadCountCovariateBinCollection covariateBinCollection) {
        switch (type) {
            case SIMPLE_COUNT:
                Integer rawReadCount = Utils.nonNull(columnValues.get(SimpleReadCountData.SIMPLE_COUNT_COLUMN_NAME),
                        String.format(" %s column value is required to construct %s read count data ",
                                SimpleReadCountData.SIMPLE_COUNT_COLUMN_NAME, ReadCountType.SIMPLE_COUNT.getReadCountTypeName()));
                return new SimpleReadCountData(interval, rawReadCount);
            case BINNED:
                Map<ReadCountCovariateBin, Integer> binToCountMap = new HashMap<>();
                columnValues.keySet().stream().forEach(
                        columnName -> {
                            if (columnValues.get(columnName) != 0) {
                                Integer count = columnValues.get(columnName);
                                ReadCountCovariateBin bin = Utils.nonNull(covariateBinCollection.getCovariateBinByName(columnName),
                                        "Covariate bin column name does not correspond to any bins in the covariate bin collection");
                                binToCountMap.put(bin, count);
                            }
                        }
                );
                return new BinnedReadCountData(interval, covariateBinCollection, binToCountMap);
            default:
                throw new UserException.BadInput(String.format(" %s is not a recognized read count type.", type));
        }
    }

    /**
     * Return {@link TableColumnCollection} corresponding to a particular read count type.
     * Note that columns for {@link ReadCountType#BINNED} read count data type are set at runtime, thus a
     * corresponding {@link ReadCountCovariateBinCollection} must be provided
     *
     * @param type type of the read count data
     * @param covariateBinCollection covariate bin collection, {@code null} if read count data type is not {@link ReadCountType#BINNED}
     * @return table column collection
     */
    public static TableColumnCollection getColumnsOfReadCountType(final ReadCountType type,
                                                                  @Nullable final ReadCountCovariateBinCollection covariateBinCollection) {
        switch (type) {
            case SIMPLE_COUNT:
                return new TableColumnCollection(Arrays.asList(SimpleReadCountData.SIMPLE_COUNT_COLUMN_NAME));
            case BINNED:
                Utils.nonNull(covariateBinCollection);
                return covariateBinCollection.getTableColumns();
            default:
                throw new UserException.BadInput(String.format(" %s is not a recognized read count type.", type));
        }
    }

}
