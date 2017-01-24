package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBin;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin.ReadCountCovariateBinCollection;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by asmirnov on 7/10/17.
 */
public class BinnedReadCountData extends ReadCountData {

    private final SimpleInterval interval;
    private final Map<ReadCountCovariateBin, Integer> binToCountMap;
    private final ReadCountCovariateBinCollection covariateBinCollection;

    public BinnedReadCountData(final SimpleInterval interval, final ReadCountCovariateBinCollection covariateBinCollection) {
        this.interval = Utils.nonNull(interval);
        this.covariateBinCollection = Utils.nonNull(covariateBinCollection);
        binToCountMap = new HashMap<>();
    }

    protected BinnedReadCountData(final SimpleInterval interval,
                                  final ReadCountCovariateBinCollection covariateBinCollection,
                                  final Map<ReadCountCovariateBin, Integer> binToCountMap) {
        this.interval = Utils.nonNull(interval);
        this.covariateBinCollection = Utils.nonNull(covariateBinCollection);
        this.binToCountMap = Utils.nonNull(binToCountMap);
    }

    @Override
    public void updateReadCount(GATKRead read) {
        ReadCountCovariateBin bin = covariateBinCollection.getReadCountCovariateBin(read);
        binToCountMap.put(bin, binToCountMap.getOrDefault(bin, 0) + 1);
    }

    @Override
    public void appendCountsTo(DataLine dataLine) {
        covariateBinCollection.getCovariateBinList().stream().forEach(
                bin -> dataLine.append(binToCountMap.getOrDefault(bin, 0))
        );
    }

    @Override
    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public TableColumnCollection getReadCountDataColumns() {
        return null;
    }
}
