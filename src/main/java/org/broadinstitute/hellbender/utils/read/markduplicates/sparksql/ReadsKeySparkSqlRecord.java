package org.broadinstitute.hellbender.utils.read.markduplicates.sparksql;

import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

/** Java bean containing the raw fields of ReadsKey for use with Spark Dataset Encoders. */
public class ReadsKeySparkSqlRecord {
    private Long firstReadKeyValue;
    private Long secondReadKeyValue;

    public ReadsKeySparkSqlRecord() {
    }

    public ReadsKeySparkSqlRecord(Long firstReadKeyValue, Long secondReadKeyValue) {
        this.firstReadKeyValue = firstReadKeyValue;
        this.secondReadKeyValue = secondReadKeyValue;
    }

    public static ReadsKeySparkSqlRecord fromReadsKey(ReadsKey readsKey) {
        if (readsKey instanceof ReadsKey.KeyForFragment) {
            return new ReadsKeySparkSqlRecord(((ReadsKey.KeyForFragment) readsKey).getKeyValue(), null);
        } else if (readsKey instanceof ReadsKey.KeyForPair) {
            return new ReadsKeySparkSqlRecord(((ReadsKey.KeyForPair) readsKey).getFirstReadKeyValue(), ((ReadsKey.KeyForPair) readsKey).getSecondReadKeyValue());
        }
        throw new IllegalArgumentException("Unsupported reads key: " + readsKey);
    }

    public Long getFirstReadKeyValue() {
        return firstReadKeyValue;
    }

    public void setFirstReadKeyValue(Long firstReadKeyValue) {
        this.firstReadKeyValue = firstReadKeyValue;
    }

    public Long getSecondReadKeyValue() {
        return secondReadKeyValue;
    }

    public void setSecondReadKeyValue(Long secondReadKeyValue) {
        this.secondReadKeyValue = secondReadKeyValue;
    }
}
