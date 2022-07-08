package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.JobStatistics;

public class StorageAPIAvroReaderAndBigQueryStatistics {
    public final StorageAPIAvroReader storageAPIAvroReader;
    public final JobStatistics.QueryStatistics queryStatistics;

    public StorageAPIAvroReaderAndBigQueryStatistics(final StorageAPIAvroReader storageAPIAvroReader, final JobStatistics.QueryStatistics queryStatistics) {
        this.storageAPIAvroReader = storageAPIAvroReader;
        this.queryStatistics = queryStatistics;
    }


}
