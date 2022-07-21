package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.JobStatistics;
import com.google.cloud.bigquery.TableResult;

public class BigQueryResultAndStatistics {
    public final TableResult result;
    public final JobStatistics.QueryStatistics queryStatistics;

    public BigQueryResultAndStatistics(final TableResult result, final JobStatistics.QueryStatistics queryStatistics) {
        this.result = result;
        this.queryStatistics = queryStatistics;
    }


}
