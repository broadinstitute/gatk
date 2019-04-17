package org.broadinstitute.hellbender.utils.Bigquery;

import com.google.cloud.bigquery.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.UUID;

/**
 * Utility class for dealing with BigQuery connections / tables / queries /etc.
 *
 * Created by jonn on 4/17/19.
 */
public class BigQueryUtils {
    private static final Logger logger = LogManager.getLogger(BigQueryUtils.class);

    //==================================================================================================================
    // Static Methods:

    /**
     * @return A {@link BigQuery} object that can be used to interact with a BigQuery data set.
     */
    public static BigQuery getBigQueryEndPoint() {
        return BigQueryOptions.getDefaultInstance().getService();
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final String queryString) {
        return executeQuery( getBigQueryEndPoint(), queryString );
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param bigQuery The {@link BigQuery} instance against which to execute the given {@code queryString}.
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final BigQuery bigQuery, final String queryString) {

        // Create a query configuration we can run based on our query string:
        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder( queryString )
                        // Use standard SQL syntax for queries.
                        // See: https://cloud.google.com/bigquery/sql-reference/
                        .setUseLegacySql(false)
                        .build();

        return submitQueryAndWaitForResults( bigQuery, queryConfig );
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param bigQuery The {@link BigQuery} instance against which to execute the given {@code queryString}.  Must contain the table name in the `FROM` clause for the table from which to retrieve data.
     * @param projectID The BigQuery {@code project ID} containing the {@code dataSet} and table from which to query data.
     * @param dataSet The BigQuery {@code dataSet} containing the table from which to query data.
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table ID in the `FROM` clause for the table from which to retrieve data.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final BigQuery bigQuery,
                                           final String projectID,
                                           final String dataSet,
                                           final String queryString) {

        // Create a query configuration we can run based on our query string:
        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder( queryString )
                        .setUseLegacySql(false)
                        .setDefaultDataset(DatasetId.of(projectID, dataSet))
                        .build();

        return submitQueryAndWaitForResults( bigQuery, queryConfig );
    }

    //==================================================================================================================
    // Helper Methods:

    /**
     * Executes the given {@code queryJobConfiguration} on the given {@code bigQuery} instance.
     * @param bigQuery The instance of {@link BigQuery} to use to connect to BigQuery.
     * @param queryJobConfiguration The {@link QueryJobConfiguration} object containing all required information to retrieve data from a BigQuery table.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    private static TableResult submitQueryAndWaitForResults( final BigQuery bigQuery,
                                                             final QueryJobConfiguration queryJobConfiguration ) {
        // Create a job ID so that we can safely retry:
        final JobId jobId = JobId.of(UUID.randomUUID().toString());

        logger.info("Sending query to server...");
        Job   queryJob       = bigQuery.create(JobInfo.newBuilder(queryJobConfiguration).setJobId(jobId).build());

        // Wait for the query to complete.
        try {
            logger.info("Waiting for query to complete...");
            queryJob = queryJob.waitFor();
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        // Check for errors:
        if (queryJob == null) {
            throw new GATKException("Query job no longer exists");
        } else if (queryJob.getStatus().getError() != null) {

            // Get all the errors we found and log them:
            for ( final BigQueryError bigQueryError : queryJob.getStatus().getExecutionErrors() ) {
                logger.error( "Encountered BigQuery Execution Error: " + bigQueryError.toString() );
            }

            // Since we found an error, we should stop and alert the user:
            throw new GATKException(queryJob.getStatus().getError().toString());
        }

        // Get the results.
        logger.info("Retrieving query results...");
        final QueryResponse response = bigQuery.getQueryResults(jobId);
        final TableResult result;
        try {
            result = queryJob.getQueryResults();
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        return result;
    }
}
