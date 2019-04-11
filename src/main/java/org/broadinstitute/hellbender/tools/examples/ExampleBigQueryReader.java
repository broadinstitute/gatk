package org.broadinstitute.hellbender.tools.examples;

// NOTE:
// Adapted from: https://github.com/googlearchive/bigquery-samples-java/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java

import com.google.cloud.bigquery.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.UUID;

/**
 * An example class that communicates with BigQuery using the google bigquery library.
 * Created by jonn on 4/9/19.
 */
@CommandLineProgramProperties(
        summary = "Example tool that communicates with BigQuery using the google bigquery library.  Will print the first several rows from the given BigQuery table.",
        oneLineSummary = "Example tool that prints the first few rows of a given BigQuery table.",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = false
)
public class ExampleBigQueryReader extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExampleBigQueryReader.class);

    //==================================================================================================================
    // Public Static Members:

    /**
     * The name of the argument for the project ID of the BigQuery instance.
     */
    private static final String PROJECT_ID_ARG_LONG_NAME = "project-id";

    /**
     * The name of the argument for the fully qualified input table ID of the BigQuery instance.
     */
    private static final String FQ_TABLE_ID_ARG_LONG_NAME = "fq-table-id";

    /**
     * The name of the argument for the bucket of the BigQuery instance.
     */
    private static final String BUCKET_ARG_LONG_NAME = "bucket";

    /**
     * The name of the argument for the bucket of the BigQuery instance.
     */
    private static final String NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME = "num-records";

    /**
     * The name of the argument for the client secrets file to use for authentication with google servers.
     */
    private static final String CLIENT_SECRETS_LOCATION_ARG_NAME = "client-secrets";

    private static final String DEFAULT_PROJECT_ID          = "bigquery-public-data";
    private static final String DEFAULT_FQ_TABLE_ID         = "bigquery-public-data:noaa_lightning.lightning_1987";
    private static final String DEFAULT_BUCKET              = "bigquery-public-data";
    private static final int    DEFAULT_RECORDS_TO_RETRIEVE = 10;
    private static final String DEFAULT_CLIENT_SECRETS_PATH = "~/client_secret.json";

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    @Argument(
            fullName = PROJECT_ID_ARG_LONG_NAME,
            doc = "The project ID of the BigQuery instance.  Defaults to " + DEFAULT_PROJECT_ID)
    private String projectId = DEFAULT_PROJECT_ID;

    @Argument(
            fullName = FQ_TABLE_ID_ARG_LONG_NAME,
            doc = "The fully-qualified table ID of the table containing data from which to query in the BigQuery instance.  Defaults to " + DEFAULT_FQ_TABLE_ID)
    private String fqTableId = DEFAULT_FQ_TABLE_ID;

    @Argument(
            fullName = BUCKET_ARG_LONG_NAME,
            doc = "The bucket containing the BigQuery instance.  Defaults to " + DEFAULT_BUCKET)
    private String bucket = DEFAULT_BUCKET;

    @Argument(
            fullName = NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME,
            doc = "The number of records to retrieve from the BigQuery table.  Defaults to " + DEFAULT_RECORDS_TO_RETRIEVE)
    private int numRecordsToRetrieve = DEFAULT_RECORDS_TO_RETRIEVE;

    @Argument(
            fullName = CLIENT_SECRETS_LOCATION_ARG_NAME,
            doc = "The location of the client secrets json file used for authentication with google servers.  " +
                    "For more information on how to create a client_secret.json file, see the \"Download credentials for API access\" section here:" +
                    "https://cloud.google.com/genomics/docs/how-tos/getting-started ." +
                    "Defaults to " + DEFAULT_CLIENT_SECRETS_PATH
    )
    private String clientSecretsLocation = DEFAULT_CLIENT_SECRETS_PATH;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public void traverse() {

        final BigQuery bigQuery = BigQueryOptions.getDefaultInstance().getService();

        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder(
                        "SELECT "
                                + "CONCAT('https://stackoverflow.com/questions/', CAST(id as STRING)) as url, "
                                + "view_count "
                                + "FROM `bigquery-public-data.stackoverflow.posts_questions` "
                                + "WHERE tags like '%google-bigquery%' "
                                + "ORDER BY favorite_count DESC LIMIT 10")
                        // Use standard SQL syntax for queries.
                        // See: https://cloud.google.com/bigquery/sql-reference/
                        .setUseLegacySql(false)
                        .build();

        // Create a job ID so that we can safely retry.
        final JobId jobId    = JobId.of(UUID.randomUUID().toString());
        Job   queryJob = bigQuery.create(JobInfo.newBuilder(queryConfig).setJobId(jobId).build());

        // Wait for the query to complete.
        try {
            queryJob = queryJob.waitFor();
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        // Check for errors
        if (queryJob == null) {
            throw new GATKException("Query job no longer exists");
        } else if (queryJob.getStatus().getError() != null) {
            // You can also look at queryJob.getStatus().getExecutionErrors() for all
            // errors, not just the latest one.
            throw new RuntimeException(queryJob.getStatus().getError().toString());
        }

        // Get the results.
        final QueryResponse response = bigQuery.getQueryResults(jobId);
        final TableResult result;
        try {
            result = queryJob.getQueryResults();
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        // Log all pages of the results.
        for ( final FieldValueList row : result.iterateAll()) {
            final String url       = row.get("url").getStringValue();
            final long   viewCount = row.get("view_count").getLongValue();
            logger.info("url: %s views: %d", url, viewCount);
        }
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helpers:

}

