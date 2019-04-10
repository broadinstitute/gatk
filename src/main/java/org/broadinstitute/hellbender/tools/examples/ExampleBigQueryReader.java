package org.broadinstitute.hellbender.tools.examples;

import com.google.api.client.http.HttpTransport;
import com.google.api.client.http.javanet.NetHttpTransport;
import com.google.api.client.json.JsonFactory;
import com.google.api.client.json.jackson2.JacksonFactory;
import com.google.api.services.bigquery.Bigquery;
import com.google.api.services.bigquery.model.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.util.List;

/**
 * An example of a class that communicates with BigQuery using the BigQuery Connector with Spark.
 * Created by jonn on 4/9/19.
 */
public class ExampleBigQueryReader extends GATKTool {
    private static final Logger logger = LogManager.getLogger(ExampleBigQueryReader.class);

    //==================================================================================================================
    // Public Static Members:

    /**
     * The name of the argument for the project ID of the BigQuery instance.
     */
    public static final String PROJECT_ID_ARG_LONG_NAME = "project-id";

    /**
     * The name of the argument for the fully qualified input table ID of the BigQuery instance.
     */
    public static final String FQ_TABLE_ID_ARG_LONG_NAME = "fq-table-id";

    /**
     * The name of the argument for the bucket of the BigQuery instance.
     */
    public static final String BUCKET_ARG_LONG_NAME = "bucket";

    /**
     * The name of the argument for the bucket of the BigQuery instance.
     */
    public static final String NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME = "num-records";

    //==================================================================================================================
    // Private Static Members:

    /** Global instances of HTTP transport and JSON factory objects. */
    private static final HttpTransport TRANSPORT    = new NetHttpTransport();
    private static final JsonFactory   JSON_FACTORY = new JacksonFactory();

    //==================================================================================================================
    // Private Members:

    @Argument(
            fullName = PROJECT_ID_ARG_LONG_NAME,
            doc = "The project ID of the BigQuery instance")
    private String projectId;

    @Argument(
            fullName = FQ_TABLE_ID_ARG_LONG_NAME,
            doc = "The fully-qualified table ID of the table containing data from which to query in the BigQuery instance")
    private String fqTableId;

    @Argument(
            fullName = BUCKET_ARG_LONG_NAME,
            doc = "The bucket containing the BigQuery instance")
    private String bucket;

    @Argument(
            fullName = NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME,
            doc = "The number of records to retrieve from the BigQuery table")
    private String numRecordsToRetrieve;


    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public void traverse() {

        try {
            logger.debug("About to create client service for Big Query table.");
            final Bigquery bigQuery = createAuthorizedClient();

            // Start a Query Job
            final String querySql = "SELECT * FROM " + fqTableId + " LIMIT " + numRecordsToRetrieve;
            final JobReference jobId = startQuery(bigQuery, projectId, querySql);

            // Poll for Query Results, return result output
            final Job completedJob = checkQueryResults(bigQuery, projectId, jobId);

            // Return and display the results of the Query Job
            displayQueryResults(bigQuery, projectId, completedJob);

        }
        catch (final IOException ex) {
            throw new GATKException("Error in the connection with BigQuery!", ex);
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while wating for results from BigQuery!", ex);
        }
    }

    //==================================================================================================================
    // Static Methods:

    /**
     * Creates an authorized BigQuery client service using the OAuth 2.0 protocol
     *
     * Note: Adapted from https://github.com/googlearchive/bigquery-samples-java/blob/master/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java
     *
     * This method first creates a BigQuery authorization URL, then prompts the
     * user to visit this URL in a web browser to authorize access. The
     * application will wait for the user to paste the resulting authorization
     * code at the command line prompt.
     *
     * @return an authorized {@link Bigquery} client object
     *
     * @throws IOException
     */
    public static Bigquery createAuthorizedClient() throws IOException {
        return new Bigquery(TRANSPORT, JSON_FACTORY, null);
    }

    //==================================================================================================================
    // Instance Methods:

    /**
     * Display all BigQuery datasets associated with a project
     *
     * Note: Adapted from https://github.com/googlearchive/bigquery-samples-java/blob/master/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java
     *
     * @param bigQuery  an authorized {@link Bigquery} client object
     * @param projectId a string containing the current project ID
     *
     * @throws IOException
     */
    public static void listDatasets(final Bigquery bigQuery,
                                    final String projectId) throws IOException {

        final Bigquery.Datasets.List datasetRequest = bigQuery.datasets().list(projectId);
        final DatasetList            datasetList    = datasetRequest.execute();

        if (datasetList.getDatasets() != null) {
            final List<DatasetList.Datasets> datasets = datasetList.getDatasets();

            logger.info("Available datasets");
            logger.info("----------------");
            logger.info(datasets.toString());

            for (final DatasetList.Datasets dataset : datasets) {
                logger.info(dataset.getDatasetReference().getDatasetId());
            }
        }
    }

    /**
     * Creates a Query {@link JobReference} for a particular query on a dataset.
     *
     * Note: Adapted from https://github.com/googlearchive/bigquery-samples-java/blob/master/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java
     *
     * @param bigQuery  an authorized {@link Bigquery} client object
     * @param projectId a {@link String} containing the current project ID
     * @param querySql  the actual query string
     *
     * @return a reference to the inserted query job
     *
     * @throws IOException
     */
    public static JobReference startQuery(final Bigquery bigQuery,
                                          final String projectId,
                                          final String querySql) throws IOException {
        logger.info("Inserting Query Job: " + querySql);

        final Job                   job         = new Job();
        final JobConfiguration      config      = new JobConfiguration();
        final JobConfigurationQuery queryConfig = new JobConfigurationQuery();
        config.setQuery(queryConfig);

        job.setConfiguration(config);
        queryConfig.setQuery(querySql);

        final Bigquery.Jobs.Insert insert = bigQuery.jobs().insert(projectId, job);
        insert.setProjectId(projectId);

        final JobReference jobId = insert.execute().getJobReference();

        logger.info("Job ID of Query Job is: " + jobId.getJobId());

        return jobId;
    }

    /**
     * Waits for the status of a BigQuery job to be "DONE", then returns the completed {@link Job}.
     *
     * Note: Adapted from https://github.com/googlearchive/bigquery-samples-java/blob/master/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java
     *
     * @param bigQuery  an authorized {@link Bigquery} client object
     * @param projectId a {@link String} containing the current project ID
     * @param jobId     a {@link JobReference} to an inserted query Job
     *
     * @return a reference to the completed {@link Job}.
     *
     * @throws IOException
     * @throws InterruptedException
     */
    private static Job checkQueryResults(final Bigquery bigQuery,
                                         final String projectId,
                                         final JobReference jobId) throws IOException, InterruptedException {

        // Variables to keep track of total query time
        final long startTime = System.currentTimeMillis();
        long elapsedTime;

        while (true) {
            final Job pollJob = bigQuery.jobs().get(projectId, jobId.getJobId()).execute();

            elapsedTime = System.currentTimeMillis() - startTime;
            logger.info("Job status (% 8dms) % 5s: %s", elapsedTime, jobId.getJobId(), pollJob.getStatus().getState());

            if (pollJob.getStatus().getState().equals("DONE")) {
                return pollJob;
            }

            // Pause execution for one second before polling job status again, to
            // reduce unnecessary calls to the BigQUery API and lower overall
            // application bandwidth.
            Thread.sleep(1000);
        }
    }

    /**
     * Displays the results of the given completed query {@link Job} in the given {@link Bigquery} client object
     * with the given {@code projectId}.
     *
     * Note: Adapted from https://github.com/googlearchive/bigquery-samples-java/blob/master/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java
     *
     * @param bigQuery  an authorized {@link Bigquery} client object
     * @param projectId a {@link String} containing the current project ID
     * @param completedJob a query {@link Job} that has already completed
     *
     * @throws IOException
     */
    private static void displayQueryResults(final Bigquery bigQuery,
                                            final String projectId,
                                            final Job completedJob) throws IOException {

        final GetQueryResultsResponse queryResult = bigQuery.jobs()
                .getQueryResults(
                        projectId, completedJob
                                .getJobReference()
                                .getJobId()
                ).execute();

        final List<TableRow> rows = queryResult.getRows();

        logger.info("Query Results:");
        logger.info("------------");
        for (final TableRow row : rows) {
            for (final TableCell field : row.getF()) {
                logger.info("%-50s", field.getV());
            }
            logger.info("");
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}

