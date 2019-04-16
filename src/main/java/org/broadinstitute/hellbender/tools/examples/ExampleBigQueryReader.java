package org.broadinstitute.hellbender.tools.examples;

// NOTE:
// Adapted from: https://github.com/googlearchive/bigquery-samples-java/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java

import com.google.cloud.bigquery.*;
import org.apache.ivy.util.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * An example class that communicates with BigQuery using the google bigquery library.
 *
 * By default it will pull from the public NOAA lightning strike database.
 *
 * You may specify the project and fully-qualified table name from which to pull records, as well as the number of
 * records to pull.  For example, the following are all valid invocations of this tool:
 *
 * `gatk ExampleBigQueryReader`
 * `gatk ExampleBigQueryReader`
 *
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
     * The name of the argument for the dataset of the BigQuery instance.
     */
    private static final String DATASET_ARG_LONG_NAME = "dataset";

    /**
     * The name of the argument for the table ID of the BigQuery instance.
     */
    private static final String TABLE_ID_ARG_LONG_NAME = "table-id";

    /**
     * The name of the argument for the bucket of the BigQuery instance.
     */
    private static final String NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME = "num-records";

    private static final String DEFAULT_PROJECT_ID          = "bigquery-public-data";
    private static final String DEFAULT_DATASET             = "noaa_lightning";
    private static final String DEFAULT_TABLE_ID            = "lightning_1987";
    private static final int    DEFAULT_RECORDS_TO_RETRIEVE = 10;

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    @Argument(
            fullName = PROJECT_ID_ARG_LONG_NAME,
            doc = "The project ID of the BigQuery instance.  Defaults to " + DEFAULT_PROJECT_ID)
    private String projectId = DEFAULT_PROJECT_ID;

    @Argument(
            fullName = TABLE_ID_ARG_LONG_NAME,
            doc = "The table ID of the table containing data from which to query in the BigQuery instance.  Defaults to " + DEFAULT_TABLE_ID)
    private String tableId = DEFAULT_TABLE_ID;

    @Argument(
            fullName = DATASET_ARG_LONG_NAME,
            doc = "The dataset containing the table and data from which to query in the BigQuery instance.  Defaults to " + DEFAULT_DATASET)
    private String dataset = DEFAULT_DATASET;

    @Argument(
            fullName = NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME,
            doc = "The number of records to retrieve from the BigQuery table.  Defaults to " + DEFAULT_RECORDS_TO_RETRIEVE)
    private int numRecordsToRetrieve = DEFAULT_RECORDS_TO_RETRIEVE;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public void traverse() {

        final BigQuery bigQuery = BigQueryOptions.getDefaultInstance().getService();

        logger.info( "Constructing query for the first " + numRecordsToRetrieve + " entries in: " + createFQTN() );

        logger.debug( "Query: " + createQueryString() );

        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder( createQueryString() )
                        // Use standard SQL syntax for queries.
                        // See: https://cloud.google.com/bigquery/sql-reference/
                        .setUseLegacySql(false)
                        .build();

        // Create a job ID so that we can safely retry.
        final JobId jobId    = JobId.of(UUID.randomUUID().toString());

        logger.info("Sending query to server...");
        Job   queryJob       = bigQuery.create(JobInfo.newBuilder(queryConfig).setJobId(jobId).build());

        // Wait for the query to complete.
        try {
            logger.info("Waiting for query to complete...");
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
        logger.info("Retrieving query results...");
        final QueryResponse response = bigQuery.getQueryResults(jobId);
        final TableResult result;
        try {
            result = queryJob.getQueryResults();
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        // Log all pages of the results;
        prettyLogResultData(result);
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * @return The fully-qualified table name from which to query data.
     */
    private String createFQTN() {
        return projectId + "." + dataset + "." + tableId;
    }

    /**
     * @return A valid BigTable query string (similar to SQL in syntax) to use to retrieve data.
     */
    private String createQueryString() {
        return "SELECT * FROM `" + createFQTN() + "` LIMIT " + numRecordsToRetrieve;
    }

    /**
     * Logs the given results in a pretty table.
     * @param result A {@link TableResult} object containing the results of a query that generated some data.
     */
    private void prettyLogResultData( final TableResult result ){
        final Schema schema = result.getSchema();

        // Go through all rows and get the length of each column:
        final List<Integer> columnLengths = new ArrayList<>(schema.getFields().size());

        // Start with schema names:
        for ( final Field field : schema.getFields() ) {
            columnLengths.add( field.getName().length() );
        }

        // Check each row:
        for ( final FieldValueList row : result.iterateAll() ) {
            for ( int i = 0; i < row.size() ; ++i ) {
                if ( columnLengths.get(i) < row.get(i).getStringValue().length() ) {
                    columnLengths.set(i, row.get(i).getStringValue().length());
                }
            }
        }

        // Create a separator string for each column:
        final String headerFooter = "+" + columnLengths.stream().map(
                l -> StringUtils.repeat("-", l+2) + "+"
        ).collect(Collectors.joining(""));

        // Now we can log our schema header and rows:
        logger.info( headerFooter );
        logger.info( "|" +
                IntStream.range(0, columnLengths.size()).boxed().map(
                        i -> String.format(" %-"+ columnLengths.get(i) +"s |", schema.getFields().get(i).getName())
                ).collect(Collectors.joining()) );
        logger.info( headerFooter );

        // Log our data:
        for ( final FieldValueList row : result.iterateAll() ) {
            logger.info( "|" +
                    IntStream.range(0, row.size()).boxed().map(
                            i -> String.format(" %-"+ columnLengths.get(i) +"s |", row.get(i).getStringValue())
                    ).collect(Collectors.joining()) );
        }

        logger.info( headerFooter );
    }

    //==================================================================================================================
    // Helpers:

}

