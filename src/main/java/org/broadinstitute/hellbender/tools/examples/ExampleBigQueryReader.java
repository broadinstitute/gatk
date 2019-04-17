package org.broadinstitute.hellbender.tools.examples;

// NOTE:
// Adapted from: https://github.com/googlearchive/bigquery-samples-java/src/main/java/com/google/cloud/bigquery/samples/BigQueryJavaGettingStarted.java

import com.google.cloud.bigquery.TableResult;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;

/**
 * An example class that communicates with BigQuery using the google bigquery library.
 *
 * By default it will pull from the public NOAA lightning strike database.
 *
 * You may specify the project and fully-qualified table name from which to pull records, as well as the number of
 * records to pull.  For example, the following are all valid invocations of this tool:
 *
 * `gatk ExampleBigQueryReader`
 * `gatk ExampleBigQueryReader --project-id broad-dsp-spec-ops --dataset gcp_joint_genotyping --table-id chr2_sample_100_new_way`
 *
 * Created by jonn on 4/9/19.
 */
@CommandLineProgramProperties(
        summary = "Example tool that communicates with BigQuery using the google bigquery library.  Will print the first several rows from the given BigQuery table.",
        oneLineSummary = "Example tool that prints the first few rows of a given BigQuery table.",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = false
)
@DocumentedFeature
public class ExampleBigQueryReader extends CommandLineProgram {
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
            doc = "The table ID of the table containing data from which to query in the BigQuery instance.")
    private String tableId = DEFAULT_TABLE_ID;

    @Argument(
            fullName = DATASET_ARG_LONG_NAME,
            doc = "The dataset containing the table and data from which to query in the BigQuery instance.")
    private String dataset = DEFAULT_DATASET;

    @Argument(
            fullName = NUM_RECORDS_TO_RETRIEVE_ARG_LONG_NAME,
            doc = "The number of records to retrieve from the BigQuery table.")
    private int numRecordsToRetrieve = DEFAULT_RECORDS_TO_RETRIEVE;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public Object doWork() {

        logger.info( "Constructing query for the first " + numRecordsToRetrieve + " entries in: " + createFQTN() );

        final String queryString = createQueryString();
        logger.debug( "Query: " + createQueryString() );

        // Execute the query against BigQuery using the default BigQuery connection information:
        final TableResult result = BigQueryUtils.executeQuery( queryString );

        // Log all pages of the results;
        BigQueryUtils.logResultDataPretty( result, logger );

        // No-op return value:
        return null;
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

    //==================================================================================================================
    // Helpers:

}

