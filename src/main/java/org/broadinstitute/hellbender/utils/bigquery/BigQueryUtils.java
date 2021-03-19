package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.*;
import com.google.cloud.bigquery.storage.v1.*;
import com.google.common.base.Preconditions;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.BinaryDecoder;
import org.apache.avro.io.DatumReader;
import org.apache.avro.io.DecoderFactory;
import org.apache.ivy.util.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

/**
 * Utility class for dealing with BigQuery connections / tables / queries /etc.
 *
 * Created by jonn on 4/17/19.
 */
public final class BigQueryUtils {
    private BigQueryUtils(){}

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
     * @param executionProjectId The google project that should be used to execute this query
     *
     * @return A {@link BigQuery} object that can be used to interact with a BigQuery data set.
     */
    public static BigQuery getBigQueryEndPoint(String executionProjectId) {
        return (executionProjectId != null) ? BigQueryOptions.newBuilder().setProjectId(executionProjectId).build().getService() : getBigQueryEndPoint();
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @param labels The {@link BigQuery} label to add the job run.  Must use Map<String, String>. Can be null to indicate no labels.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final String queryString, final Map<String, String> labels) {
        return executeQuery(getBigQueryEndPoint(), queryString, false, labels );
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @param runQueryInBatchMode If true, run the query in batch mode, which is lower priority but has no limit on the number of concurrent queries
     * @param labels The {@link BigQuery} label to add the job run.  Must use Map<String, String>. Can be null to indicate no labels.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final String queryString, final boolean runQueryInBatchMode, final Map<String, String> labels) {
        return executeQuery(getBigQueryEndPoint(), queryString, runQueryInBatchMode, labels);
    }

    /**
     * Executes the given {@code queryString} on the provided BigQuery instance.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param bigQuery The {@link BigQuery} instance against which to execute the given {@code queryString}.
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @param runQueryInBatchMode If true, run the query in batch mode, which is lower priority but has no limit on the number of concurrent queries
     * @param labels The {@link BigQuery} label to add the job run.  Must use Map<String, String>. Can be null to indicate no labels.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final BigQuery bigQuery, final String queryString, final boolean runQueryInBatchMode, final Map<String, String> labels) {

        // Create a query configuration we can run based on our query string:
        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder( queryString )
                        .setUseLegacySql(false)
                        .setPriority(runQueryInBatchMode ? QueryJobConfiguration.Priority.BATCH : QueryJobConfiguration.Priority.INTERACTIVE)
                        .setLabels(labels)
                        .build();

        logger.info("Executing Query: \n\n" + queryString);
        final TableResult result = submitQueryAndWaitForResults( bigQuery, queryConfig );
        logger.info("Query returned " + result.getTotalRows() + " results.");
        return result;
    }

    /**
     * Executes the given {@code queryString} on the provided BigQuery instance.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param bigQuery The {@link BigQuery} instance against which to execute the given {@code queryString}.  Must contain the table name in the `FROM` clause for the table from which to retrieve data.
     * @param projectID The BigQuery {@code project ID} containing the {@code dataSet} and table from which to query data.
     * @param dataSet The BigQuery {@code dataSet} containing the table from which to query data.
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table ID in the `FROM` clause for the table from which to retrieve data.
     * @param labels The {@link BigQuery} label to add the job run.  Must use Map<String, String>. Can be null to indicate no labels.
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final BigQuery bigQuery,
                                           final String projectID,
                                           final String dataSet,
                                           final String queryString,
                                           final Map<String, String> labels) {

        // Create a query configuration we can run based on our query string:
        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder( queryString )
                        .setUseLegacySql(false)
                        .setDefaultDataset(DatasetId.of(projectID, dataSet))
                        .setLabels(labels)
                        .build();

        return submitQueryAndWaitForResults( bigQuery, queryConfig );
    }

    /**
     * Creates a {@link String} containing the given results in a pretty ascii-art table.
     * @param result A {@link TableResult} object containing the results of a query that generated some data.
     * @return A {@link String} containing the contents of the given query result pretty-printed in an ascii-art table.
     */
    public static String getResultDataPrettyString(final TableResult result){
        final Schema schema = result.getSchema();

        final List<Integer> columnWidths = calculateColumnWidths( result );
        final boolean rowsAllPrimitive =
                StreamSupport.stream(result.iterateAll().spliterator(), false)
                        .flatMap( row -> row.stream().map(v -> v.getAttribute() == FieldValue.Attribute.PRIMITIVE) )
                        .allMatch( b -> b );

        // Create a separator string for the header and boarders:
        final String headerFooter = "+" + columnWidths.stream().map(
                l -> StringUtils.repeat("=", l+2) + "+"
        ).collect(Collectors.joining(""));

        // Create a Row Separator:
        final String rowSeparator = headerFooter.replace('=', '-');

        // Create a string builder to keep the pretty table:
        final StringBuilder sb = new StringBuilder();

        // Now we can write our schema header and rows:
        addHeaderToStringBuilder(schema, columnWidths, headerFooter, sb);

        // Write our data to the string builder:
        for ( final FieldValueList row : result.iterateAll() ) {

            // If the row fields are all simple, then we can do the simple thing:
            if ( rowsAllPrimitive ) {
                addPrimitiveRowToStringBuilder(row, columnWidths, sb);
            }
            else {
                addComplexRowToStringBuilder(row, schema, columnWidths, sb);
                sb.append(rowSeparator);
            }
            sb.append("\n");
        }

        sb.append( headerFooter );
        sb.append("\n");

        return sb.toString();
    }

    private static void addHeaderToStringBuilder(final Schema schema, final List<Integer> columnWidths, final String headerFooter, final StringBuilder sb) {
        sb.append( headerFooter );
        sb.append("\n");
        sb.append("|");
        sb.append(
                IntStream.range(0, columnWidths.size()).boxed().map(
                        i -> String.format(" %-"+ columnWidths.get(i) +"s |", schema.getFields().get(i).getName())
                ).collect(Collectors.joining()) );
        sb.append("\n");
        sb.append( headerFooter );
        sb.append("\n");
    }

    private static void addComplexRowToStringBuilder(final FieldValueList row, final Schema schema, final List<Integer> columnWidths, final StringBuilder sb) {

        // TODO: Clean this up... Probably make a getStringValue(FIELD) method to use here, addPrimitiveRowToStringBuilder, and calculateColumnWidths

        // For fields that have multiple values, we need to do something special.
        // In fact, we need to go through each value of each row and track how many fields it has.
        int maxNumValuesInRow = 1;
        for ( int i = 0; i < row.size(); ++i ) {
            final FieldValue value = row.get(schema.getFields().get(i).getName());
            if ( !value.isNull() && (value.getAttribute() != FieldValue.Attribute.PRIMITIVE) ) {
                if (maxNumValuesInRow <= value.getRecordValue().size()) {
                    maxNumValuesInRow = value.getRecordValue().size();
                }
            }
        }

        for ( int currentFieldNum = 0; currentFieldNum < maxNumValuesInRow ; ++currentFieldNum ) {
            sb.append("|");
            for ( int i = 0; i < row.size(); ++i ) {
                final FieldValue value = row.get(i);
                if ( value.isNull() ) {
                    sb.append(getEmptyColumnString(columnWidths, i));
                }
                else if (value.getAttribute() == FieldValue.Attribute.PRIMITIVE ) {
                    if ( currentFieldNum == 0 ) {
                        sb.append( String.format(" %-" + columnWidths.get(i) + "s |", value.getStringValue()) );
                    }
                    else {
                        sb.append(getEmptyColumnString(columnWidths, i));
                    }
                }
                else {
                    if ( value.getRepeatedValue().size() == 0 ) {
                        sb.append(getEmptyColumnString(columnWidths, i));
                    }
                    else if ( currentFieldNum < value.getRepeatedValue().size() ) {
                        if ( value.getRepeatedValue().get(currentFieldNum).getAttribute() == FieldValue.Attribute.PRIMITIVE ) {
                            sb.append(String.format(" %-" + columnWidths.get(i) + "s |",
                                    value.getRepeatedValue().get(currentFieldNum).getStringValue())
                            );
                        }
                        else {
                            // This is kind of gross, but it seems to be the only way to get a particular
                            // value of a field that is in an array:
                            sb.append(String.format(" %-" + columnWidths.get(i) + "s |",
                                    value.getRepeatedValue().get(currentFieldNum).getRecordValue().get(0).getStringValue())
                            );
                        }
                    }
                    else {
                        sb.append(getEmptyColumnString(columnWidths, i));
                    }
                }
            }
            sb.append("\n");
        }
    }

    private static String getEmptyColumnString(final List<Integer> columnWidths, final int i) {
        return String.format(" %-" + columnWidths.get(i) + "s |", "");
    }

    private static void addPrimitiveRowToStringBuilder(final FieldValueList row, final List<Integer> columnWidths, final StringBuilder sb) {
        sb.append("|");
        sb.append(IntStream.range(0, row.size()).boxed().map(
                i -> String.format(" %-" + columnWidths.get(i) + "s |", row.get(i).getStringValue())
        ).collect(Collectors.joining()));
    }

    /**
     * Creates a {@link String} containing the given results in a pretty ascii-art table.
     * @param result A {@link TableResult} object containing the results of a query that generated some data.
     * @param theLogger A {@link Logger} object with which to log the results contained in {@code result}.
     */
    public static void logResultDataPretty( final TableResult result, final Logger theLogger ){
        for ( final String line : getResultDataPrettyString(result).split("\n") ) {
            theLogger.info( line );
        }
    }

    //==================================================================================================================
    // Helper Methods:

    private static List<Integer> calculateColumnWidths( final TableResult result ) {
        // Go through all rows and get the length of each column:
        final List<Integer> columnWidths = new ArrayList<>(result.getSchema().getFields().size());

        // Start with schema names:
        for ( final Field field : result.getSchema().getFields() ) {
            columnWidths.add( field.getName().length() );
        }

        // Check each row and each row's array values (if applicable):
        for ( final FieldValueList row : result.iterateAll() ) {
            for ( int i = 0; i < row.size() ; ++i ) {
                // Only get the row size if it's not null:
                if ( !row.get(i).isNull() ) {
                    if ( row.get(i).getAttribute() == FieldValue.Attribute.PRIMITIVE ) {
                        if ( columnWidths.get(i) < row.get(i).getStringValue().length() ) {
                            columnWidths.set(i, row.get(i).getStringValue().length());
                        }
                    }
                    else {
                        for ( int j = 0; j < row.get(i).getRepeatedValue().size(); ++j ) {
                            final String stringValue;
                            if ( row.get(i).getRepeatedValue().get(j).getAttribute() == FieldValue.Attribute.PRIMITIVE ) {
                                stringValue = row.get(i).getRepeatedValue().get(j).getStringValue();
                            }
                            else {
                                stringValue = row.get(i).getRepeatedValue().get(j).getRecordValue().get(0).getStringValue();
                            }
                            if ( columnWidths.get(i) < stringValue.length() ) {
                                columnWidths.set(i, stringValue.length());
                            }
                        }
                    }
                }
            }
        }

        return columnWidths;
    }

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

            long bytesProcessed = ((JobStatistics.QueryStatistics) queryJob.getStatistics()).getTotalBytesProcessed();
            logger.info(String.format("%.2f MB actually scanned", bytesProcessed / 1000000.0));

        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        return result;
    }

    private static long getQueryCostBytesProcessedEstimate(String queryString) {
        final QueryJobConfiguration dryRunQueryConfig =
                QueryJobConfiguration.newBuilder( queryString )
                        .setUseLegacySql(false)
                        .setDryRun(true)
                        .setUseQueryCache(false)
                        .setPriority(QueryJobConfiguration.Priority.INTERACTIVE)
                        .build();

        Job dryRunJob = getBigQueryEndPoint().create(JobInfo.newBuilder(dryRunQueryConfig).build());
        long bytesProcessed = ((JobStatistics.QueryStatistics) dryRunJob.getStatistics()).getTotalBytesProcessed();
        return bytesProcessed;
    }

    public static StorageAPIAvroReader executeQueryWithStorageAPI(final String queryString, final List<String> fieldsToRetrieve, final String projectID, Map<String, String> labels) {

        return executeQueryWithStorageAPI(queryString, fieldsToRetrieve, projectID, false, labels);
    }

    public static StorageAPIAvroReader executeQueryWithStorageAPI(final String queryString, final List<String> fieldsToRetrieve, final String projectID, final boolean runQueryInBatchMode,  Map<String, String> labels) {
        final String tempTableDataset = "temp_tables";
        final String tempTableName = UUID.randomUUID().toString().replace('-', '_');
        final String tempTableFullyQualified = String.format("%s.%s.%s", projectID, tempTableDataset, tempTableName);

        long bytesProcessed = getQueryCostBytesProcessedEstimate(queryString);
        logger.info(String.format("Estimated %s MB scanned", bytesProcessed/1000000));

        final String queryStringIntoTempTable = "CREATE TABLE `" + tempTableFullyQualified + "`\n" +
                "OPTIONS(\n" +
                "  expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 1 DAY)\n" +
                ") AS\n" +
                queryString;

        executeQuery(queryStringIntoTempTable, runQueryInBatchMode, labels);

        final Table tableInfo = getBigQueryEndPoint().getTable( TableId.of(projectID, tempTableDataset, tempTableName) );
        logger.info(String.format("Query temp table created with %s rows and %s bytes in size", tableInfo.getNumRows(), tableInfo.getNumBytes()));

        TableReference tr = new TableReference(tempTableFullyQualified, fieldsToRetrieve);

        return new StorageAPIAvroReader(tr);
    }
}
