package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.gax.rpc.ServerStream;
import com.google.cloud.bigquery.*;
import com.google.cloud.bigquery.storage.v1beta1.*;
import com.google.common.base.Preconditions;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.avro.file.DataFileReader;
import org.apache.avro.file.DataFileStream;
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
import org.broadinstitute.hellbender.tools.evoquer.GATKAvroReader;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

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
        return executeQuery(getBigQueryEndPoint(), queryString, false);
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @param runQueryInBatchMode If true, run the query in batch mode, which is lower priority but has no limit on the number of concurrent queries
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final String queryString, final boolean runQueryInBatchMode) {
        return executeQuery(getBigQueryEndPoint(), queryString, runQueryInBatchMode);
    }

    /**
     * Executes the given {@code queryString} on the default instance of {@link BigQuery} as created by {@link #getBigQueryEndPoint()}.
     * Will block until results are returned.
     * For more information on querying BigQuery tables, see: https://cloud.google.com/bigquery/sql-reference/
     * @param bigQuery The {@link BigQuery} instance against which to execute the given {@code queryString}.
     * @param queryString The {@link BigQuery} query string to execute.  Must use standard SQL syntax.  Must contain the project ID, data set, and table name in the `FROM` clause for the table from which to retrieve data.
     * @param runQueryInBatchMode If true, run the query in batch mode, which is lower priority but has no limit on the number of concurrent queries
     * @return A {@link TableResult} object containing the results of the query executed.
     */
    public static TableResult executeQuery(final BigQuery bigQuery, final String queryString, final boolean runQueryInBatchMode) {

        // Create a query configuration we can run based on our query string:
        final QueryJobConfiguration queryConfig =
                QueryJobConfiguration.newBuilder( queryString )
                        .setUseLegacySql(false)
                        .setPriority(runQueryInBatchMode ? QueryJobConfiguration.Priority.BATCH : QueryJobConfiguration.Priority.INTERACTIVE)
                        .build();

        logger.info("Executing Query: \n\n" + queryString);
        final TableResult result = submitQueryAndWaitForResults( bigQuery, queryConfig );
        logger.info("Query returned " + result.getTotalRows() + " results.");
        return result;
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
        }
        catch (final InterruptedException ex) {
            throw new GATKException("Interrupted while waiting for query job to complete", ex);
        }

        return result;
    }

    public static StorageAPIAvroReader executeQueryWithStorageAPI(final String queryString) {

        return executeQueryWithStorageAPI(queryString, false);
    }

    public static StorageAPIAvroReader executeQueryWithStorageAPI(final String queryString, final boolean runQueryInBatchMode) {
        final String tempTableProject = "broad-dsp-spec-ops";
        final String tempTableDataset = "temp_tables";
        final String tempTableName = UUID.randomUUID().toString().replace('-', '_');
        final String tempTableFullyQualified = String.format("%s.%s.%s", tempTableProject, tempTableDataset, tempTableName);

        final String queryStringIntoTempTable = "CREATE TABLE `" + tempTableFullyQualified + "`\n" +
                "OPTIONS(\n" +
                "  expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 1 DAY)\n" +
                ") AS\n" +
                queryString;

        final TableResult result = executeQuery(queryStringIntoTempTable, runQueryInBatchMode);

        final List<String> allResultFields = Arrays.asList("position", "values");

        return new StorageAPIAvroReader(tempTableProject, tempTableDataset, tempTableName, allResultFields);
    }

    public static class StorageAPIAvroReader implements GATKAvroReader {

        private static int rowCount = 0;

        private BigQueryStorageClient client;

        private Iterator<Storage.ReadRowsResponse> serverStream;

        private org.apache.avro.Schema schema;

        private DatumReader<GenericRecord> datumReader;

        // Decoder object will be reused to avoid re-allocation and too much garbage collection.
        private BinaryDecoder decoder = null;

        private AvroProto.AvroRows currentAvroRows;

        // GenericRecord object will be reused.
        private GenericRecord nextRow = null;

        public StorageAPIAvroReader( final String tableProject,
                                     final String tableDataset,
                                     final String tableName,
                                     final List<String> fields) {

            try {
                this.client = BigQueryStorageClient.create();

                final String parent = String.format("projects/%s", tableProject);

                final TableReferenceProto.TableReference tableReference = TableReferenceProto.TableReference.newBuilder()
                        .setProjectId(tableProject)
                        .setDatasetId(tableDataset)
                        .setTableId(tableName)
                        .build();

                final ReadOptions.TableReadOptions tableReadOptions =
                        ReadOptions.TableReadOptions.newBuilder()
                                .addAllSelectedFields(fields)
                                .build();

                final Storage.CreateReadSessionRequest.Builder builder = Storage.CreateReadSessionRequest.newBuilder()
                        .setParent(parent)
                        .setTableReference(tableReference)
                        .setReadOptions(tableReadOptions)
                        .setRequestedStreams(1)
                        .setFormat(Storage.DataFormat.AVRO);

                final Storage.ReadSession session = client.createReadSession(builder.build());
                Preconditions.checkState(session.getStreamsCount() > 0);

                this.schema = new org.apache.avro.Schema.Parser().parse(session.getAvroSchema().getSchema());

                this.datumReader = new GenericDatumReader<>(
                        new org.apache.avro.Schema.Parser().parse(session.getAvroSchema().getSchema()));

                // Use the first stream to perform reading.
                Storage.StreamPosition readPosition = Storage.StreamPosition.newBuilder()
                        .setStream(session.getStreams(0))
                        .build();

                Storage.ReadRowsRequest readRowsRequest = Storage.ReadRowsRequest.newBuilder()
                        .setReadPosition(readPosition)
                        .build();

                this.serverStream = client.readRowsCallable().call(readRowsRequest).iterator();

                loadNextRow();
            } catch ( IOException e ) {
                throw new GATKException("I/O Error", e);
            }
        }

        private void loadNextRow() {
            try {
                if ( decoder != null && ! decoder.isEnd() ) {
                    nextRow = datumReader.read(null, decoder);
                } else {
                    fetchNextAvroRows();

                    if ( decoder != null && ! decoder.isEnd() ) {
                        nextRow = datumReader.read(null, decoder);
                    } else {
                        nextRow = null; // end of traversal
                    }
                }
            } catch ( IOException e ) {
                throw new GATKException("I/O error", e);
            }
        }

        private void fetchNextAvroRows() {
            if ( serverStream.hasNext() ) {
                currentAvroRows = serverStream.next().getAvroRows();
                decoder = DecoderFactory.get()
                        .binaryDecoder(currentAvroRows.getSerializedBinaryRows().toByteArray(), decoder);
            } else {
                currentAvroRows = null;
                decoder = null;
            }
        }

        @Override
        public org.apache.avro.Schema getSchema() {
            return schema;
        }

        /*public void processRows( AvroProto.AvroRows avroRows ) throws IOException {
            decoder = DecoderFactory.get()
                    .binaryDecoder(avroRows.getSerializedBinaryRows().toByteArray(), decoder);

            while ( !decoder.isEnd() ) {
                rowCount++;
                // Reusing object row
                row = datumReader.read(row, decoder);

                System.out.println("Schema:");
                row.getSchema().getFields().forEach(field -> System.out.println(field.name()));
                System.out.println("Values schema:");
                final List<org.apache.avro.Schema.Field> valuesFields = row.getSchema().getField("values").schema().getElementType().getFields(); //.getFields().forEach(field -> System.out.println(field.name()));

                System.out.println("Position: " + row.get("position"));
                row.getSchema().getFields().forEach(field -> System.out.println(field.name() + " " + row.get(field.name())));

                System.out.println("values class: " + row.get("values").getClass());
                final GenericData.Array<?> valArray = (GenericData.Array) row.get("values");
                for ( final Object valArrayItem : valArray ) {
                    System.out.println(valArrayItem.getClass());
                    final GenericData.Record valRec = (GenericData.Record) valArrayItem;
                    for ( final org.apache.avro.Schema.Field valField : valRec.getSchema().getFields() ) {
                        System.out.println(valField.name() + ": " + valRec.get(valField.name()));
                    }
                }

                System.exit(0);
                if ( rowCount % 100 == 0 ) {
                    System.out.println("Processed " + rowCount + " rows");
                }
            }
        } */

        @Override
        public Iterator<GenericRecord> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            return nextRow != null;
        }

        @Override
        public GenericRecord next() {
            if ( ! hasNext() ) {
                throw new NoSuchElementException("next() called when ! hasNext()");
            }

            final GenericRecord recordToReturn = nextRow;
            loadNextRow();
            return recordToReturn;
        }

        @Override
        public void close() {
            client.close();
        }
    }

    private static class SimpleRowReader {

        private static int rowCount = 0;

        private final DatumReader<GenericRecord> datumReader;

        // Decoder object will be reused to avoid re-allocation and too much garbage collection.
        private BinaryDecoder decoder = null;

        // GenericRecord object will be reused.
        private GenericRecord row = null;

        public SimpleRowReader( org.apache.avro.Schema schema ) {
            Preconditions.checkNotNull(schema);
            datumReader = new GenericDatumReader<>(schema);
        }

        /**
         * Sample method for processing AVRO rows which only validates decoding.
         *
         * @param avroRows object returned from the ReadRowsResponse.
         */
        public void processRows( AvroProto.AvroRows avroRows ) throws IOException {
            decoder = DecoderFactory.get()
                    .binaryDecoder(avroRows.getSerializedBinaryRows().toByteArray(), decoder);

            while ( !decoder.isEnd() ) {
                rowCount++;
                // Reusing object row
                row = datumReader.read(row, decoder);

                System.out.println("Schema:");
                row.getSchema().getFields().forEach(field -> System.out.println(field.name()));
                System.out.println("Values schema:");
                final List<org.apache.avro.Schema.Field> valuesFields = row.getSchema().getField("values").schema().getElementType().getFields(); //.getFields().forEach(field -> System.out.println(field.name()));

                System.out.println("Position: " + row.get("position"));
                row.getSchema().getFields().forEach(field -> System.out.println(field.name() + " " + row.get(field.name())));

                System.out.println("values class: " + row.get("values").getClass());
                final GenericData.Array<?> valArray = (GenericData.Array) row.get("values");
                for ( final Object valArrayItem : valArray ) {
                    System.out.println(valArrayItem.getClass());
                    final GenericData.Record valRec = (GenericData.Record) valArrayItem;
                    for ( final org.apache.avro.Schema.Field valField : valRec.getSchema().getFields() ) {
                        System.out.println(valField.name() + ": " + valRec.get(valField.name()));
                    }
                }

                System.exit(0);
                if ( rowCount % 100 == 0 ) {
                    System.out.println("Processed " + rowCount + " rows");
                }
            }
        }
    }
}
