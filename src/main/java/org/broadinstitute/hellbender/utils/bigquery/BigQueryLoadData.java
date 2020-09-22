package org.broadinstitute.hellbender.utils.bigquery;

import com.google.api.client.util.Charsets;
import com.google.cloud.bigquery.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import shaded.cloud_nio.com.google.common.io.Files;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.util.concurrent.TimeoutException;

public class BigQueryLoadData {
    private static final Logger logger = LogManager.getLogger(BigQueryLoadData.class);

    /**
     * Example of creating a dataset.
     */
    // [TARGET create(DatasetInfo, DatasetOption...)]
    // [VARIABLE "my_dataset_name"]
    public static Dataset createDataset(BigQuery bigquery, String datasetName) {
        // [START bigquery_create_dataset]
        Dataset dataset = null;
        DatasetInfo datasetInfo = DatasetInfo.newBuilder(datasetName).build();
        try {
            // the dataset was created
            dataset = bigquery.create(datasetInfo);
        } catch (BigQueryException e) {
            // the dataset was not created
        }
        // [END bigquery_create_dataset]
        return dataset;
    }

    /**
     * Example of creating a table.
     */
    // [TARGET create(TableInfo, TableOption...)]
    // [VARIABLE "my_dataset_name"]
    // [VARIABLE "my_table_name"]
    // [VARIABLE "string_field"]
    public static Table createTable(BigQuery bigquery, String datasetName, String tableName, String fieldName) {
        // [START bigquery_create_table]
        TableId tableId = TableId.of(datasetName, tableName);
        // Table field definition
        Field field = Field.of(fieldName, LegacySQLTypeName.STRING);
        // Table schema definition
        Schema schema = Schema.of(field);
        TableDefinition tableDefinition = StandardTableDefinition.of(schema);
        TableInfo tableInfo = TableInfo.newBuilder(tableId, tableDefinition).build();
        Table table = bigquery.create(tableInfo);
        // [END bigquery_create_table]
        return table;
    }

    /**
     * Example of creating a channel with which to write to a table.
     */
    // [TARGET writer(WriteChannelConfiguration)]
    // [VARIABLE "my_dataset_name"]
    // [VARIABLE "my_table_name"]
    // [VARIABLE "StringValue1\nStringValue2\n"]
    public static long writeFileToTable(BigQuery bigquery, String datasetName, String tableName, File csvFile) {
        long numRows = 0;
        try {
            TableId tableId = TableId.of(datasetName, tableName);
            WriteChannelConfiguration writeChannelConfiguration =
                    WriteChannelConfiguration.newBuilder(tableId).setFormatOptions(FormatOptions.csv()).build();
            TableDataWriteChannel writer = bigquery.writer(writeChannelConfiguration);
            // Write data to writer
            try (OutputStream stream = Channels.newOutputStream(writer)) {
                Files.copy(csvFile, stream);
            }
            // Get load job
            Job job = writer.getJob();
            job = job.waitFor();
            JobStatistics.LoadStatistics stats = job.getStatistics();
            numRows = stats.getOutputRows();
        } catch (Exception ex) {
            logger.warn("big query load job interrupted", ex);

        }
        return numRows;
        // [END ]
    }



    public static long loadData(BigQuery bigquery, String datasetName, String tableName, String location, File csvFile) {
        long rowsLoaded = 0;
        TableId tableId = TableId.of(datasetName, tableName);
        WriteChannelConfiguration writeChannelConfiguration =
                WriteChannelConfiguration.newBuilder(tableId).setFormatOptions(FormatOptions.csv()).build();
        // The location must be specified; other fields can be auto-detected.
        JobId jobId = JobId.newBuilder().setLocation(location).build();
        TableDataWriteChannel writer = bigquery.writer(jobId, writeChannelConfiguration);
        // Write data to writer
        try {
            try (OutputStream stream = Channels.newOutputStream(writer)) {
                Files.copy(csvFile, stream);
            }
            // Get load job
            Job job = writer.getJob();
            job = job.waitFor();
            JobStatistics.LoadStatistics stats = job.getStatistics();
            rowsLoaded = stats.getOutputRows();
        } catch (InterruptedException | IOException ex) {
            logger.warn("big query load job interrupted", ex);
        }
        return rowsLoaded;
    }

    /**
     * Example of deleting a dataset from its id, even if non-empty.
     */
    // [TARGET delete(String, DatasetDeleteOption...)]
    // [VARIABLE "my_dataset_name"]
    public static boolean deleteDataset(BigQuery bigquery, String datasetName) {
        // [START ]
        boolean deleted = bigquery.delete(datasetName, BigQuery.DatasetDeleteOption.deleteContents());
        if (deleted) {
            // the dataset was deleted
        } else {
            // the dataset was not found
        }
        // [END ]
        return deleted;
    }

    /**
     * Example of deleting a table.
     */
    // [TARGET delete(String, String)]
    // [VARIABLE "my_dataset_name"]
    // [VARIABLE "my_table_name"]
    public static boolean deleteTable(BigQuery bigquery, TableId tableId) {
        // [START ]
        boolean deleted = bigquery.delete(tableId);
        if (deleted) {
            // the table was deleted
        } else {
            // the table was not found
        }
        // [END ]
        return deleted;
    }



}
