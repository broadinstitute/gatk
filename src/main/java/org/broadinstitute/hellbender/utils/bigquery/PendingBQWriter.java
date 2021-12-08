package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.storage.v1beta2.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PendingBQWriter extends CommittedBQWriter {
    static final Logger logger = LogManager.getLogger(PendingBQWriter.class);

    public PendingBQWriter(String projectId, String datasetName, String tableName) {
        super(projectId, datasetName, tableName, WriteStream.Type.PENDING);
    }

    public void flushBuffer() {
        try {
            writeJsonArray(0);
        } catch (Exception ex) {
            logger.error("Caught exception writing last records on close", ex);
        }
    }

    public void commitWriteStreams() {
        FinalizeWriteStreamResponse finalizeResponse =
                bqWriteClient.finalizeWriteStream(writeStream.getName());
        logger.info("Rows written: " + finalizeResponse.getRowCount());

        BatchCommitWriteStreamsRequest commitRequest =
                BatchCommitWriteStreamsRequest.newBuilder()
                        .setParent(parentTable.toString())
                        .addWriteStreams(writeStream.getName())
                        .build();
        BatchCommitWriteStreamsResponse commitResponse =
                bqWriteClient.batchCommitWriteStreams(commitRequest);
        // If the response does not have a commit time, it means the commit operation failed.
        if (commitResponse.hasCommitTime() == false) {
            for (StorageError err : commitResponse.getStreamErrorsList()) {
                logger.error(err.getErrorMessage());
            }
            throw new RuntimeException("Error committing the streams");
        }
        logger.info("Appended and committed records successfully.");
    }


    public void close() {
        super.close();
    }
}
