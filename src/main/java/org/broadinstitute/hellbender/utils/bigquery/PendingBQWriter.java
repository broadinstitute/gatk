package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.storage.v1beta2.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

public class PendingBQWriter extends CommittedBQWriter {
    static final Logger logger = LogManager.getLogger(PendingBQWriter.class);

    public PendingBQWriter(String projectId, String datasetName, String tableName) {
        super(projectId, datasetName, tableName, WriteStream.Type.PENDING);
    }

    public void flushBuffer() {
        try {
            if (jsonArr.length() > 0) {
                writeJsonArray();
            }
        } catch (Exception ex) {
            throw new GATKException("Caught exception writing last records on close of " + writeStream.getName(), ex);
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
        if (!commitResponse.hasCommitTime()) {
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
