package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.storage.v1beta2.*;
import com.google.protobuf.Descriptors;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class PendingBQWriter extends CommittedBQWriter {
    static final Logger logger = LogManager.getLogger(PendingBQWriter.class);

    public PendingBQWriter(BigQueryWriteClient bqWriteClient, String projectId, String datasetName, String tableName) throws Descriptors.DescriptorValidationException, InterruptedException, IOException {
        super(bqWriteClient, projectId, datasetName, tableName, WriteStream.Type.PENDING);
    }

    public void flushBuffer() {
        try {
            writeJsonArray(0);
        } catch (Exception ex) {
            logger.error("Caught exception writing last records on close", ex);
        }
    }

    public void commitWriteStreams() {
//        try {
//            writeJsonArray(0);
//        } catch (Exception ex) {
//            logger.error("Caught exception writing last records on close", ex);
//        }

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
