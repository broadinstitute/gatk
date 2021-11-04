package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.io.DatumWriter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.bigquery.PendingBQWriter;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class RefRangesBQWriter extends RefRangesWriter {
    private PendingBQWriter bqWriter;

    public RefRangesBQWriter(String projectId, String datasetName, String tableName) throws IOException{
        bqWriter = new PendingBQWriter(projectId, datasetName, tableName);
    }

    @Override
    public void write(long location, long sampleId, int length, String state) throws IOException {
        try {
            bqWriter.addJsonRow(createJsonRow(location, sampleId, length, state));
        } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
            throw new IOException("BQ exception", ex);
        }
    }

    private JSONObject createJsonRow(long location, long sampleId, int length, String state) {
        JSONObject record = new JSONObject();
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("length", length);
        record.put("state", state);
        return record;
    }

    public void commitData() {
        bqWriter.flushBuffer();
        bqWriter.commitWriteStreams();
    }

    public void close() throws IOException {
        bqWriter.close();
    }
}
