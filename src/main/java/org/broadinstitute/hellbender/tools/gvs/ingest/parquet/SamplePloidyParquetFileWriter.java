package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.tools.gvs.ingest.SamplePloidyWriter;

import java.io.IOException;

public class SamplePloidyParquetFileWriter extends SamplePloidyWriter {
    public SamplePloidyParquetFileWriter(Path file, MessageType schema, CompressionCodecName codecName) {
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public void write(Long sampleId, Long chromosome, Integer ploidy) throws IOException {

    }
}
