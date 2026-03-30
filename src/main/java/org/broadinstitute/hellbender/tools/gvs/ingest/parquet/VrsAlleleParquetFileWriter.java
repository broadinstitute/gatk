package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.tools.gvs.common.VrsAlleleRecord;
import org.broadinstitute.hellbender.tools.gvs.ingest.VrsAlleleWriter;

import java.io.IOException;

/**
 * Parquet writer for the {@code vrs_allele} table.
 */
public class VrsAlleleParquetFileWriter extends AbstractParquetFileWriter implements VrsAlleleWriter {

    public VrsAlleleParquetFileWriter(Path file, MessageType schema, CompressionCodecName codecName) throws IOException {
        super(file, schema, codecName);
    }

    @Override
    public void write(VrsAlleleRecord record) throws IOException {
        super.write(toJson(record));
    }
}
