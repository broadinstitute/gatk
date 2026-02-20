package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.tools.gvs.ingest.VetWriter;
import org.json.JSONObject;

import java.io.IOException;

/**
 * Parquet writer for variant data.
 * <p>
 * This writer stores variant information in Parquet format.
 */
public class VetParquetFileWriter extends AbstractParquetFileWriter implements VetWriter {

    private final boolean skipLoadingVqsrFields;
    private final boolean forceLoadingFromNonAlleleSpecific;

    public VetParquetFileWriter(Path file, boolean skipLoadingVqsrFields, boolean forceLoadingFromNonAlleleSpecific, MessageType schema, CompressionCodecName codecName) throws IOException {
        super(file, schema, codecName);
        this.skipLoadingVqsrFields = skipLoadingVqsrFields;
        this.forceLoadingFromNonAlleleSpecific = forceLoadingFromNonAlleleSpecific;
    }

    /**
     * Writes a JSON object to the Parquet file.
     *
     * @param object the JSON object to write
     * @throws IOException if an error occurs during writing
     */
    public void write(JSONObject object) throws IOException {
        super.write(object);
    }

    @Override
    public void write(long location, VariantContext variant, long sampleId) throws IOException {
        JSONObject jsonObject = createJson(location, variant, sampleId, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific);
        this.write(jsonObject);
    }
}
