package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import org.broadinstitute.hellbender.tools.gvs.ingest.VetWriter;
import org.json.JSONObject;

import java.io.IOException;
import java.util.List;

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

    @Override
    public void write(long location, VariantContext variant, long sampleId, List<String> vrsAlleleIds) throws IOException {
        JSONObject jsonObject = createJson(location, variant, sampleId, skipLoadingVqsrFields, forceLoadingFromNonAlleleSpecific, vrsAlleleIds);
        this.write(jsonObject);
    }
}
