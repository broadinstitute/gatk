package org.broadinstitute.hellbender.tools.gvs.ingest.bq;

import com.google.cloud.bigquery.BigQuery;
import com.google.cloud.bigquery.FieldList;
import com.google.cloud.bigquery.Schema;
import com.google.cloud.bigquery.Table;
import com.google.cloud.bigquery.TableId;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.IngestUtils;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.tools.gvs.ingest.RefRangesWriter;
import org.broadinstitute.hellbender.utils.gvs.bigquery.BigQueryUtils;
import org.json.JSONObject;

import java.io.IOException;

public class RefRangesBQWriter extends AbstractBQWriter implements RefRangesWriter {

    private final static String REF_RANGES_FILETYPE_PREFIX = "ref_ranges_";

    public RefRangesBQWriter(String projectId, String datasetName, String tableName) throws IOException {
        super(projectId, datasetName, tableName);
    }

    public static void sanityCheckRefRangesSchemaForCompressedReferences(String projectID, String datasetName, Long sampleId, boolean storeCompressedReferences) {
        final BigQuery bigquery = BigQueryUtils.getBigQueryEndPoint(projectID);
        final int sampleTableNumber = IngestUtils.getTableNumber(sampleId, IngestConstants.partitionPerTable);
        final String tableName = REF_RANGES_FILETYPE_PREFIX + String.format("%03d", sampleTableNumber);
        final TableId tableId = TableId.of(projectID, datasetName, tableName);

        final Table table = bigquery.getTable(tableId);
        final Schema schema = table.getDefinition().getSchema();
        final FieldList schemaFields = schema == null ? null : schema.getFields();
        final String COMPRESSED_REF_RANGES_COLUMN = "packed_ref_data";

        boolean foundCompressedColumn = false;
        if (schemaFields != null) {
            for (com.google.cloud.bigquery.Field field : schemaFields) {
                if (field.getName().equals(COMPRESSED_REF_RANGES_COLUMN)) {
                    foundCompressedColumn = true;
                    break;
                }
            }
        }
        if (foundCompressedColumn ^ storeCompressedReferences) {
            throw new UserException("reference ranges table " + projectID + "." + datasetName + "." + tableName +
                    " has a schema incompatible with the invocation of this tool. " +
                    "This tool was invoked with storeCompressedReferences = " + storeCompressedReferences + ", but a check for " +
                    "compressed references column '" + COMPRESSED_REF_RANGES_COLUMN + "' returned " + foundCompressedColumn + ".");
        }
    }

    public static boolean doRowsExistFor(CommonCode.OutputType outputType, String projectId, String datasetName, String tableNumber, Long sampleId) {
        if (outputType != CommonCode.OutputType.BQ) return false;
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, REF_RANGES_FILETYPE_PREFIX + tableNumber, SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
    }

    @Override
    public void write(long location, long sampleId, int length, String state, Integer dp) throws IOException {
        write(createJsonRow(location, sampleId, length, state, dp));
    }

    @Override
    public void writeCompressed(long packedData, long sampleId, Integer dp) throws IOException {
        write(createCompressedJsonRow(packedData, sampleId, dp));
    }

    private JSONObject createJsonRow(long location, long sampleId, int length, String state, Integer dp) {
        JSONObject record = new JSONObject();
        record.put("location", location);
        record.put("sample_id", sampleId);
        record.put("length", length);
        record.put("state", state);
        if (dp != null) {
            record.put("dp", dp);
        }
        return record;
    }

    private JSONObject createCompressedJsonRow(long packedData, long sampleId, Integer dp) {
        JSONObject record = new JSONObject();
        record.put("packed_ref_data", packedData);
        record.put("sample_id", sampleId);
        if (dp != null) {
            record.put("dp", dp);
        }
        return record;
    }
}
