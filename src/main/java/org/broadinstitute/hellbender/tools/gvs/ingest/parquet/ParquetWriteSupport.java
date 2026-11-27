package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.parquet.hadoop.api.WriteSupport;
import org.apache.parquet.io.api.Binary;
import org.apache.parquet.io.api.RecordConsumer;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.Type;
import org.jetbrains.annotations.NotNull;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.HashMap;

/**
 * A generic WriteSupport implementation for writing JSONObject records to Parquet files
 * in the GVS system.
 *
 * <p>Supports primitive fields only, required, optional or repeated: INT64, FLOAT, BOOLEAN and
 * BINARY (UTF8 string). A repeated primitive is the 1-level list encoding, as in
 * {@code repeated binary vrs_allele_ids (UTF8)}, which BigQuery loads as a REPEATED column.
 * Group fields, including the 3-level LIST encoding, are rejected.
 */
public class ParquetWriteSupport extends WriteSupport<JSONObject> {
    MessageType schema;
    RecordConsumer recordConsumer;

    // support specifying encodings and compression?
    public ParquetWriteSupport(@NotNull MessageType schema) {
        this.schema = schema;
    }

    /**
     * Initializes the WriteSupport with schema information.
     *
     * Note: This method signature init(Configuration) is deprecated in Parquet 1.13.1,
     * but it's still the required method to implement WriteSupport. The Configuration
     * parameter is not used in our implementation as the schema is provided via constructor.
     * We suppress the deprecation warning because this is currently the only way to
     * implement WriteSupport in this version of Parquet.
     *
     * @param config Configuration object (unused in our implementation)
     * @return WriteContext containing the schema and metadata
     */
    @Deprecated
    @Override
    public WriteContext init(Configuration config) {
        return new WriteContext(schema, new HashMap<>());
    }

    @Override
    public void prepareForWrite(RecordConsumer recordConsumer) {
        this.recordConsumer = recordConsumer;
    }

    /**
     * Writes one record. Iterates top-level schema fields; handles scalar and repeated primitives.
     *
     * @param record one record to write to the previously provided record consumer
     */
    @Override
    public void write(JSONObject record) {
        recordConsumer.startMessage();

        for (int fieldIndex = 0; fieldIndex < schema.getFieldCount(); fieldIndex++) {
            Type fieldType = schema.getType(fieldIndex);
            String fieldName = fieldType.getName();

            if (!record.has(fieldName) || record.get(fieldName) == JSONObject.NULL) {
                continue;
            }

            if (fieldType.isPrimitive()) {
                // Primitive field: support scalar and repeated values.
                if (fieldType.isRepetition(Type.Repetition.REPEATED)) {
                    JSONArray jsonArray = record.getJSONArray(fieldName);
                    if (jsonArray.isEmpty()) {
                        continue;
                    }

                    recordConsumer.startField(fieldName, fieldIndex);
                    for (int i = 0; i < jsonArray.length(); i++) {
                        switch (fieldType.asPrimitiveType().getPrimitiveTypeName()) {
                            case INT64  -> recordConsumer.addLong(jsonArray.getLong(i));
                            case FLOAT  -> recordConsumer.addFloat((float) jsonArray.getDouble(i));
                            case BOOLEAN -> recordConsumer.addBoolean(jsonArray.getBoolean(i));
                            case BINARY -> recordConsumer.addBinary(Binary.fromString(jsonArray.getString(i)));
                            default     -> throw new UnsupportedOperationException(
                                    "Unsupported primitive type for column '" + fieldName + "': "
                                    + fieldType.asPrimitiveType().getPrimitiveTypeName());
                        }
                    }
                    recordConsumer.endField(fieldName, fieldIndex);
                } else {
                    // Scalar primitive field
                    recordConsumer.startField(fieldName, fieldIndex);
                    switch (fieldType.asPrimitiveType().getPrimitiveTypeName()) {
                        case INT64  -> recordConsumer.addLong(record.getLong(fieldName));
                        case FLOAT  -> recordConsumer.addFloat(record.getFloat(fieldName));
                        case BOOLEAN -> recordConsumer.addBoolean(record.getBoolean(fieldName));
                        case BINARY -> recordConsumer.addBinary(Binary.fromString(record.getString(fieldName)));
                        default     -> throw new UnsupportedOperationException(
                                "Unsupported primitive type for column '" + fieldName + "': "
                                + fieldType.asPrimitiveType().getPrimitiveTypeName());
                    }
                    recordConsumer.endField(fieldName, fieldIndex);
                }
            } else {
                throw new UnsupportedOperationException("Group fields are not supported; column '" + fieldName
                        + "' is a group. Declare a list as a repeated primitive instead.");
            }
        }

        recordConsumer.endMessage();
    }
}
