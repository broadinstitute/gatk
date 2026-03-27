package org.broadinstitute.hellbender.tools.gvs.ingest.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.parquet.hadoop.api.WriteSupport;
import org.apache.parquet.io.api.Binary;
import org.apache.parquet.io.api.RecordConsumer;
import org.apache.parquet.schema.GroupType;
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
 * <p>Supports two field shapes:
 * <ul>
 *   <li>Scalar primitives: INT64, FLOAT, BINARY (UTF8 string)</li>
 *   <li>Repeated string lists using the standard Parquet 3-level list encoding:
 *     <pre>
 *     optional group &lt;name&gt; (LIST) {
 *       repeated group list {
 *         optional binary element (UTF8);
 *       }
 *     }
 *     </pre>
 *   </li>
 * </ul>
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
     * Writes one record. Iterates top-level schema fields; handles scalar primitives
     * and repeated string lists (3-level Parquet list encoding).
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
                // Scalar primitive field
                recordConsumer.startField(fieldName, fieldIndex);
                switch (fieldType.asPrimitiveType().getPrimitiveTypeName()) {
                    case INT64  -> recordConsumer.addLong(record.getLong(fieldName));
                    case FLOAT  -> recordConsumer.addFloat(record.getFloat(fieldName));
                    case BINARY -> recordConsumer.addBinary(Binary.fromString(record.getString(fieldName)));
                    default     -> throw new UnsupportedOperationException(
                            "Unsupported primitive type for column '" + fieldName + "': "
                            + fieldType.asPrimitiveType().getPrimitiveTypeName());
                }
                recordConsumer.endField(fieldName, fieldIndex);
            } else {
                // Group field — only repeated string lists are supported here.
                // Standard Parquet 3-level list encoding:
                //   optional group <name> (LIST) {
                //     repeated group list {
                //       optional binary element (UTF8);
                //     }
                //   }
                JSONArray jsonArray = record.getJSONArray(fieldName);
                if (jsonArray.isEmpty()) {
                    continue; // skip empty lists (same as null)
                }

                GroupType listGroup = fieldType.asGroupType();
                // level 1: the outer LIST group
                recordConsumer.startField(fieldName, fieldIndex);
                recordConsumer.startGroup();

                // level 2: the repeated "list" group
                GroupType repeatedGroup = listGroup.getType(0).asGroupType();
                recordConsumer.startField(repeatedGroup.getName(), 0);

                for (int i = 0; i < jsonArray.length(); i++) {
                    recordConsumer.startGroup();
                    // level 3: the "element" primitive
                    recordConsumer.startField(repeatedGroup.getType(0).getName(), 0);
                    recordConsumer.addBinary(Binary.fromString(jsonArray.getString(i)));
                    recordConsumer.endField(repeatedGroup.getType(0).getName(), 0);
                    recordConsumer.endGroup();
                }

                recordConsumer.endField(repeatedGroup.getName(), 0);
                recordConsumer.endGroup();
                recordConsumer.endField(fieldName, fieldIndex);
            }
        }

        recordConsumer.endMessage();
    }
}
