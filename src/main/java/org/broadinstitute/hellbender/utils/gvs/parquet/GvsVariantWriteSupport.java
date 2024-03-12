package org.broadinstitute.hellbender.utils.gvs.parquet;

import org.apache.hadoop.conf.Configuration;
import org.apache.parquet.column.ColumnDescriptor;
import org.apache.parquet.hadoop.api.WriteSupport;
import org.apache.parquet.io.api.Binary;
import org.apache.parquet.io.api.RecordConsumer;
import org.apache.parquet.schema.MessageType;
import org.json.JSONObject;

import java.util.HashMap;
import java.util.List;

public class GvsVariantWriteSupport extends WriteSupport<JSONObject> {
    MessageType schema;
    RecordConsumer recordConsumer;
    List<ColumnDescriptor> cols;

    // support specifying encodings and compression?
    public GvsVariantWriteSupport(MessageType schema) {
        this.schema = schema;
        this.cols = schema.getColumns();
    }

    @Override
    public WriteContext init(Configuration config) {
        return new WriteContext(schema, new HashMap<>());
    }

    @Override
    public void prepareForWrite(RecordConsumer recordConsumer) {
        this.recordConsumer = recordConsumer;
    }

    /**
     * Current implementation is a FLAT one.  Will want to separate it out and structure it a little better later!
     * @param vetRow one record to write to the previously provided record consumer
     */
    @Override
    public void write(JSONObject vetRow) {
        recordConsumer.startMessage();
        // let's iterate through the possible values we have in the JSON
        for (int field = 0; field < cols.size(); ++field) {
            ColumnDescriptor col = cols.get(field);
            String columnName = col.getPrimitiveType().getName();
            // if this isn't here, we're supposed to just skip right over it
            if (vetRow.has(columnName) && vetRow.get(columnName) != JSONObject.NULL) {
                recordConsumer.startField(columnName, field);
                switch(col.getPrimitiveType().getPrimitiveTypeName()) {
                    case INT64 -> recordConsumer.addLong(vetRow.getLong(columnName));
                    case FLOAT -> recordConsumer.addFloat(vetRow.getFloat(columnName));
                    case BINARY -> recordConsumer.addBinary(Binary.fromString(vetRow.getString(columnName)));
                    default -> throw new UnsupportedOperationException("Haven't implemented other types yet! Can't process column "+columnName+ "with type "+col.getPrimitiveType().getName());
                }
                recordConsumer.endField(columnName, field);
            }
        }

        recordConsumer.endMessage();;
    }

}
