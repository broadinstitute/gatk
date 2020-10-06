package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.FieldValueList;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;

/** A BigQuery record for a row of results from a query. */
public class QueryRecord implements GenericRecord {
    private final FieldValueList fields;

    public QueryRecord(FieldValueList fields) {
        this.fields = fields;
    }
    @Override
    public void put(String key, Object v) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Object get(String key) {
        return fields.get(key).getStringValue();
    }

    @Override
    public void put(int i, Object v) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Object get(int i) {
        return null;
    }

    @Override
    public Schema getSchema() {
        return null;
    }
}
