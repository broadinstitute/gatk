package org.broadinstitute.hellbender.utils.bigquery;

import com.google.cloud.bigquery.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Iterator;

public class QueryAPIRowReader implements Iterable<FieldValueList>, Iterator<FieldValueList> {
    private static final Logger logger = LogManager.getLogger(QueryAPIRowReader.class);

    private Iterator<FieldValueList> rowIterator;

    public QueryAPIRowReader(TableResult tableResult) {
        rowIterator = tableResult.iterateAll().iterator();

    }
    @Override
    public Iterator<FieldValueList> iterator() {
        return rowIterator;
    }

    @Override
    public boolean hasNext() {
        return rowIterator.hasNext();
    }

    @Override
    public FieldValueList next() {
        return rowIterator.next();
    }

}
