package org.broadinstitute.hellbender.utils.tsv;

import java.util.HashMap;
import java.util.Map;

public class TableReaderOptions {
    boolean headerIsLastComment = false;

    Map<String, String> columnRenamer = new HashMap<>();

    public TableReaderOptions() {}

    public TableReaderOptions(final boolean headerIsLastComment, final Map<String, String> columnRenamer) {
        this.headerIsLastComment = headerIsLastComment;
        this.columnRenamer = columnRenamer;
    }
}
