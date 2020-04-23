package org.broadinstitute.hellbender.utils.bigquery;

import java.util.List;

public class TableReference {
    public final String tableProject;
    public final String tableDataset;
    public final String tableName;
    public final List<String> fields;

    public TableReference(String fullyQualifiedTableName, List<String> fields) {
        String[] vals = fullyQualifiedTableName.split("\\.");
        tableProject = vals[0];
        tableDataset = vals[1];
        tableName = vals[2];
        this.fields = fields;
    }

    public TableReference(String project, String dataset, String tableName, List<String> fields) {
        this.tableProject = project;
        this.tableDataset = dataset;
        this.tableName = tableName;
        this.fields = fields;
    }

    public String getFQTableName() {
        return tableProject +"." + tableDataset + "." + tableName;
    }
}
