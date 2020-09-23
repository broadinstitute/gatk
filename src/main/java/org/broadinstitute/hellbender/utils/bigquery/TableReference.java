package org.broadinstitute.hellbender.utils.bigquery;

import com.google.common.collect.ImmutableList;

import java.util.List;

/**
 * A reference to a BigQuery table by project, dataset, and table name, along with the contained fields.
 */
public class TableReference {
    public final String tableProject;
    public final String tableDataset;
    public final String tableName;
    public final ImmutableList<String> fields;

    /**
     * TableReference from fully qualified table name.
     * @param fullyQualifiedTableName Formatted as PROJECT.DATASET.TABLENAME, separated by "."
     * @param fields Fields contained in the table.
     */
    public TableReference(String fullyQualifiedTableName, List<String> fields) {
        String[] vals = fullyQualifiedTableName.split("\\.");
        if (vals.length != 3) {
            throw new IllegalArgumentException(String.format("Fully Qualified Table should be PROJECT.DATASET.TABLENAME, but only %s values were found in %s.", vals.length, fullyQualifiedTableName));
        }
        if ( vals[0].equals("") || vals[1].equals("") || vals[2].equals("") ) {
            throw new IllegalArgumentException(String.format("Fully Qualified Table should be PROJECT.DATASET.TABLENAME, but some values not found in %s.", fullyQualifiedTableName));
        }
        tableProject = vals[0];
        tableDataset = vals[1];
        tableName = vals[2];
        this.fields = ImmutableList.copyOf(fields);
    }

    public TableReference(String project, String dataset, String tableName, List<String> fields) {
        this.tableProject = project;
        this.tableDataset = dataset;
        this.tableName = tableName;
        this.fields = ImmutableList.copyOf(fields);
    }

    public String getFQTableName() {
        return tableProject +"." + tableDataset + "." + tableName;
    }
}
