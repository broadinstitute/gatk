package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;


public enum OrientationBiasFilterSummaryTableColumn {

    SAMPLE("Sample"), ARTIFACT_MODE("Artifact_Mode"),  ARTIFACT_MODE_COMPLEMENT("Artifact_Mode_Complement"), OBQ(OrientationBiasFilterConstants.PRE_ADAPTER_METRIC_FIELD_NAME),
    NUM_OB("Num_Artifact_Mode"), NUM_FILTERED("Num_Artifact_Mode_Filtered"), NUM_NOT_AM("Num_Not_Artifact_Mode"), NUM_VARIANTS("Num_NonRef_Passing_Variants");

    private final String columnName;

    OrientationBiasFilterSummaryTableColumn(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() { return columnName; }

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
}
