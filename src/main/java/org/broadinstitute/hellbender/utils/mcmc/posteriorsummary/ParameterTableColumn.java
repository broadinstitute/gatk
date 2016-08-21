package org.broadinstitute.hellbender.utils.mcmc.posteriorsummary;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum ParameterTableColumn {
    PARAMETER_NAME("Parameter"),
    PARAMETER_POSTERIOR_MODE("Post_Mode"),
    PARAMETER_POSTERIOR_LOWER("Post_Lo"),
    PARAMETER_POSTERIOR_UPPER("Post_Hi"),
    PARAMETER_POSTERIOR_0("Post_0"),
    PARAMETER_POSTERIOR_10("Post_10"),
    PARAMETER_POSTERIOR_20("Post_20"),
    PARAMETER_POSTERIOR_30("Post_30"),
    PARAMETER_POSTERIOR_40("Post_40"),
    PARAMETER_POSTERIOR_50("Post_50"),
    PARAMETER_POSTERIOR_60("Post_60"),
    PARAMETER_POSTERIOR_70("Post_70"),
    PARAMETER_POSTERIOR_80("Post_80"),
    PARAMETER_POSTERIOR_90("Post_90"),
    PARAMETER_POSTERIOR_100("Post_100");

    private final String columnName;  //store the column names

    ParameterTableColumn(final String columnName) {  this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
}
