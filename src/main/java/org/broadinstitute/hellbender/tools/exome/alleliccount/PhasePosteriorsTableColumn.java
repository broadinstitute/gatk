package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * Table columns for phase posteriors.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum PhasePosteriorsTableColumn {
    REF_MINOR_PROB, ALT_MINOR_PROB, OUTLIER_PROB;

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());

    public static TableColumnCollection appendPhasePosteriorColumns(final TableColumnCollection allelicCountTableColumns) {
        return new TableColumnCollection(ListUtils.union(allelicCountTableColumns.names(), COLUMNS.names()));
    }
}
