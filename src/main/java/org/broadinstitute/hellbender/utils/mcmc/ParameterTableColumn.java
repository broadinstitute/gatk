package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum ParameterTableColumn {
    PARAMETER_NAME,
    POSTERIOR_MODE,
    POSTERIOR_LOWER,
    POSTERIOR_UPPER,
    POSTERIOR_10,
    POSTERIOR_20,
    POSTERIOR_30,
    POSTERIOR_40,
    POSTERIOR_50,
    POSTERIOR_60,
    POSTERIOR_70,
    POSTERIOR_80,
    POSTERIOR_90;


    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
}
