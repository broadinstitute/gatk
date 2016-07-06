package org.broadinstitute.hellbender.tools.pon.allelic;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * Table columns of an {@link AllelicPanelOfNormals} tab-separated file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public enum AllelicPanelOfNormalsTableColumn {
    CONTIG, POSITION, ALPHA, BETA;

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
}
