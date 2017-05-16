package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Mandatory table columns for contig germline ploidy annotation tab-separated files.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public enum ContigGermlinePloidyAnnotationTableColumn {
    /**
     * Contig identifier ("1", "chr1", "X", "chrX", etc.)
     */
    CONTIG,

    /**
     * Contig class (see {@link ContigClass})
     */
    CLASS;

    public static final TableColumnCollection MANDATORY_CONTIG_ANNOTATION_COLUMNS = new TableColumnCollection(
            CONTIG, CLASS);

    public static final Set<String> MANDATORY_CONTIG_ANNOTATION_COLUMNS_SET = Collections.unmodifiableSet(
            new HashSet<>(MANDATORY_CONTIG_ANNOTATION_COLUMNS.names()));

}
