package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * Table columns of an allelic count tab-separated file.  All columns are mandatory.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
enum AllelicCountTableColumn {
    CONTIG,
    POSITION,
    REF_COUNT,
    ALT_COUNT,
    REF_NUCLEOTIDE,
    ALT_NUCLEOTIDE;

    static final TableColumnCollection COLUMNS = new TableColumnCollection(
            CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE);
}
