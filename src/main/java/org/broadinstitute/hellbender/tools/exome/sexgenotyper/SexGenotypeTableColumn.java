package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Mandatory columns for a sample sex genotype table.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public enum SexGenotypeTableColumn {

    /**
     * Sample name; if the table is used as an input for the target coverage modeller, sample names must match
     * those in {@link org.broadinstitute.hellbender.tools.exome.ReadCountCollection}
     */
    SAMPLE_NAME,

    /**
     * Sex genotype; if the table is used as an input for target coverage modeller, genotype names must match
     * those in {@link GermlinePloidyAnnotatedTargetCollection}
     */
    SEX_GENOTYPE;

    public static final TableColumnCollection MANDATORY_SEX_GENOTYPE_COLUMNS = new TableColumnCollection(
            SAMPLE_NAME, SEX_GENOTYPE);

    public static final Set<String> MANDATORY_SEX_GENOTYPE_COLUMNS_SET = Collections.unmodifiableSet(
            new HashSet<>(MANDATORY_SEX_GENOTYPE_COLUMNS.names()));

}
