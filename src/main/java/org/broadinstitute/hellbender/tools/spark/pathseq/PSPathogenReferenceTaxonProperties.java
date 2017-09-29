package org.broadinstitute.hellbender.tools.spark.pathseq;

import java.util.HashSet;
import java.util.Set;

/**
 * Helper class for ClassifyReads that stores the name, taxonomic class and parent, reference length,
 * and reference contig names of a given taxon in the pathogen reference.
 */
public final class PSPathogenReferenceTaxonProperties {
    public String name = null;
    public String rank = null;
    public String parentTaxId = null;
    public long length = 0;
    public final Set<String> accessions;

    public PSPathogenReferenceTaxonProperties() {
        accessions = new HashSet<>(1);
    }
}
