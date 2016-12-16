package org.broadinstitute.hellbender.tools.spark.pathseq;

import java.util.Map;

/**
 * Helper class for holding taxonomy data used by ClassifyReads
 */
public class PSTaxonomyDatabase {
    public final PSTree tree;
    public final Map<String, String> contigToTaxIDMap; //Reference contig name to taxonomic ID

    public PSTaxonomyDatabase(final PSTree tree, final Map<String, String> map) {
        this.tree = tree;
        this.contigToTaxIDMap = map;
    }
}
