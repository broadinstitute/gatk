package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Helper class for ClassifyReads that stores the name, taxonomic class and parent, reference length,
 * and reference contig names of a given taxon in the pathogen reference.
 */
public final class PSPathogenReferenceTaxonProperties {
    public String name = null;
    public String rank = null;
    public String parentTaxId = null;
    private long length = 0;
    private final Map<String, Long> accessions;

    public PSPathogenReferenceTaxonProperties() {
        accessions = new HashMap<>(SVUtils.hashMapCapacity(1));
    }

    public void addAccession(final String name, final long length) {
        accessions.put(name, length);
        this.length += length;
    }

    public boolean hasAccession(final String name) {
        return accessions.containsKey(name);
    }

    public Set<String> getAccessions() {
        return accessions.keySet();
    }

    public long getAccessionLength(final String name) {
        if (accessions.containsKey(name)) {
            return accessions.get(name);
        }
        return 0;
    }

    public long getTotalLength() {
        return length;
    }
}
