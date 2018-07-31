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
    private String name;
    private String rank = null;
    private int parentTaxId = PSTree.NULL_NODE;
    private long length = 0;
    private final Map<String, Long> accessions;

    public PSPathogenReferenceTaxonProperties() {
        accessions = new HashMap<>(SVUtils.hashMapCapacity(1));
    }

    public PSPathogenReferenceTaxonProperties(final String name) {
        this();
        this.name = name;
    }

    public void addAccession(final String name, final long length) {
        accessions.put(name, length);
        this.length += length;
    }

    public boolean hasAccession(final String name) {
        return accessions.containsKey(name);
    }

    public long getAccessionLength(final String name) {
        if (accessions.containsKey(name)) {
            return accessions.get(name);
        }
        return 0;
    }

    public void setName(final String name) { this.name = name;}
    public void setRank(final String rank) { this.rank = rank;}
    public void setParent(final int parentTaxId) { this.parentTaxId = parentTaxId;}

    public String getName() { return name; }
    public String getRank() { return rank; }
    public int getParent() { return parentTaxId; }
    public Set<String> getAccessions() {
        return accessions.keySet();
    }
    public long getTotalLength() {
        return length;
    }
}
