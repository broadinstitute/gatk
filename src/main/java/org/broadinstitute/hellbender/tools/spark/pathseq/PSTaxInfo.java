package org.broadinstitute.hellbender.tools.spark.pathseq;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Helper class for ClassifyReads that stores taxon data
 */
public final class PSTaxInfo {
    public String name = null;
    public String rank = null;
    public String parent_tax = null;
    public long length = 0;
    public final Collection<String> ref_names;

    public PSTaxInfo() {
        ref_names = new ArrayList<>();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        PSTaxInfo psTaxInfo = (PSTaxInfo) o;

        if (length != psTaxInfo.length) return false;
        if (name != null ? !name.equals(psTaxInfo.name) : psTaxInfo.name != null) return false;
        if (rank != null ? !rank.equals(psTaxInfo.rank) : psTaxInfo.rank != null) return false;
        if (parent_tax != null ? !parent_tax.equals(psTaxInfo.parent_tax) : psTaxInfo.parent_tax != null) return false;
        return ref_names.containsAll(psTaxInfo.ref_names) && ref_names.size() == psTaxInfo.ref_names.size();
    }

    @Override
    public int hashCode() {
        int result = name != null ? name.hashCode() : 0;
        result = 31 * result + (rank != null ? rank.hashCode() : 0);
        result = 31 * result + (parent_tax != null ? parent_tax.hashCode() : 0);
        result = 31 * result + (int) (length ^ (length >>> 32));
        result = 31 * result + ref_names.hashCode();
        return result;
    }

    @Override
    public String toString() {
        String str = "PSTaxInfo{" + name + "," + rank + "," + parent_tax + "," + length + ",[";
        for (final String ref : ref_names) {
            str += ref + ",";
        }
        return str + "]}";
    }
}
