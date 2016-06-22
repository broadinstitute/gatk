package org.broadinstitute.hellbender.engine;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;

import java.io.Serializable;
import java.util.List;

/**
 * ReadContextData is additional data that's useful when processing reads. Currently the extra data is
 * ReferenceBases and Variants, either of which can be null. This is intended to be a simple data storage class.
 */
public class ReadContextData implements Serializable {
    private static final long serialVersionUID = 1L;

    private final ReferenceBases referenceBases;
    private final List<GATKVariant> variants;

    public ReadContextData( final ReferenceBases referenceBases, final Iterable<GATKVariant> variants ) {
        this.referenceBases = referenceBases;
        this.variants = Lists.newArrayList(variants);
    }

    public ReadContextData( final ReferenceBases referenceBases, final List<GATKVariant> variants ) {
        this.referenceBases = referenceBases;
        this.variants = variants;
    }

    public ReferenceBases getOverlappingReferenceBases() {
        return referenceBases;
    }

    public Iterable<GATKVariant> getOverlappingVariants() {
        return variants;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ReadContextData that = (ReadContextData) o;

        if (!referenceBases.equals(that.referenceBases)) return false;
        return variants.equals(that.variants);

    }

    @Override
    public int hashCode() {
        int result = referenceBases.hashCode();
        result = 31 * result + variants.hashCode();
        return result;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("<");
        for (GATKVariant v : variants) {
            builder.append(v).append(",");
        }
        builder.append(">");
        return "ReadContextData{" +
                "referenceBases=" + referenceBases.toString() +
                ", variants=" + builder.toString() +
                '}';
    }
}
