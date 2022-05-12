package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;

import java.util.List;

/**
 * This class exists to allow VariantContext objects to be compared based only on their location and set of alleles,
 * providing a more liberal equals method so that VariantContext objects can be placed into a Set
 * which retains only VCs that have non-redundant location and Allele lists.
 */
public class LocationAndAlleles {
    private final int loc;
    private final List<Allele> alleles;

    public LocationAndAlleles(final int loc, final List<Allele> alleles) {
        this.loc = loc;
        this.alleles = alleles;
    }

    public int getLoc() {
        return loc;
    }

    public List<Allele> getAlleles() {
        return alleles;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final LocationAndAlleles that = (LocationAndAlleles) o;

        if (loc != that.loc) return false;
        return alleles != null ? alleles.equals(that.alleles) : that.alleles == null;
    }

    @Override
    public int hashCode() {
        return 31 * loc + (alleles != null ? alleles.hashCode() : 0);
    }


}