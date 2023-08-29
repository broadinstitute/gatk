package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IndexedSet;

import java.util.Collection;

/**
 * Allele list implementation using an indexed-set.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IndexedAlleleList<A extends Allele> implements AlleleList<A> {

    private final IndexedSet<A> alleles;

    /**
     * Constructs a new empty allele-list
     */
    @SuppressWarnings("unchecked")
    public IndexedAlleleList() {
        alleles = new IndexedSet<>();
    }

    /**
     * Constructs a new allele-list from an array of alleles.
     *
     * <p>
     *     Repeats in the input array will be ignored (keeping the first one). The order of alleles in the
     *     resulting list is the same as in the natural traversal of the input collection.
     *
     * </p>
     * @param alleles the original allele array
     *
     * @throws IllegalArgumentException if {@code alleles} is {@code null} or contains {@code null}s.
     */
    @SafeVarargs
    @SuppressWarnings({"unchecked", "varargs"})
    public IndexedAlleleList(final A ... alleles) {
        this.alleles = new IndexedSet<>(alleles);
    }

    /**
     * Constructs a new allele-list from a collection of alleles.
     *
     * <p>
     *     Repeats in the input collection will be ignored (keeping the first one). The order of alleles in the
     *     resulting list is the same as in the natural traversal of the input collection.
     *
     * </p>
     * @param alleles the original allele collection
     *
     * @throws IllegalArgumentException if {@code alleles} is {@code null} or contains {@code null}s.
     */
    public IndexedAlleleList(final Collection<A> alleles) {
        this.alleles = new IndexedSet<>(alleles);
    }

    @Override
    public int numberOfAlleles() {
        return alleles.size();
    }

    @Override
    public int indexOfAllele(final Allele allele) {
        Utils.nonNull(allele);
        return alleles.indexOf(allele);
    }

    @Override
    public A getAllele(final int index) {
        Utils.validIndex(index, alleles.size());
        return alleles.get(index);
    }
}
