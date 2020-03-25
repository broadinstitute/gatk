package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.Comparator;

public class FRDBQDUtils {


    // TODO this needs to be evaluated as it is a hacky way of evaluating this property of the alleles to be GT called
    // currently this accounts (probably) for indels TODO (though one suspects A* alleles might cause yet more problems....)...
    public static  <A extends Allele> boolean containsInsertionOrDeletion(final AlleleList<A> allelesToTest) {
        final int refAlleleSize = allelesToTest.getAllele(allelesToTest.indexOfReference()).length();
        return allelesToTest.asListOfAlleles().stream().anyMatch(a -> !a.isSymbolic() && a.length()!=refAlleleSize);
    }

    /**
     * These two methods are taken from from the DRAGEN genotyping model for BQD
     * @param paddedReference       reference to check for homopolymer span
     * @param offsetForRefIntoEvent offset of the base upon which to make a call
     */
    public static double computeForwardHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent) {
        int length = 1;
        byte hBase = paddedReference[offsetForRefIntoEvent - length];
        while(length < 4) {
            if (hBase != paddedReference[offsetForRefIntoEvent - length]) {
                length--;
                break;
            }
            length++;
        }
        return 5.0 * length;
    }
    public static double computeReverseHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent) {
        int length = 1;
        final byte hBase = paddedReference[offsetForRefIntoEvent + length];
        while (length < 5) {
            if (hBase != paddedReference[offsetForRefIntoEvent + length]) {
                length--;
                break;
            }
            length++;
        }
        return 5.0 * length;
    }

}
