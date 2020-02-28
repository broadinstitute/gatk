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
        byte hBase = paddedReference[offsetForRefIntoEvent];
        int length = 1;
        while(length < 5) {
            if (hBase != paddedReference[offsetForRefIntoEvent + length]) {
                length--;
                break;
            }
            length++;
        }
        return 5.0 * length;
    }
    public static double computeReverseHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent) {
        final byte hBase = paddedReference[offsetForRefIntoEvent];
        int length = 1;
        while(length < 5) {
            if (hBase != paddedReference[offsetForRefIntoEvent - length]) {
                length--;
                break;
            }
            length++;
        }
        return 5.0 * length;
    }

    // Orders the reads based on the number of bases there are to the left of the fatherEndComparisonLocation as aligned according to the cigar
    public static class ReadFeatherEndForwardComparitor implements Comparator<Pair<Pair<GATKRead,Integer>, Integer>>, Serializable {
        //TODO efficeiency

        /**
         * Evaluate first the number of bases to the end of the read and follow that by the base quality
         */
        @Override
        public int compare(final Pair<Pair<GATKRead,Integer>, Integer> read1, final Pair<Pair<GATKRead,Integer>, Integer> read2) {
            //TODO possibly some secondary sorting asserting
            final int read1VariantCoordinate = read1.getLeft().getRight();
            final int read2VariantCoordinate = read2.getLeft().getRight();
            int diffVal = (read1.getLeft().getLeft().getLength() - read1VariantCoordinate)
                    - (read2.getLeft().getLeft().getLength() - read2VariantCoordinate);
            if (diffVal == 0) {
                diffVal = (read1VariantCoordinate != -1 ? read1.getLeft().getLeft().getBaseQuality(read1VariantCoordinate) : 0)
                        - (read2VariantCoordinate != -1 ? read2.getLeft().getLeft().getBaseQuality(read2VariantCoordinate) : 0);
            }

            return diffVal;
        }
    }

    // Orders the reads based on the number of bases in the read that occur before the fatherEndComparisonLocation as aligned according to the cigar
    public static class ReadFeatherEndRevereseComparitor implements Comparator<Pair<Pair<GATKRead,Integer>, Integer>>, Serializable {

        /**
         * Evaluate first the number of bases to the end of the read and follow that by the base quality
         */
        @Override
        public int compare(final Pair<Pair<GATKRead,Integer>, Integer> read1, final Pair<Pair<GATKRead,Integer>, Integer> read2) {
            final int read1VariantCoordinate = read1.getLeft().getRight();
            final int read2VariantCoordinate = read2.getLeft().getRight();

            int diffVal = read1VariantCoordinate - read2VariantCoordinate;
            if (diffVal==0) {
                //TODO verify this lines up with the sort in DRAGBQD
                diffVal = (read1VariantCoordinate != -1 ? read1.getLeft().getLeft().getBaseQuality(read1VariantCoordinate) : 0)
                        - (read2VariantCoordinate != -1 ? read2.getLeft().getLeft().getBaseQuality(read2VariantCoordinate) : 0);
            }
            return diffVal;
        }

    }

}
