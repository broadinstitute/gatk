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

    // Orders the reads based on the number of bases there are to the left of the fatherEndComparisonLocation as aligned according to the cigar
    public static class ReadFeatherEndForwardComparitor implements Comparator<Pair<GATKRead, Integer>>, Serializable {
        final SimpleInterval fatherEndComparisonLocation;

        public ReadFeatherEndForwardComparitor(final SimpleInterval intervalOfVariant) {
            this.fatherEndComparisonLocation = intervalOfVariant;
        }

        //TODO efficeiency

        /**
         * Evaluate first the number of bases to the end of the read and follow that by the base quality
         */
        @Override
        public int compare(final Pair<GATKRead, Integer> read1, final Pair<GATKRead, Integer> read2) {
            //TODO possibly some secondary sorting asserting
            final int read1VariantCoordinate = ReadUtils.getReadCoordinateForReferenceCoordinate(read1.getLeft(), fatherEndComparisonLocation.getStart()).getLeft();
            final int read2VariantCoordinate = ReadUtils.getReadCoordinateForReferenceCoordinate(read2.getLeft(), fatherEndComparisonLocation.getStart()).getLeft();
            int diffVal = (read1.getLeft().getLength() - read1VariantCoordinate)
                    - (read2.getLeft().getLength() - read2VariantCoordinate);
            if (diffVal == 0) {
                diffVal = read1.getLeft().getBaseQuality(read1VariantCoordinate) - read2.getLeft().getBaseQuality(read2VariantCoordinate);
            }

            return diffVal;
        }
    }

    // Orders the reads based on the number of bases in the read that occur before the fatherEndComparisonLocation as aligned according to the cigar
    public static class ReadFeatherEndRevereseComparitor implements Comparator<Pair<GATKRead, Integer>>, Serializable {
        final SimpleInterval fatherEndComparisonLocation;

        public ReadFeatherEndRevereseComparitor(final SimpleInterval intervalOfVariant) {
            this.fatherEndComparisonLocation = intervalOfVariant;
        }

        /**
         * Evaluate first the number of bases to the end of the read and follow that by the base quality
         */
        @Override
        public int compare(final Pair<GATKRead, Integer> read1, final Pair<GATKRead, Integer> read2) {
            //TODO possibly some secondary sorting asserting
            final int read1VariantCoordinate = ReadUtils.getReadCoordinateForReferenceCoordinate(read1.getLeft(), fatherEndComparisonLocation.getStart()).getLeft();
            final int read2VariantCoordinate = ReadUtils.getReadCoordinateForReferenceCoordinate(read2.getLeft(), fatherEndComparisonLocation.getStart()).getLeft();

            int diffVal = read1VariantCoordinate - read2VariantCoordinate;
            if (diffVal==0) {
                //TODO verify this lines up with the sort in DRAGBQD
                diffVal = read1.getLeft().getBaseQuality(read1VariantCoordinate) - read2.getLeft().getBaseQuality(read2VariantCoordinate);
            }
            return diffVal;
        }

    }

}
