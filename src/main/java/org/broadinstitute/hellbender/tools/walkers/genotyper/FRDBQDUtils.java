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

        @Override
        public int compare(final Pair<GATKRead, Integer> read1, final Pair<GATKRead, Integer> read2) {
            //TODO possibly some secondary sorting asserting
            return (read1.getLeft().getLength() - ReadUtils.getReadCoordinateForReferenceCoordinate(read1.getLeft(), fatherEndComparisonLocation.getStart()).getLeft())
                    - (read2.getLeft().getLength() - ReadUtils.getReadCoordinateForReferenceCoordinate(read2.getLeft(), fatherEndComparisonLocation.getStart()).getLeft()) ;
        }
    }

    // Orders the reads based on the number of bases in the read that occur before the fatherEndComparisonLocation as aligned according to the cigar
    public static class ReadFeatherEndRevereseComparitor implements Comparator<Pair<GATKRead, Integer>>, Serializable {
        final SimpleInterval fatherEndComparisonLocation;

        public ReadFeatherEndRevereseComparitor(final SimpleInterval intervalOfVariant) {
            this.fatherEndComparisonLocation = intervalOfVariant;
        }

        @Override
        public int compare(final Pair<GATKRead, Integer> read1, final Pair<GATKRead, Integer> read2) {
            //TODO possibly some secondary sorting asserting
            return ReadUtils.getReadCoordinateForReferenceCoordinate(read1.getLeft(), fatherEndComparisonLocation.getStart()).getLeft()
                    - ReadUtils.getReadCoordinateForReferenceCoordinate(read2.getLeft(), fatherEndComparisonLocation.getStart()).getLeft() ;
        }

    }

}
