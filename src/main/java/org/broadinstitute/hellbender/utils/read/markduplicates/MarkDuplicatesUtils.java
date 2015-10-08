package org.broadinstitute.hellbender.utils.read.markduplicates;

import org.broadinstitute.hellbender.utils.read.GATKRead;


public class MarkDuplicatesUtils {

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    public static final int MIN_BASE_QUAL = 15;

    /**
     * How to assign a score to the read in MarkDuplicates (so that we pick the best one to be the non-duplicate).
     */
    //Note: copied from htsjdk.samtools.DuplicateScoringStrategy
    public static int scoreForRead(final GATKRead record) {
        if (record == null) {
            return 0;
        } else {
            int sum = 0;
            for ( byte b : record.getBaseQualities() ) {
                int i = (int)b;
                if ( i >= MIN_BASE_QUAL ) {
                    sum += i;
                }
            }
            return sum;
        }
    }
}
