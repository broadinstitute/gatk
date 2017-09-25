package org.broadinstitute.hellbender.utils.read.markduplicates;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.ToIntFunction;


/**
 * This class helps us compute and compare duplicate scores, which are used for selecting the non-duplicate
 * during duplicate marking (see MarkDuplicates2).
 *
 * Adapted to GATK from Picard's DuplicateScoringStrategy.
 */
public enum MarkDuplicatesScoringStrategy {

    SUM_OF_BASE_QUALITIES(MarkDuplicatesScoringStrategy::sumOfBaseQualities),
    TOTAL_MAPPED_REFERENCE_LENGTH(MarkDuplicatesScoringStrategy::totalMappedReferenceLength);

    private final ToIntFunction<GATKRead> scoring;

    MarkDuplicatesScoringStrategy(final ToIntFunction<GATKRead> scoring){
        this.scoring = scoring;
    }

    /**
     * Given a read, compute the reads's score for mark duplicates.
     */
    public int score(final GATKRead read){
        Utils.nonNull(read);
        return scoring.applyAsInt(read);
    }

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    public static final int MIN_BASE_QUAL = 15;

    /**
     * Computes the sum of base qualities of the given read. 
     * Base qualities must be at least {@link MarkDuplicatesScoringStrategy#MIN_BASE_QUAL} to be counted.
     * Returns 0 if null is passed as an argument.
     */
    public static int sumOfBaseQualities(final GATKRead read) {
        if (read == null) {
            return 0;
        } else {
            int sum = 0;
            for ( final byte b : read.getBaseQualities() ) {
                int i = (int)b;
                if ( i >= MIN_BASE_QUAL ) {
                    sum += i;
                }
            }
            return sum;
        }
    }

    public static int totalMappedReferenceLength(final GATKRead read){
        return read.isUnmapped() ? 0 : read.getCigar().getReferenceLength();
    }
}
