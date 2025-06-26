package org.broadinstitute.hellbender.tools.gvs.common;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;

/**
 * Utility class for variant sequence minimization operations.
 * Provides methods to minimize genomic sequences by removing redundant parts
 * from the right sides of both reference and allele sequences.
 */
public class VariantMinimizationUtils {
    static final Logger logger = LogManager.getLogger(VariantMinimizationUtils.class);

    /**
     * Minimizes variant genomic sequences by removing redundant parts from the right sides
     * of both the reference and allele sequences. The minimum length of either the minimized
     * allele or the reference sequence will be 1.
     *
     * @param reference the reference sequence
     * @param allele the allele sequence
     * @return a Pair where the left element is the minimized reference and the right element is the minimized allele
     * @throws UserException if either input is null or empty
     */
    public static Pair<String, String> minimize(final String reference, final String allele) {
        if (reference == null || reference.isEmpty()) {
            throw new UserException("Reference sequence cannot be null or empty");
        }
        if (allele == null || allele.isEmpty()) {
            throw new UserException("Allele sequence cannot be null or empty");
        }

        int refLength = reference.length();
        int alleleLength = allele.length();

        // Find the longest common suffix by comparing characters from the end
        int commonSuffixLength = 0;
        int maxComparable = Math.min(refLength, alleleLength);

        // Work backwards from the end, comparing characters at the same positions from the end
        for (int i = 1; i <= maxComparable; i++) {
            char refChar = reference.charAt(refLength - i);
            char alleleChar = allele.charAt(alleleLength - i);

            if (refChar == alleleChar) {
                commonSuffixLength = i;
            } else {
                break; // Stop at first mismatch
            }
        }

        // Determine how much we can actually remove while keeping at least 1 character in each
        int maxRemovableFromRef = refLength - 1;
        int maxRemovableFromAllele = alleleLength - 1;
        int actualRemoval = Math.min(commonSuffixLength, Math.min(maxRemovableFromRef, maxRemovableFromAllele));

        // Remove the determined amount from both sequences
        String minimizedRef = reference.substring(0, refLength - actualRemoval);
        String minimizedAllele = allele.substring(0, alleleLength - actualRemoval);

        return Pair.of(minimizedRef, minimizedAllele);
    }
}
