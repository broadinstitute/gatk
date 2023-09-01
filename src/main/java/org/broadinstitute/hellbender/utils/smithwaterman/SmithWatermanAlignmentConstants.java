package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

/**
 * This class collects the various {@link SWParameters} that are used for various alignment procedures.
 * It is likely that these parameters have not been rigorously optimized.
 * Documentation is also somewhat lacking, but we preserve any original comments that may have accompanied each set.
 * See also some relevant issues and comments:
 *  <p>
 *      <a href="http://github.com/broadinstitute/gatk/issues/2498">http://github.com/broadinstitute/gatk/issues/2498</a>
 *      <a href="http://github.com/broadinstitute/gatk/issues/5564">http://github.com/broadinstitute/gatk/issues/5564</a>
 *      <a href="http://github.com/broadinstitute/gatk/pull/4858#discussion_r194048530">http://github.com/broadinstitute/gatk/pull/4858#discussion_r194048530</a>
 *  </p>
 */
public final class SmithWatermanAlignmentConstants {
    private SmithWatermanAlignmentConstants() {}

    /**
     * {@code ORIGINAL_DEFAULT} is only used in test code. It is worth noting that these tests are somewhat insensitive
     * to the particular values used; (e.g., the majority pass if {@link SmithWatermanAlignmentConstants#STANDARD_NGS}
     * is substituted, and all pass if {@link SmithWatermanAlignmentConstants#NEW_SW_PARAMETERS} is substituted).
     *
     * Original comments:
     *      match=1, mismatch = -1/3, gap=-(1+k/3)
     */
    public static final SWParameters ORIGINAL_DEFAULT = new SWParameters(3, -1, -4, -3);

    /**
     * {@code STANDARD_NGS} is the default for {@link AbstractReadThreadingGraph} methods for the recovery of dangling heads/tails.
     *
     * Original comments:
     *      none
     */
    public static final SWParameters STANDARD_NGS = new SWParameters(25, -50, -110, -6);

    /**
     * {@code NEW_SW_PARAMETERS} is the default for {@link CigarUtils#calculateCigar} for haplotype-to-reference alignment.
     * It was added in <a href="https://github.com/broadinstitute/gatk/pull/586">https://github.com/broadinstitute/gatk/pull/586</a>
     * at the same time as the {@link CigarUtils#calculateCigar} method. The original comments indicate that these values
     * were chosen via an optimization procedure, but no record or documentation of this procedure exists.
     * As indicated below in the original comments of {@link SmithWatermanAlignmentConstants#ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS},
     * the values chosen here heavily favor indels over substitutions.
     *
     * Original comments:
     *      used in the bubble state machine to apply Smith-Waterman to the bubble sequence
     *      these values were chosen via optimization against the NA12878 knowledge base
     */
    public static final SWParameters NEW_SW_PARAMETERS = new SWParameters(200, -150, -260, -11);

    /**
     * {@code ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS} is the default for read-to-haplotype alignment and was added in
     * <a href="https://github.com/broadinstitute/gatk/pull/4858">https://github.com/broadinstitute/gatk/pull/4858</a>
     * (superseding the use of {@link SmithWatermanAlignmentConstants#NEW_SW_PARAMETERS} in {@link AlignmentUtils#createReadAlignedToRef} in the
     * read-to-haplotype alignment step of that method).
     *
     * Original comments:
     *      In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
     *      The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
     *      Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
     *      would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
     *      Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
     *      we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
     */
    public static final SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS = new SWParameters(10, -15, -30, -5);
}
