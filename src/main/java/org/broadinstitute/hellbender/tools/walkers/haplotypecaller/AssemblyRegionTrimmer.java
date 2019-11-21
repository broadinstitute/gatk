package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;

/**
 * Helper component to manage active region trimming
 *
 * <p/>
 * It receives the user arguments that controls trimming and also performs the trimming region calculation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AssemblyRegionTrimmer {

    /**
     * Holds the extension to be used
     */
    private int usableExtension;

    private ReadThreadingAssemblerArgumentCollection assemblyArgs;

    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Holds a reference the trimmer logger.
     */
    private static final Logger logger = LogManager.getLogger(AssemblyRegionTrimmer.class);

    /**
     * Initializes the trimmer.
     *
     * <p/>
     * This method should be called once and only once before any trimming is performed.
     *
     * @param assemblyArgs user arguments for the trimmer
     * @param sequenceDictionary dictionary to determine the bounds of contigs
     * @throws IllegalStateException if this trim calculator has already been initialized.
     * @throws IllegalArgumentException if the input location parser is {@code null}.
     * @throws CommandLineException.BadArgumentValue if any of the user argument values is invalid.
     */
    public AssemblyRegionTrimmer(final ReadThreadingAssemblerArgumentCollection assemblyArgs, final SAMSequenceDictionary sequenceDictionary) {
        this.assemblyArgs = Utils.nonNull(assemblyArgs);;
        this.sequenceDictionary = sequenceDictionary;
        checkUserArguments();
        usableExtension = this.assemblyArgs.extension;
    }

    /**
     * Checks user trimming argument values
     *
     * @throws CommandLineException.BadArgumentValue if there is some problem with any of the arguments values.
     */
    private void checkUserArguments() {
        if ( assemblyArgs.snpPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("paddingAroundSNPs", "" + assemblyArgs.snpPadding + "< 0");
        }
        if ( assemblyArgs.indelPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("paddingAroundIndels", "" + assemblyArgs.indelPadding + "< 0");
        }
        if ( assemblyArgs.extension < 0) {
            throw new CommandLineException.BadArgumentValue("maxDiscARExtension", "" + assemblyArgs.extension + "< 0");
        }
    }

    /**
     * Holds the result of trimming.
     */
    public static final class Result {

        /**
         * Holds the input active region.
         */
        protected final AssemblyRegion originalRegion;

        /**
         * Holds the smaller range that contain all relevant callable variants in the
         * input active region (not considering the extension).
         */
        protected final SimpleInterval variantSpan;

        /**
         * The trimmed variant region span including the extension.
         */
        protected final SimpleInterval paddedSpan;

        /**
         * Holds the flanking spans that do not contain the callable variants.
         * <p/>
         * The first element of the pair is the left (up-stream) non-variant flank, whereas the second element is
         * the right (down-stream) non-variant flank.
         */
        protected final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks;

        /**
         * Holds the collection of callable events within the variant trimming region.
         */
        protected final List<VariantContext> variants;

        /**
         * Holds variant-containing callable region.
         * <p/>
         * This is lazy-initialized using {@link #variantSpan}.
         */
        protected AssemblyRegion variantRegion;


        /**
         * Non-variant left flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getLeft()} () getFirst()}.
         */
        private AssemblyRegion leftFlankRegion;

        /**
         * Non-variant right flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getLeft()} () getSecond()}.
         */
        private AssemblyRegion rightFlankRegion;

        /**
         * Creates a trimming result given all its properties.
         * @param originalRegion the original active region.
         * @param overlappingEvents contained callable variation events.
         * @param nonVariantFlanks pair of non-variant flank spans around the variant containing span.
         * @param paddedSpan final trimmed variant span including the extension.
         * @param variantSpan variant containing span without padding.
         */
        protected Result(final AssemblyRegion originalRegion,
                         final List<VariantContext> overlappingEvents,
                         final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks,
                         final SimpleInterval paddedSpan,
                         final SimpleInterval variantSpan) {
            this.originalRegion = originalRegion;
            this.nonVariantFlanks = nonVariantFlanks;
            variants = overlappingEvents;
            this.variantSpan = variantSpan;
            this.paddedSpan = paddedSpan;

            Utils.validateArg(paddedSpan == null || variantSpan == null || paddedSpan.contains(variantSpan), "the extended callable span must include the callable span");
        }


        /**
         * Checks whether there is any variation present in the target region.
         *
         * @return {@code true} if there is any variant, {@code false} otherwise.
         */
        public boolean isVariationPresent() {
            return ! variants.isEmpty();
        }

        /**
         * Returns the trimmed variant containing region
         *
         * @throws IllegalStateException if there is no variation detected.
         *
         * @return never {@code null}.
         */
        public AssemblyRegion getVariantRegion() {
            if (variantRegion == null && paddedSpan != null) {
                variantRegion = originalRegion.trim(variantSpan, paddedSpan);
            } else if (paddedSpan == null) {
                throw new IllegalStateException("there is no variation thus no variant region");
            }
            return variantRegion;
        }

        /**
         * Checks whether there is a non-empty left flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial left flank region, {@code false} otherwise.
         */
        public boolean hasLeftFlankingRegion() {
            return nonVariantFlanks.getLeft() != null;
        }

        /**
         * Checks whether there is a non-empty right flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial right flank region, {@code false} otherwise.
         */
        public boolean hasRightFlankingRegion() {
            return nonVariantFlanks.getRight() != null;
        }

        /**
         *  Returns the trimmed out left non-variant region.
         *  <p/>
         *  Notice that in case of no variation, the whole original region is considered the left flanking region.
         *
         *  @throws IllegalStateException if there is not such as left flanking region.
         */
        public AssemblyRegion nonVariantLeftFlankRegion() {
            if (leftFlankRegion == null && nonVariantFlanks.getLeft() != null) {
                leftFlankRegion = originalRegion.trim(nonVariantFlanks.getLeft(), originalRegion.getPadding());
            } else if (nonVariantFlanks.getLeft() == null) {
                throw new IllegalStateException("there is no left flank non-variant trimmed out region");
            }
            return leftFlankRegion;
        }

        /**
         *  Returns the trimmed out right non-variant region.
         */
        public AssemblyRegion nonVariantRightFlankRegion() {
            if (rightFlankRegion == null && nonVariantFlanks.getRight() != null) {
                rightFlankRegion = originalRegion.trim(nonVariantFlanks.getRight(), originalRegion.getPadding());
            } else if (nonVariantFlanks.getRight() == null) {
                throw new IllegalStateException("there is no right flank non-variant trimmed out region");
            }
            return rightFlankRegion;
        }

        /**
         * Creates a result indicating that no variation was found.
         */
        protected static Result noVariation(final AssemblyRegion targetRegion) {
            final Result result = new Result(targetRegion, Collections.emptyList(), Pair.of(targetRegion.getSpan(), null),
                    null, null);
            result.leftFlankRegion = targetRegion;
            return result;
        }
    }

    /**
     * Returns a trimming result object from which the variant trimmed region and flanking non-variant sections
     * can be recovered latter.
     *
     * @param originalRegion the genome location range to trim.
     * @param allVariantsWithinExtendedRegion list of variants contained in the trimming location. Variants therein
     *                                        not overlapping with {@code originalRegion} are simply ignored.
     * @return never {@code null}.
     */
    public Result trim(final AssemblyRegion originalRegion,
                       final SortedSet<VariantContext> allVariantsWithinExtendedRegion) {


        if ( allVariantsWithinExtendedRegion.isEmpty() ) // no variants,
        {
            return Result.noVariation(originalRegion);
        }

        final List<VariantContext> withinActiveRegion = new LinkedList<>();
        final SimpleInterval originalRegionRange = originalRegion.getSpan();
        boolean foundNonSnp = false;
        SimpleInterval variantSpan = null;
        for ( final VariantContext vc : allVariantsWithinExtendedRegion ) {
            final SimpleInterval vcLoc = new SimpleInterval(vc);
            if ( originalRegionRange.overlaps(vcLoc) ) {
                foundNonSnp = foundNonSnp || ! vc.isSNP();
                variantSpan = variantSpan == null ? vcLoc : variantSpan.spanWith(vcLoc);
                withinActiveRegion.add(vc);
            }
        }
        final int padding = foundNonSnp ? assemblyArgs.indelPadding : assemblyArgs.snpPadding;

        // we don't actually have anything in the region after skipping out variants that don't overlap
        // the region's full location
        if ( variantSpan == null ) {
            return Result.noVariation(originalRegion);
        }

        final SimpleInterval finalSpan = originalRegionRange.expandWithinContig(usableExtension, sequenceDictionary).intersect(variantSpan.expandWithinContig(padding, sequenceDictionary)).mergeWithContiguous(variantSpan);

        // Make double sure that, if we are emitting GVCF we won't call non-variable positions beyond the target active region span.
        // In regular call we don't do so so we don't care and we want to maintain behavior, so the conditional.
        final SimpleInterval callableSpan = variantSpan.intersect(originalRegionRange);

        final Pair<SimpleInterval, SimpleInterval> nonVariantRegions = nonVariantTargetRegions(originalRegion, callableSpan);

        if ( assemblyArgs.debugAssembly ) {
            logger.info("events       : " + withinActiveRegion);
            logger.info("region       : " + originalRegion);
            logger.info("variantSpan : " + callableSpan);
            logger.info("finalSpan    : " + finalSpan);
        }

        return new Result(originalRegion, withinActiveRegion,nonVariantRegions,finalSpan, variantSpan);
    }

    /**
     * Calculates the list of region to trim away.
     * @param targetRegion region for which to generate the flanking regions.
     * @param variantSpan the span of the core region containing relevant variation and required padding.
     * @return never {@code null}; 0, 1 or 2 element list.
     */
    private Pair<SimpleInterval, SimpleInterval> nonVariantTargetRegions(final AssemblyRegion targetRegion, final SimpleInterval variantSpan) {
        final SimpleInterval targetRegionRange = targetRegion.getSpan();
        final int finalStart = variantSpan.getStart();
        final int finalStop = variantSpan.getEnd();

        final int targetStart = targetRegionRange.getStart();
        final int targetStop = targetRegionRange.getEnd();

        final boolean preTrimmingRequired = targetStart < finalStart;
        final boolean postTrimmingRequired = targetStop > finalStop;
        if (preTrimmingRequired) {
            final String contig = targetRegionRange.getContig();
            return postTrimmingRequired ? Pair.of(new SimpleInterval(contig, targetStart, finalStart - 1), new SimpleInterval(contig, finalStop + 1, targetStop)) :
                                          Pair.of(new SimpleInterval(contig, targetStart, finalStart - 1), null);
        } else if (postTrimmingRequired) {
            return Pair.of(null, new SimpleInterval(targetRegionRange.getContig(), finalStop + 1, targetStop));
        } else {
            return Pair.of(null, null);
        }
    }
}