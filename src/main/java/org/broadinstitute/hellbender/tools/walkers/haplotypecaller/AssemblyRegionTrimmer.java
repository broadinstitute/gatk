package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.TandemRepeat;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Event;

import java.util.LinkedList;
import java.util.List;
import java.util.Optional;
import java.util.SortedSet;
import java.util.stream.Collectors;

/**
 * Helper component to manage active region trimming
 *
 * <p/>
 * It receives the user arguments that controls trimming and also performs the trimming region calculation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AssemblyRegionTrimmer {

    private AssemblyRegionArgumentCollection assemblyRegionArgs;

    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Initializes the trimmer.
     *
     * <p/>
     * This method should be called once and only once before any trimming is performed.
     *
     * @param assemblyRegionArgs user arguments for the trimmer
     * @param sequenceDictionary dictionary to determine the bounds of contigs
     * @throws IllegalStateException if this trim calculator has already been initialized.
     * @throws IllegalArgumentException if the input location parser is {@code null}.
     * @throws CommandLineException.BadArgumentValue if any of the user argument values is invalid.
     */
    public AssemblyRegionTrimmer(final AssemblyRegionArgumentCollection assemblyRegionArgs, final SAMSequenceDictionary sequenceDictionary) {
        this.assemblyRegionArgs = Utils.nonNull(assemblyRegionArgs);
        this.sequenceDictionary = sequenceDictionary;
        assemblyRegionArgs.validate();
    }

    /**
     * Holds the result of trimming.
     */
    public final class Result {

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
         * Creates a trimming result given all its properties.
         * @param originalRegion the original active region.
         * @param variantSpan variant containing span without padding.
         * @param paddedSpan final trimmed variant span including the extension.
         */
        protected Result(final AssemblyRegion originalRegion, final SimpleInterval variantSpan, final SimpleInterval paddedSpan) {
            this.originalRegion = originalRegion;
            this.variantSpan = variantSpan;
            this.paddedSpan = paddedSpan;

            Utils.validateArg(paddedSpan == null || variantSpan == null || paddedSpan.contains(variantSpan), "the padded span must include the variant span");
        }

        /**
         * Checks whether there is any variation present in the target region.
         *
         * @return {@code true} if there is any variant, {@code false} otherwise.
         */
        public boolean isVariationPresent() {
            return variantSpan != null;
        }

        /**
         * Returns the trimmed variant containing region
         *
         * @throws IllegalStateException if there is no variation detected.
         *
         * @return never {@code null}.
         */
        public AssemblyRegion getVariantRegion() {
            Utils.validate(isVariationPresent(), "There is no variation present.");
            return originalRegion.trim(variantSpan, paddedSpan);
        }

        /**
         * Returns the trimmed variant containing a specific region.
         *
         * @return never {@code null}.
         */
        public AssemblyRegion getVariantRegion(SimpleInterval variantSpan, SimpleInterval paddedSpan) {

            return originalRegion.trim(variantSpan, paddedSpan);
        }

        /**
         *  Returns the trimmed out left non-variant region.
         *  <p/>
         *  Notice that in case of no variation, the whole original region is considered the left flanking region.
         *
         *  @throws IllegalStateException if there is not such as left flanking region.
         */
        public Optional<AssemblyRegion> nonVariantLeftFlankRegion() {
            if (!isVariationPresent()) {
                return Optional.of(originalRegion);
            } else if (originalRegion.getStart() < variantSpan.getStart()) {
                final SimpleInterval leftFlank = new SimpleInterval(originalRegion.getContig(), originalRegion.getStart(), variantSpan.getStart() - 1);
                return Optional.of(originalRegion.trim(leftFlank, assemblyRegionArgs.assemblyRegionPadding));
            } else {
                return Optional.empty();
            }
        }

        /**
         *  Returns the trimmed out right non-variant region.
         */
        public Optional<AssemblyRegion> nonVariantRightFlankRegion() {
            if (variantSpan.getEnd() < originalRegion.getEnd()) {
                final SimpleInterval rightFlank = new SimpleInterval(originalRegion.getContig(), variantSpan.getEnd() + 1, originalRegion.getEnd());
                return Optional.of(originalRegion.trim(rightFlank, assemblyRegionArgs.assemblyRegionPadding));
            } else {
                return Optional.empty();
            }
        }
    }

    /**
     * Creates a result indicating that no variation was found.
     */
    protected Result noVariation(final AssemblyRegion targetRegion) {
        return new Result(targetRegion, null, null);
    }

    /**
     * Returns a trimming result object from which the variant trimmed region and flanking non-variant sections
     * can be recovered latter.
     *
     * @param region the genome location range to trim.
     * @param events list of variants contained in the trimming location. Variants therein
     *                                        not overlapping with {@code region} are simply ignored.
     * @param referenceContext
     * @return never {@code null}.
     */
    public Result trim(final AssemblyRegion region, final SortedSet<Event> events, ReferenceContext referenceContext) {
        if (assemblyRegionArgs.enableLegacyAssemblyRegionTrimming) {
            return trimLegacy(region, events);
        }

        final List<Event> eventsInRegion = events.stream().filter(region::overlaps).collect(Collectors.toList());

        if ( eventsInRegion.isEmpty() ) {
            return noVariation(region);
        }

        int minStart = eventsInRegion.stream().mapToInt(Event::getStart).min().getAsInt();
        int maxEnd = eventsInRegion.stream().mapToInt(Event::getEnd).max().getAsInt();
        final SimpleInterval variantSpan = new SimpleInterval(region.getContig(), minStart, maxEnd).intersect(region);

        for (final Event event : eventsInRegion) {
            int padding = assemblyRegionArgs.snpPaddingForGenotyping;
            if (event.isIndel()) {
                padding = assemblyRegionArgs.indelPaddingForGenotyping;
                final Pair<List<Integer>, byte[]> numRepeatsAndUnit = TandemRepeat.getNumTandemRepeatUnits(referenceContext, event);
                if (numRepeatsAndUnit != null && numRepeatsAndUnit.getRight() != null) {
                    final int repeatLength = numRepeatsAndUnit.getRight() == null ? 0 : numRepeatsAndUnit.getRight().length;
                    final int mostRepeats = numRepeatsAndUnit.getLeft().stream().max(Integer::compareTo).orElse(0);
                    final int longestSTR = mostRepeats * repeatLength;
                    padding = assemblyRegionArgs.strPaddingForGenotyping + longestSTR;
                }
            }

            minStart = Math.min(minStart, Math.max(event.getStart() - padding,1));
            maxEnd = Math.max(maxEnd, event.getEnd() + padding);
        }

        final SimpleInterval paddedVariantSpan = new SimpleInterval(region.getContig(), minStart, maxEnd).intersect(region.getPaddedSpan());

        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("Padded and trimmed the region to this span: "+ paddedVariantSpan);
        }
        return new Result(region, variantSpan, paddedVariantSpan);
    }


    /**
     * Returns a trimming result object from which the variant trimmed region and flanking non-variant sections
     * can be recovered latter.
     *
     * @param originalRegion the genome location range to trim.
     * @param allEventsInExtendedRegion list of variants contained in the trimming location. Variants therein
     *                                        not overlapping with {@code originalRegion} are simply ignored.
     * @return never {@code null}.
     */
    public Result trimLegacy(final AssemblyRegion originalRegion, final SortedSet<Event> allEventsInExtendedRegion) {

        if ( allEventsInExtendedRegion.isEmpty() ) {
            return noVariation(originalRegion);
        }

        final List<Event> withinActiveRegion = new LinkedList<>();
        final SimpleInterval originalRegionRange = originalRegion.getSpan();
        boolean foundNonSnp = false;
        SimpleInterval variantSpan = null;
        for ( final Event event : allEventsInExtendedRegion ) {
            final SimpleInterval vcLoc = new SimpleInterval(event);
            if ( originalRegionRange.overlaps(vcLoc) ) {
                foundNonSnp = foundNonSnp || ! event.isSNP();
                variantSpan = variantSpan == null ? vcLoc : variantSpan.spanWith(vcLoc);
                withinActiveRegion.add(event);
            }
        }
        final int padding = foundNonSnp ? assemblyRegionArgs.indelPaddingForGenotyping : assemblyRegionArgs.snpPaddingForGenotyping;

        // we don't actually have anything in the region after skipping out variants that don't overlap
        // the region's full location
        if ( variantSpan == null ) {
            return noVariation(originalRegion);
        }

        final SimpleInterval maximumSpan = originalRegionRange.expandWithinContig(assemblyRegionArgs.maxExtensionIntoRegionPadding, sequenceDictionary);
        final SimpleInterval idealSpan = variantSpan.expandWithinContig(padding, sequenceDictionary);
        final SimpleInterval finalSpan = maximumSpan.intersect(idealSpan).mergeWithContiguous(variantSpan);
//      TODO disable this code with ERC
//        // Make double sure that, if we are emitting GVCF we won't call non-variable positions beyond the target active region span.
//        // In regular call we don't do so so we don't care and we want to maintain behavior, so the conditional.
//        final SimpleInterval callableSpan = emitReferenceConfidence ? variantSpan.intersect(originalRegionRange) : variantSpan;
        final SimpleInterval callableSpan = variantSpan;

        // TODO add equivalent debug garbage to the real assembly region trimming code
        if (HaplotypeCallerGenotypingDebugger.isEnabled()) {
            HaplotypeCallerGenotypingDebugger.println("events       : " + withinActiveRegion);
            HaplotypeCallerGenotypingDebugger.println("region       : " + originalRegion);
            HaplotypeCallerGenotypingDebugger.println("callableSpan : " + callableSpan);
            HaplotypeCallerGenotypingDebugger.println("padding      : " + padding);
            HaplotypeCallerGenotypingDebugger.println("idealSpan    : " + idealSpan);
            HaplotypeCallerGenotypingDebugger.println("maximumSpan  : " + maximumSpan);
            HaplotypeCallerGenotypingDebugger.println("finalSpan    : " + finalSpan);
        }


        return new Result(originalRegion, variantSpan, finalSpan);
    }

}
