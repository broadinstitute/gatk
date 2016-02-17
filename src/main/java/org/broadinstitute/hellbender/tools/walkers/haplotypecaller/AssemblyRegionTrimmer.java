/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Helper component to manage active region trimming
 *
 * <p/>
 * It receives the user arguments that controls trimming and also performs the trimming region calculation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ActiveRegionTrimmer {

    /**
     * Genome location parser use in order to create and manipulate genomic intervals.
     */
    private GenomeLocParser locParser;

    /**
     * Holds the debug flag. If {@code true} the trimmer will output debugging messages into the log.
     */
    private boolean debug;

    /**
     * Holds the extension to be used based on whether GGA mode is on or off.
     */
    private int usableExtension;

    /**
     * Records whether the trimming intervals are going to be used to emit reference confidence, {@code true},
     * or regular HC output {@code false}.
     */
    private boolean emitReferenceConfidence;

    @Advanced
    @Argument(fullName="dontTrimActiveRegions", shortName="dontTrimActiveRegions", doc="If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping", required = false)
    protected boolean dontTrimActiveRegions = false;

    /**
     * the maximum extent into the full active region extension that we're willing to go in genotyping our events
     */
    @Hidden
    @Argument(fullName="maxDiscARExtension", shortName="maxDiscARExtension", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping our events for discovery", required=false)
    protected int discoverExtension = 25;

    @Hidden
    @Argument(fullName="maxGGAARExtension", shortName="maxGGAARExtension", doc = "the maximum extent into the full active region extension that we're willing to go in genotyping our events for GGA mode", required=false)
    protected int ggaExtension = 300;

    /**
     * Include at least this many bases around an event for calling it
     */
    @Hidden
    @Argument(fullName="paddingAroundIndels", shortName="paddingAroundIndels", doc = "Include at least this many bases around an event for calling indels", required=false)
    public int indelPadding = 150;

    @Hidden
    @Argument(fullName="paddingAroundSNPs", shortName="paddingAroundSNPs", doc = "Include at least this many bases around an event for calling snps", required=false)
    public int snpPadding = 20;

    /**
     * Holds a reference the trimmer logger.
     */
    private final static Logger logger = Logger.getLogger(ActiveRegionTrimmer.class);

    /**
     * Initializes the trimmer.
     *
     * <p/>
     * This method should be called once and only once before any trimming is performed.
     *
     *
     * @param glp the genome-location-parser to be used when operating with genomic locations.
     * @param debug whether to show extra debug log messages.
     * @param isGGA whether the trimming region calculator should act as if we are in GGA mode or not.
     * @param emitReferenceConfidence indicates whether we plan to use this trimmer to generate trimmed regions
     *                                to be used for emitting reference confidence.
     *
     * @throws IllegalStateException if this trim calculator has already been initialized.
     * @throws IllegalArgumentException if the input location parser is {@code null}.
     * @throws UserException.BadArgumentValue if any of the user argument values is invalid.
     */
    public void initialize(final GenomeLocParser glp, final boolean debug, final boolean isGGA, final boolean emitReferenceConfidence) {
        if (locParser != null)
            throw new IllegalStateException(getClass().getSimpleName() + " instance initialized twice");
        if (glp == null)
            throw new IllegalArgumentException("input genome-loc-parser cannot be null");
        checkUserArguments();
        locParser = glp;
        this.debug = debug;
        usableExtension = isGGA ? ggaExtension : discoverExtension;
        this.emitReferenceConfidence = emitReferenceConfidence;
    }

    /**
     * Checks user trimming argument values
     *
     * @throws UserException.BadArgumentValue if there is some problem with any of the arguments values.
     */
    private void checkUserArguments() {
        if ( snpPadding < 0 ) throw new UserException.BadArgumentValue("paddingAroundSNPs","" + snpPadding + "< 0");
        if ( indelPadding < 0 ) throw new UserException.BadArgumentValue("paddingAroundIndels","" + indelPadding + "< 0");
        if ( discoverExtension < 0) throw new UserException.BadArgumentValue("maxDiscARExtension","" + discoverExtension + "< 0");
        if ( ggaExtension < 0) throw new UserException.BadArgumentValue("maxGGAAREExtension","" + ggaExtension + "< 0");
    }

    /**
     * Holds the result of trimming.
     *
     *
     *
     */
    public static class Result {

        /**
         * Indicates whether trimming is required per data and user request.
         */
        protected final boolean needsTrimming;

        /**
         * Holds the input active region.
         */
        protected final ActiveRegion originalRegion;

        /**
         * Holds the smaller range that contain all relevant callable variants in the
         * input active region (not considering the extension).
         *
         */
        protected final GenomeLoc callableSpan;

        /**
         * Maximum available range for the trimmed variant region.
         */
        protected final GenomeLoc maximumSpan;

        /**
         * The trimmed variant region span including the extension.
         */
        protected final GenomeLoc extendedSpan;


        /**
         * The ideal trimmer variant region span including the extension.
         */
        protected final GenomeLoc idealSpan;

        /**
         * Returns the ideal trimming span.
         *
         * <p/>
         * The ideal span is the one containing all callable variation overlapping the original active region span
         * (without extension) and the applicable padding {@link #getPadding()} in both sides.
         *
         *
         * @return never {@code null}.
         */
        @SuppressWarnings("unused")
        public GenomeLoc getIdealSpan() {
            return idealSpan;
        }

        /**
         * Holds the flanking spans that do not contain the callable variants.
         * <p/>
         * The first element of the pair is the left (up-stream) non-variant flank, whereas the second element is
         * the right (down-stream) non-variant flank.
         */
        protected final Pair<GenomeLoc,GenomeLoc> nonVariantFlanks;

        /**
         * Holds the collection of callable events within the variant trimming region.
         */
        protected final List<VariantContext> callableEvents;

        /**
         * Required padding around the variant trimming region.
         */
        protected final int padding;


        /**
         * Returns the required padding around callable variation.
         *
         * <p/>
         * Notice that due to the limiting span of the original active region (including its extension) it
         * is possible that the resulting final trimmed variant region span does not satisfies the padding. However
         * that should be rare.
         *
         * @return 0 or greater.
         */
        @SuppressWarnings("unused")
        public int getPadding() {
            return padding;
        }

        /**
         * Holds the maximum extension around the original active region span considered for the trimmed
         * variation region.
         */
        protected final int usableExtension;

        /**
         * Returns the maximum extension around the original active region span considered for the trimmed
         * variation region.
         *
         * <p/>
         * From time to time, the trimmed region may require a span beyond the input original active region's.
         * For example when there is a callable event close ot one of its ends and the required padding makes it
         * round beyond that limit.
         *
         * <p/>
         * Notice that due to the limiting span of the original active region (including its extended region) it
         * is possible that the resulting final trimmed variant region span goes beyond this extension including more of
         * the original active region own extension.
         *
         * @return 0 or greater.
         */
        @SuppressWarnings("unused")
        public int getUsableExtension() {
            return usableExtension;
        }

        /**
         * Holds variant-containing callable region.
         * <p/>
         * This is lazy-initialized using {@link #callableSpan}.
         */
        protected ActiveRegion callableRegion;


        /**
         * Non-variant left flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getFirst() getFirst()}.
         */
        private ActiveRegion leftFlankRegion;

        /**
         * Non-variant right flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getFirst() getSecond()}.
         */
        private ActiveRegion rightFlankRegion;

        /**
         * Whether the variant trimmed region is going to be used for emitting reference confidence records.
         */
        private final boolean emitReferenceConfidence;

        /**
         * Creates a trimming result given all its properties.
         *
         * @param emitReferenceConfidence whether reference confidence output modes are on.
         * @param needsTrimming whether there is any trimming needed at all.
         * @param originalRegion the original active region.
         * @param padding padding around contained callable variation events.
         * @param extension the extension applied to the trimmed variant span.
         * @param overlappingEvents contained callable variation events.
         * @param nonVariantFlanks pair of non-variant flank spans around the variant containing span.
         * @param extendedSpan final trimmed variant span including the extension.
         * @param idealSpan the ideal span, that contains.
         * @param maximumSpan maximum possible trimmed span based on the input original active region extended span.
         * @param callableSpan variant containing span without padding.
         */
        protected Result(final boolean emitReferenceConfidence, final boolean needsTrimming, final ActiveRegion originalRegion,
                      final int padding, final int extension,
                      final List<VariantContext> overlappingEvents, final Pair<GenomeLoc,GenomeLoc> nonVariantFlanks,
                      final GenomeLoc extendedSpan,
                      final GenomeLoc idealSpan,
                      final GenomeLoc maximumSpan,
                      final GenomeLoc callableSpan) {
            this.emitReferenceConfidence = emitReferenceConfidence;
            this.needsTrimming = needsTrimming;
            this.originalRegion = originalRegion;
            this.nonVariantFlanks = nonVariantFlanks;
            this.padding = padding;
            this.usableExtension = extension;
            this.callableEvents = overlappingEvents;
            this.callableSpan = callableSpan;
            this.idealSpan = idealSpan;
            this.maximumSpan = maximumSpan;
            this.extendedSpan = extendedSpan;

            if (!extendedSpan.isUnmapped() && !callableSpan.isUnmapped() && !extendedSpan.containsP(callableSpan))
                throw new IllegalArgumentException("the extended callable span must include the callable span");
        }


        /**
         * Checks whether there is any variation present in the target region.
         *
         * @return {@code true} if there is any variant, {@code false} otherwise.
         */
        public boolean isVariationPresent() {
            return ! callableEvents.isEmpty();
        }

        /**
         * Checks whether the active region needs trimming.
         */
        public boolean needsTrimming() {
            return needsTrimming;
        }

        /**
         * Returns the trimmed variant containing region
         *
         * @throws IllegalStateException if there is no variation detected.
         *
         * @return never {@code null}.
         */
        public ActiveRegion getCallableRegion() {
            if (callableRegion == null && !extendedSpan.isUnmapped())
                //TODO this conditional is a patch to retain the current standard HC run behaviour
                //TODO we should simply remove this difference between trimming with or without GVCF
                //TODO embracing slight changes in the standard HC output
                callableRegion = emitReferenceConfidence ? originalRegion.trim(callableSpan, extendedSpan) : originalRegion.trim(extendedSpan);
            else if (extendedSpan.isUnmapped())
                throw new IllegalStateException("there is no variation thus no variant region");
            return callableRegion;
        }

        /**
         * Checks whether there is a non-empty left flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial left flank region, {@code false} otherwise.
         */
        public boolean hasLeftFlankingRegion() {
            return ! nonVariantFlanks.getFirst().isUnmapped();
        }

        /**
         * Checks whether there is a non-empty right flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial right flank region, {@code false} otherwise.
         */
        public boolean hasRightFlankingRegion() {
            return ! nonVariantFlanks.getSecond().isUnmapped();
        }

        /**
         *  Returns the trimmed out left non-variant region.
         *  <p/>
         *  Notice that in case of no variation, the whole original region is considered the left flanking region.
         *
         *  @throws IllegalStateException if there is not such as left flanking region.
         */
        public ActiveRegion nonVariantLeftFlankRegion() {
            if (leftFlankRegion == null && ! nonVariantFlanks.getFirst().isUnmapped())
                leftFlankRegion = originalRegion.trim(nonVariantFlanks.getFirst(),originalRegion.getExtension());
            else if (nonVariantFlanks.getFirst().isUnmapped())
                throw new IllegalStateException("there is no left flank non-variant trimmed out region");
            return leftFlankRegion;
        }

        /**
         *  Returns the trimmed out right non-variant region.
         */
        public ActiveRegion nonVariantRightFlankRegion() {
            if (rightFlankRegion == null && ! nonVariantFlanks.getSecond().isUnmapped())
                rightFlankRegion = originalRegion.trim(nonVariantFlanks.getSecond(),originalRegion.getExtension());
            else if (nonVariantFlanks.getSecond().isUnmapped())
                throw new IllegalStateException("there is no right flank non-variant trimmed out region");
            return rightFlankRegion;
        }

        /**
         * Creates a result indicating that there was no trimming to be done.
         */
        protected static Result noTrimming(final boolean emitReferenceConfidence,
                                           final ActiveRegion targetRegion, final int padding,
                                           final int usableExtension,final List<VariantContext> events) {
            final GenomeLoc targetRegionLoc = targetRegion.getLocation();
            final Result result = new Result(emitReferenceConfidence,false,targetRegion,padding,usableExtension,events,new Pair<>(GenomeLoc.UNMAPPED,GenomeLoc.UNMAPPED),
                    targetRegionLoc,targetRegionLoc,targetRegionLoc,targetRegionLoc);
            result.callableRegion = targetRegion;
            return result;
        }

        /**
         * Creates a result indicating that no variation was found.
         */
        protected static Result noVariation(final boolean emitReferenceConfidence, final ActiveRegion targetRegion,
                                            final int padding, final int usableExtension) {
            final Result result = new Result(emitReferenceConfidence,false,targetRegion,padding,usableExtension,
                    Collections.<VariantContext>emptyList(),new Pair<>(targetRegion.getLocation(),GenomeLoc.UNMAPPED),
                    GenomeLoc.UNMAPPED,GenomeLoc.UNMAPPED,GenomeLoc.UNMAPPED,GenomeLoc.UNMAPPED);
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
    public Result trim(final ActiveRegion originalRegion,
                                            final TreeSet<VariantContext> allVariantsWithinExtendedRegion) {


        if ( allVariantsWithinExtendedRegion.isEmpty() ) // no variants,
            return Result.noVariation(emitReferenceConfidence,originalRegion,snpPadding, usableExtension);

        final List<VariantContext> withinActiveRegion = new LinkedList<>();
        final GenomeLoc originalRegionRange = originalRegion.getLocation();
        boolean foundNonSnp = false;
        GenomeLoc variantSpan = null;
        for ( final VariantContext vc : allVariantsWithinExtendedRegion ) {
            final GenomeLoc vcLoc = locParser.createGenomeLoc(vc);
            if ( originalRegionRange.overlapsP(vcLoc) ) {
                foundNonSnp = foundNonSnp || ! vc.isSNP();
                variantSpan = variantSpan == null ? vcLoc : variantSpan.endpointSpan(vcLoc);
                withinActiveRegion.add(vc);
            }
        }
        final int padding = foundNonSnp ? indelPadding : snpPadding;

        // we don't actually have anything in the region after skipping out variants that don't overlap
        // the region's full location
        if ( variantSpan == null )
            return Result.noVariation(emitReferenceConfidence,originalRegion,padding, usableExtension);

        if ( dontTrimActiveRegions)
            return Result.noTrimming(emitReferenceConfidence,originalRegion, padding, usableExtension, withinActiveRegion);

        final GenomeLoc maximumSpan = locParser.createPaddedGenomeLoc(originalRegionRange, usableExtension);
        final GenomeLoc idealSpan = locParser.createPaddedGenomeLoc(variantSpan, padding);
        final GenomeLoc finalSpan = maximumSpan.intersect(idealSpan).union(variantSpan);

        // Make double sure that, if we are emitting GVCF we won't call non-variable positions beyond the target active region span.
        // In regular call we don't do so so we don't care and we want to maintain behavior, so the conditional.
        final GenomeLoc callableSpan = emitReferenceConfidence ? variantSpan.intersect(originalRegionRange) : variantSpan;

        final Pair<GenomeLoc,GenomeLoc> nonVariantRegions = nonVariantTargetRegions(originalRegion, callableSpan);

        if ( debug ) {
            logger.info("events       : " + withinActiveRegion);
            logger.info("region       : " + originalRegion);
            logger.info("callableSpan : " + callableSpan);
            logger.info("padding      : " + padding);
            logger.info("idealSpan    : " + idealSpan);
            logger.info("maximumSpan  : " + maximumSpan);
            logger.info("finalSpan    : " + finalSpan);
        }

        return new Result(emitReferenceConfidence,true,originalRegion,padding, usableExtension,withinActiveRegion,nonVariantRegions,finalSpan,idealSpan,maximumSpan,variantSpan);
    }

    /**
     * Calculates the list of region to trim away.
     * @param targetRegion region for which to generate the flanking regions.
     * @param variantSpan the span of the core region containing relevant variation and required padding.
     * @return never {@code null}; 0, 1 or 2 element list.
     */
    private Pair<GenomeLoc,GenomeLoc> nonVariantTargetRegions(final ActiveRegion targetRegion, final GenomeLoc variantSpan) {
        final GenomeLoc targetRegionRange = targetRegion.getLocation();
        final int finalStart = variantSpan.getStart();
        final int finalStop = variantSpan.getStop();

        final int targetStart = targetRegionRange.getStart();
        final int targetStop = targetRegionRange.getStop();

        final boolean preTrimmingRequired = targetStart < finalStart;
        final boolean postTrimmingRequired = targetStop > finalStop;
        if (preTrimmingRequired) {
            final String contig = targetRegionRange.getContig();
            return postTrimmingRequired ? new Pair<>(
                    locParser.createGenomeLoc(contig, targetStart, finalStart - 1),
                    locParser.createGenomeLoc(contig, finalStop + 1, targetStop)) :
                    new Pair<>(locParser.createGenomeLoc(contig, targetStart, finalStart - 1),GenomeLoc.UNMAPPED);
        } else if (postTrimmingRequired)
            return new Pair<>(GenomeLoc.UNMAPPED,locParser.createGenomeLoc(targetRegionRange.getContig(), finalStop + 1, targetStop));
        else
            return new Pair<>(GenomeLoc.UNMAPPED,GenomeLoc.UNMAPPED);
    }
}