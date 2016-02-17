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

package org.broadinstitute.gatk.tools.walkers.genotyper;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContextUtils;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.gatk.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.gga.GenotypingGivenAllelesUtils;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Base class for genotyper engines.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class GenotypingEngine<Config extends StandardCallerArgumentCollection> {

    protected final AFCalculatorProvider afCalculatorProvider   ;

    protected Logger logger;

    protected final Config configuration;

    protected VariantAnnotatorEngine annotationEngine;

    protected final int numberOfGenomes;

    protected final SampleList samples;


    private final AFPriorProvider log10AlleleFrequencyPriorsSNPs;

    private final AFPriorProvider log10AlleleFrequencyPriorsIndels;

    protected final GenomeLocParser genomeLocParser;

    /**
     * Construct a new genotyper engine, on a specific subset of samples.
     *
     * @param configuration engine configuration object.
     * @param samples subset of sample to work on identified by their names. If {@code null}, the full toolkit
     *                    sample set will be used instead.
     * @param genomeLocParser the genome-loc-parser
     *
     * @throws IllegalArgumentException if any of {@code samples}, {@code configuration} or {@code genomeLocParser} is {@code null}.
     */
    protected GenotypingEngine(final Config configuration, final SampleList samples,
                               final GenomeLocParser genomeLocParser, final AFCalculatorProvider afCalculatorProvider) {

        if (configuration == null)
            throw new IllegalArgumentException("the configuration cannot be null");
        if (samples == null)
            throw new IllegalArgumentException("the sample list provided cannot be null");
        if (afCalculatorProvider == null)
            throw new IllegalArgumentException("the AF calculator provider cannot be null");

        this.afCalculatorProvider = afCalculatorProvider;
        this.configuration = configuration;
        logger = Logger.getLogger(getClass());
        this.samples = samples;
        numberOfGenomes = this.samples.sampleCount() * configuration.genotypeArgs.samplePloidy;
        MathUtils.Log10Cache.ensureCacheContains(numberOfGenomes * 2);
        log10AlleleFrequencyPriorsSNPs = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.snpHeterozygosity, configuration.genotypeArgs.inputPrior);
        log10AlleleFrequencyPriorsIndels = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.indelHeterozygosity, configuration.genotypeArgs.inputPrior);
        this.genomeLocParser = genomeLocParser;
    }

    /**
     * Changes the logger for this genotyper engine.
     *
     * @param logger new logger.
     *
     * @throws IllegalArgumentException if {@code logger} is {@code null}.
     */
    public void setLogger(final Logger logger) {
        if (logger == null)
            throw new IllegalArgumentException("the logger cannot be null");
        this.logger = logger;
    }

    public Set<VCFInfoHeaderLine> getAppropriateVCFInfoHeaders() {
        Set<VCFInfoHeaderLine> headerInfo = new HashSet<>();
        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED )
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
        return headerInfo;
    }

    /**
     * Changes the annotation engine for this genotyping-engine.
     *
     * @param annotationEngine the new annotation engine (can be {@code null}).
     */
    public void setAnnotationEngine(final VariantAnnotatorEngine annotationEngine) {
        this.annotationEngine = annotationEngine;
    }

    /**
     * Returns a reference to the engine configuration
     *
     * @return never {@code null}.
     */
    public Config getConfiguration() {
        return configuration;
    }

    /**
     * Completes a variant context with genotype calls and associated annotations given the genotype likelihoods and
     *  the model that need to be applied.
     *
     * @param vc variant-context to complete.
     * @param model model name.
     *
     * @throws IllegalArgumentException if {@code model} or {@code vc} is {@code null}.
     *
     * @return can be {@code null} indicating that genotyping it not possible with the information provided.
     */
    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model) {
        if (vc == null)
            throw new IllegalArgumentException("vc cannot be null");
        if (model == null)
            throw new IllegalArgumentException("the model cannot be null");
        return calculateGenotypes(null,null,null,null,vc,model,false,null);
    }

    /**
     * Main entry function to calculate genotypes of a given VC with corresponding GL's that is shared across genotypers (namely UG and HC).
     *
     * @param tracker                            Tracker
     * @param refContext                         Reference context
     * @param rawContext                         Raw context
     * @param stratifiedContexts                 Stratified alignment contexts
     * @param vc                                 Input VC
     * @param model                              GL calculation model
     * @param inheritAttributesFromInputVC       Output VC will contain attributes inherited from input vc
     * @return                                   VC with assigned genotypes
     */
    protected VariantCallContext calculateGenotypes(final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                 final AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc, final GenotypeLikelihoodsCalculationModel.Model model,
                                                 final boolean inheritAttributesFromInputVC,
                                                 final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;
        // if input VC can't be genotyped, exit with either null VCC or, in case where we need to emit all sites, an empty call
        if (hasTooManyAlternativeAlleles(vc) || vc.getNSamples() == 0)
            return emptyCallContext(tracker,refContext,rawContext);

        final int defaultPloidy = configuration.genotypeArgs.samplePloidy;
        final int maxAltAlleles = configuration.genotypeArgs.MAX_ALTERNATE_ALLELES;
        final AFCalculator afCalculator = afCalculatorProvider.getInstance(vc,defaultPloidy,maxAltAlleles);
        final AFCalculationResult AFresult = afCalculator.getLog10PNonRef(vc, defaultPloidy,maxAltAlleles, getAlleleFrequencyPriors(vc,defaultPloidy,model));

        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult);

        final double PoFGT0 = Math.pow(10, AFresult.getLog10PosteriorOfAFGT0());

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        double log10Confidence =
                ! outputAlternativeAlleles.siteIsMonomorphic ||
                    configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES || configuration.annotateAllSitesWithPLs
                        ? AFresult.getLog10PosteriorOfAFEq0() + 0.0
                        : AFresult.getLog10PosteriorOfAFGT0() + 0.0 ;


        // Add 0.0 removes -0.0 occurrences.
        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic) && !forceSiteEmission()) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            double[] AFpriors = getAlleleFrequencyPriors(vc, defaultPloidy, model);
            final int INDEX_FOR_AC_EQUALS_1 = 1;
            return limitedContext ? null : estimateReferenceConfidence(vc, stratifiedContexts, AFpriors[INDEX_FOR_AC_EQUALS_1], true, PoFGT0);
        }

        // start constructing the resulting VC
        final GenomeLocParser genomeLocParser = this.genomeLocParser != null || refContext == null ? this.genomeLocParser : refContext.getGenomeLocParser();
        if (genomeLocParser == null)
            throw new IllegalStateException("this UG engine was created without a valid genomeLocParser and no refContext was provided");
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(vc);
        final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
        final VariantContextBuilder builder = new VariantContextBuilder(callSourceString(), loc.getContig(), loc.getStart(), loc.getStop(), outputAlleles);
        
        builder.log10PError(log10Confidence);
        if ( ! passesCallThreshold(phredScaledConfidence) )
            builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);

        // create the genotypes

        final GenotypesContext genotypes = afCalculator.subsetAlleles(vc, defaultPloidy, outputAlleles, true);
        builder.genotypes(genotypes);

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        final Map<String, Object> attributes = composeCallAttributes(inheritAttributesFromInputVC, vc, rawContext, stratifiedContexts, tracker, refContext,
                outputAlternativeAlleles.alternativeAlleleMLECounts(), outputAlternativeAlleles.siteIsMonomorphic, AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()),genotypes,model,perReadAlleleLikelihoodMap);

        builder.attributes(attributes);

        VariantContext vcCall = builder.make();

        if ( annotationEngine != null && !limitedContext ) { // limitedContext callers need to handle annotations on their own by calling their own annotationEngine
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadBackedPileup pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(pileup);

            vcCall = annotationEngine.annotateContext(tracker, refContext, stratifiedContexts, vcCall, perReadAlleleLikelihoodMap);
        }

        // if we are subsetting alleles (either because there were too many or because some were not polymorphic)
        // then we may need to trim the alleles (because the original VariantContext may have had to pad at the end).
        if ( outputAlleles.size() != vc.getAlleles().size() && !limitedContext ) // limitedContext callers need to handle allele trimming on their own to keep their perReadAlleleLikelihoodMap alleles in sync
            vcCall = GATKVariantContextUtils.reverseTrimAlleles(vcCall);

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, PoFGT0));
    }

    /**
     * What string to use as source of variant-context generated by this genotyper-engine.
     * @return never {@code null} nor empty.
     */
    protected abstract String callSourceString();

    /**
     * Holds information about the alternative allele subsetting based on supporting evidence, genotyping and
     * output modes.
     */
    private static class OutputAlleleSubset {
        private  final Allele[] alleles;
        private  final boolean siteIsMonomorphic;
        private  final int[] mleCounts;
        private  final int count;

        @Requires({"count >= 0", "alleles != null", "mleCounts != null", "count <= alleles.length", " count <= mleCounts"})
        private OutputAlleleSubset(final int count, final Allele[] alleles, final int[] mleCounts, final boolean siteIsMonomorphic) {
            this.siteIsMonomorphic = siteIsMonomorphic;
            this.count = count;
            this.alleles = alleles;
            this.mleCounts = mleCounts;
        }

        private List<Allele> outputAlleles(final Allele referenceAllele) {
            final ArrayList<Allele> result = new ArrayList<>(count + 1);
            result.add(referenceAllele);
            for (int i = 0; i < count; i++)
                result.add(alleles[i]);
            return result;
        }

        public List<Integer> alternativeAlleleMLECounts() {
            final List<Integer> result = new ArrayList<>(count);
            for (int i = 0; i < count; i++)
                result.add(mleCounts[i]);
            return result;
        }
    }


    /**
     * Provided the exact mode computations it returns the appropiate subset of alleles that progress to genotyping.
     * @param afcr the exact model calcualtion result.
     * @return never {@code null}.
     */
    private OutputAlleleSubset calculateOutputAlleleSubset(final AFCalculationResult afcr) {
        final List<Allele> alleles = afcr.getAllelesUsedInGenotyping();

        final int alternativeAlleleCount = alleles.size() - 1;
        Allele[] outputAlleles = new Allele[alternativeAlleleCount];
        int[] mleCounts = new int[alternativeAlleleCount];
        int outputAlleleCount = 0;
        boolean siteIsMonomorphic = true;
        for (final Allele alternativeAllele : alleles) {
            if (alternativeAllele.isReference()) continue;
            final boolean isPlausible = afcr.isPolymorphicPhredScaledQual(alternativeAllele, configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING);
            final boolean toOutput = isPlausible || forceKeepAllele(alternativeAllele);

            siteIsMonomorphic &= ! isPlausible;
            if (!toOutput) continue;
            outputAlleles[outputAlleleCount] = alternativeAllele;
            mleCounts[outputAlleleCount++] = afcr.getAlleleCountAtMLE(alternativeAllele);
        }

        return new OutputAlleleSubset(outputAlleleCount,outputAlleles,mleCounts,siteIsMonomorphic);
    }

    /**
     * Checks whether even if the allele is not well supported by the data, we should keep it for genotyping.
     *
     * @param allele target allele.
     *
     * @return {@code true} iff we need to keep this alleles even if does not seem plausible.
     */
    protected abstract boolean forceKeepAllele(final Allele allele);

    /**
     * Checks whether a variant site seems confidently called base on user threshold that the score provided
     * by the exact model.
     *
     * @param conf
     * @param PofF
     * @return {@code true} iff the variant is confidently called.
     */
    protected final boolean confidentlyCalled(final double conf, final double PofF) {
        return conf >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING ||
                (configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
                        && QualityUtils.phredScaleErrorRate(PofF) >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING);
    }


    /**
     * Checks whether the variant context has too many alternative alleles for progress to genotyping the site.
     * <p>
     *     AF calculation may get intro trouble with too many alternative alleles.
     * </p>
     *
     * @param vc the variant context to evaluate.
     *
     * @throws NullPointerException if {@code vc} is {@code null}.
     *
     * @return {@code true} iff there is too many alternative alleles based on
     * {@link GenotypeLikelihoods#MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED}.
     */
    @Requires("vc != null")
    protected final boolean hasTooManyAlternativeAlleles(final VariantContext vc) {
        // protect against too many alternate alleles that we can't even run AF on:
        if (vc.getNAlleles() <= GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
            return false;
        logger.warn("Attempting to genotype more than "+GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED +
             " alleles. Site will be skipped at location "+vc.getChr()+":"+vc.getStart());
        return true;
    }

    /**
     * Produces an empty variant-call context to output when there is no enough data provided to call anything.
     *
     * @param tracker meta-data tracker.
     * @param ref the reference context.
     * @param rawContext the read alignment at that location.
     * @return it might be null if no enough information is provided to do even an empty call. For example when
     * we have limited-context (i.e. any of the tracker, reference or alignment is {@code null}.
     */
    protected final VariantCallContext emptyCallContext(final RefMetaDataTracker tracker, final ReferenceContext ref,
                                                      final AlignmentContext rawContext) {
        if (tracker == null || ref == null || rawContext == null)
            return null;

        if (!forceSiteEmission())
            return null;

        VariantContext vc;

        if ( configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext ggaVc = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(tracker,
                    rawContext.getLocation(), false, logger, configuration.alleles);
            if (ggaVc == null)
                return null;
            vc = new VariantContextBuilder(callSourceString(), ref.getLocus().getContig(), ggaVc.getStart(),
                        ggaVc.getEnd(), ggaVc.getAlleles()).make();
        } else {
            // deal with bad/non-standard reference bases
            if ( !Allele.acceptableAlleleBases(new byte[]{ref.getBase()}) )
               return null;
            final Set<Allele> alleles = new HashSet<>(Collections.singleton(Allele.create(ref.getBase(),true)));
            vc = new VariantContextBuilder(callSourceString(), ref.getLocus().getContig(),
                   ref.getLocus().getStart(), ref.getLocus().getStart(), alleles).make();
        }

        if ( vc != null && annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadBackedPileup pileup = rawContext.getBasePileup();
            vc = annotationEngine.annotateContext(tracker, ref,  AlignmentContextUtils.splitContextBySampleName(pileup), vc);
        }

        return new VariantCallContext(vc, false);
    }

    /**
     * Indicates whether we have to emit any site no matter what.
     * <p>
     *     Note: this has been added to allow differences between UG and HC GGA modes where the latter force emmitions of all given alleles
     *     sites even if there is no enough confidence.
     * </p>
     *
     * @return {@code true} iff we force emissions.
     */
    protected abstract boolean forceSiteEmission();

    protected final VariantCallContext estimateReferenceConfidence(VariantContext vc, Map<String, AlignmentContext> contexts, double log10OfTheta, boolean ignoreCoveredSamples, double initialPofRef) {
        if ( contexts == null )
            return null;

        double log10POfRef = Math.log10(initialPofRef);

        // for each sample that we haven't examined yet
        final int sampleCount = samples.sampleCount();
        for (int i = 0; i < sampleCount; i++) {
            final String sample = samples.sampleAt(i);
            final AlignmentContext context = contexts.get(sample);
            if ( ignoreCoveredSamples && context != null )
                continue;
            final int depth = context == null ? 0 : context.getBasePileup().depthOfCoverage();
            log10POfRef += estimateLog10ReferenceConfidenceForOneSample(depth, log10OfTheta);
        }

        return new VariantCallContext(vc, QualityUtils.phredScaleLog10CorrectRate(log10POfRef) >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING, false);
    }

    /**
     * Returns the log10 prior probability for all possible allele counts from 0 to N where N is the total number of
     * genomes (total-ploidy).
     *
     * @param vc the target variant-context, use to determine the total ploidy thus the possible ACs.
     * @param defaultPloidy default ploidy to be assume if we do not have the ploidy for some sample in {@code vc}.
     * @param model the calculation model (SNP,INDEL or MIXED) whose priors are to be retrieved.
     * @throws java.lang.NullPointerException if either {@code vc} or {@code model} is {@code null}
     * @return never {@code null}, an array with exactly <code>total-ploidy(vc) + 1</code> positions.
     */
    protected final double[] getAlleleFrequencyPriors( final VariantContext vc, final int defaultPloidy, final GenotypeLikelihoodsCalculationModel.Model model ) {
        final int totalPloidy = GATKVariantContextUtils.totalPloidy(vc, defaultPloidy);
        switch (model) {
            case SNP:
            case GENERALPLOIDYSNP:
                return log10AlleleFrequencyPriorsSNPs.forTotalPloidy(totalPloidy);
            case INDEL:
            case GENERALPLOIDYINDEL:
                return log10AlleleFrequencyPriorsIndels.forTotalPloidy(totalPloidy);
            default:
                throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + model);
        }
    }

    /**
     * Compute the log10 probability of a sample with sequencing depth and no alt allele is actually truly homozygous reference
     *
     * Assumes the sample is diploid
     *
     * @param depth the depth of the sample
     * @param log10OfTheta the heterozygosity of this species (in log10-space)
     *
     * @return a valid log10 probability of the sample being hom-ref
     */
    @Ensures("MathUtils.goodLog10Probability(result)")
    protected final double estimateLog10ReferenceConfidenceForOneSample(final int depth, final double log10OfTheta) {
        final double log10PofNonRef = log10OfTheta + getRefBinomialProbLog10(depth);
        return MathUtils.log10OneMinusX(Math.pow(10.0, log10PofNonRef));
    }

    /**
     * Calculates the hom-reference binomial log10 probability given the depth.
     *
     * @param depth the query depth.
     *
     * @throws IllegalArgumentException if {@code depth} is less than 0.
     *
     * @return a valid log10 probability between 0 and {@link Double#NEGATIVE_INFINITY}.
     */
    protected final double getRefBinomialProbLog10(final int depth) {
        if (depth < 0)
            throw new IllegalArgumentException("depth cannot be less than 0");
        return MathUtils.log10BinomialProbability(depth, 0);
    }

    /**
     * Function that fills vector with allele frequency priors. By default, infinite-sites, neutral variation prior is used,
     * where Pr(AC=i) = theta/i where theta is heterozygosity
     * @param N                                Number of chromosomes
     * @param heterozygosity                   default heterozygosity to use, if inputPriors is empty
     * @param inputPriors                      Input priors to use (in which case heterozygosity is ignored)
     *
     * @throws IllegalArgumentException if {@code inputPriors} has size != {@code N} or any entry in {@code inputPriors} is not in the (0,1) range.
     *
     * @return never {@code null}.
     */
    public static AFPriorProvider composeAlleleFrequencyPriorProvider(final int N, final double heterozygosity, final List<Double> inputPriors) {

        if (!inputPriors.isEmpty()) {
            // user-specified priors
            if (inputPriors.size() != N)
                throw new UserException.BadArgumentValue("inputPrior","Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");
            for (final Double prior : inputPriors) {
                if (prior <= 0 || prior >= 1) throw new UserException.BadArgumentValue("inputPrior","inputPrior vector values must be greater than 0 and less than 1");
            }
            return new CustomAFPriorProvider(inputPriors);
        }
        else
            return new HeterozygosityAFPriorProvider(heterozygosity);
    }

    /**
     * Function that fills vector with allele frequency priors. By default, infinite-sites, neutral variation prior is used,
     * where Pr(AC=i) = theta/i where theta is heterozygosity
     * @param N                                Number of chromosomes
     * @param priors                           (output) array to be filled with priors
     * @param heterozygosity                   default heterozygosity to use, if inputPriors is empty
     * @param inputPriors                      Input priors to use (in which case heterozygosity is ignored)
     */
    public static void computeAlleleFrequencyPriors(final int N, final double[] priors, final double heterozygosity, final List<Double> inputPriors) {


        double sum = 0.0;

        if (!inputPriors.isEmpty()) {
            // user-specified priors
            if (inputPriors.size() != N)
                throw new UserException.BadArgumentValue("inputPrior","Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");

            int idx = 1;
            for (final double prior: inputPriors) {
                if (prior < 0.0)
                    throw new UserException.BadArgumentValue("Bad argument: negative values not allowed","inputPrior");
                priors[idx++] = Math.log10(prior);
                sum += prior;
            }
        }
        else {
            // for each i
            for (int i = 1; i <= N; i++) {
                final double value = heterozygosity / (double)i;
                priors[i] = Math.log10(value);
                sum += value;
            }
        }

        // protection against the case of heterozygosity too high or an excessive number of samples (which break population genetics assumptions)
        if (sum > 1.0) {
            throw new UserException.BadArgumentValue("heterozygosity","The heterozygosity value is set too high relative to the number of samples to be processed, or invalid values specified if input priors were provided - try reducing heterozygosity value or correct input priors.");
        }
        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
    }


    protected final boolean passesEmitThreshold(double conf, boolean bestGuessIsRef) {
        return (configuration.outputMode == OutputMode.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) &&
                conf >= Math.min(configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING,
                        configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING);
    }

    protected final boolean passesCallThreshold(double conf) {
        return conf >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    protected Map<String,Object> composeCallAttributes(final boolean inheritAttributesFromInputVC, final VariantContext vc,
                                                       final AlignmentContext rawContext, final Map<String, AlignmentContext> stratifiedContexts, final RefMetaDataTracker tracker, final ReferenceContext refContext, final List<Integer> alleleCountsofMLE, final boolean bestGuessIsRef,
                                                       final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes,
                                                       final GenotypeLikelihoodsCalculationModel.Model model, final Map<String, org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final HashMap<String, Object> attributes = new HashMap<>();

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        // inherit attributes from input vc if requested
        if (inheritAttributesFromInputVC)
            attributes.putAll(vc.getAttributes());
        // if the site was down-sampled, record that fact
        if ( !limitedContext && rawContext.hasPileupBeenDownsampled() )
            attributes.put(GATKVCFConstants.DOWNSAMPLED_KEY, true);

        // add the MLE AC and AF annotations
        if ( alleleCountsofMLE.size() > 0 ) {
            attributes.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
            final ArrayList<Double> MLEfrequencies = calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes);
            attributes.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
        }

        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED )
            attributes.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, vc.getAlternateAlleles().size());


        return attributes;
    }

    private ArrayList<Double> calculateMLEAlleleFrequencies(List<Integer> alleleCountsofMLE, GenotypesContext genotypes) {
        int AN = 0;
        for (final Genotype g : genotypes)
           for (final Allele a : g.getAlleles())
              if (!a.isNoCall()) AN++;

        final ArrayList<Double> MLEfrequencies = new ArrayList<>(alleleCountsofMLE.size());
        // the MLEAC is allowed to be larger than the AN (e.g. in the case of all PLs being 0, the GT is ./. but the exact model may arbitrarily choose an AC>1)
        for (final int AC : alleleCountsofMLE )
            MLEfrequencies.add(Math.min(1.0, (double)AC / (double)AN));
        return MLEfrequencies;
    }

    /**
     * Calculates the active state profile value for a single sample.
     *
     * @param log10GenotypeLikelihoods the single sample genotype likelihoods.
     * @return log10 probability from 0 to -Infinity.
     */
    public double calculateSingleSampleRefVsAnyActiveStateProfileValue(final double[] log10GenotypeLikelihoods) {
        if (log10GenotypeLikelihoods == null) {
            throw new IllegalArgumentException("the input likelihoods cannot be null");
        } else if (log10GenotypeLikelihoods.length != this.configuration.genotypeArgs.samplePloidy + 1) {
            throw new IllegalArgumentException("wrong likelihoods dimensions");
        } else {
            final double[] log10Priors = log10AlleleFrequencyPriorsSNPs.forTotalPloidy(this.configuration.genotypeArgs.samplePloidy);
            final double log10ACeq0Likelihood = log10GenotypeLikelihoods[0];
            final double log10ACeq0Prior = log10Priors[0];
            final double log10ACeq0Posterior = log10ACeq0Likelihood + log10ACeq0Prior;

            // If the Maximum a-posteriori AC is 0 then the profile value must be 0.0 as per existing code; it does
            // not matter whether a AC > 0 is at all plausible.
            boolean mapACeq0 = true;
            for (int AC = 1; AC < log10Priors.length; AC++)
                if (log10Priors[AC] + log10GenotypeLikelihoods[AC] > log10ACeq0Posterior) {
                    mapACeq0 = false;
                    break;
                }
            if (mapACeq0)
                return 0.0;

            //TODO bad way to calculate AC > 0 posterior that follows the current behaviour of ExactAFCalculator (StateTracker)
            //TODO this is the lousy part... this code just adds up lks and priors of AC != 0 before as if
            //TODO Sum(a_i * b_i) is equivalent to Sum(a_i) * Sum(b_i)
            //TODO This has to be changed not just here but also in the AFCalculators (StateTracker).
            final double log10ACgt0Likelihood = MathUtils.approximateLog10SumLog10(log10GenotypeLikelihoods, 1, log10GenotypeLikelihoods.length);
            final double log10ACgt0Prior = MathUtils.approximateLog10SumLog10(log10Priors, 1, log10Priors.length);
            final double log10ACgt0Posterior = log10ACgt0Likelihood + log10ACgt0Prior;
            final double log10PosteriorNormalizationConstant = MathUtils.approximateLog10SumLog10(log10ACeq0Posterior, log10ACgt0Posterior);
            //TODO End of lousy part.

            final double normalizedLog10ACeq0Posterior = log10ACeq0Posterior - log10PosteriorNormalizationConstant;
            // This is another condition to return a 0.0 also present in AFCalculator code as well.
            if (normalizedLog10ACeq0Posterior >= QualityUtils.qualToErrorProbLog10(configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_EMITTING))
                return 0.0;

            return 1.0 - Math.pow(10.0, normalizedLog10ACeq0Posterior);
        }
    }
}
