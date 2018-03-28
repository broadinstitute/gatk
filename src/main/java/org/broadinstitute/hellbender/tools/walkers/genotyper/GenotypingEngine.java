package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Base class for genotyper engines.
 */
public abstract class GenotypingEngine<Config extends StandardCallerArgumentCollection> {

    protected final AFCalculator newAFCalculator;

    protected final AFCalculatorProvider afCalculatorProvider;

    protected final Config configuration;

    protected VariantAnnotatorEngine annotationEngine;

    protected Logger logger;

    protected final int numberOfGenomes;

    protected final SampleList samples;

    private final AFPriorProvider log10AlleleFrequencyPriorsSNPs;

    private final AFPriorProvider log10AlleleFrequencyPriorsIndels;

    private final List<SimpleInterval> upstreamDeletionsLoc = new LinkedList<>();

    private final boolean doAlleleSpecificCalcs;

    /**
     * Construct a new genotyper engine, on a specific subset of samples.
     *
     * @param configuration engine configuration object.
     * @param samples subset of sample to work on identified by their names. If {@code null}, the full toolkit
     *                    sample set will be used instead.
     * @param doAlleleSpecificCalcs Whether the AS_QUAL key should be calculated and added to newly genotyped variants.
     *
     * @throws IllegalArgumentException if any of {@code samples}, {@code configuration} is {@code null}.
     */
    protected GenotypingEngine(final Config configuration,
                               final SampleList samples,
                               final AFCalculatorProvider afCalculatorProvider,
                               final boolean doAlleleSpecificCalcs) {
        this.configuration = Utils.nonNull(configuration, "the configuration cannot be null");
        this.samples = Utils.nonNull(samples, "the sample list cannot be null");
        this.afCalculatorProvider = Utils.nonNull(afCalculatorProvider, "the AF calculator provider cannot be null");
        this.doAlleleSpecificCalcs = doAlleleSpecificCalcs;
        logger = LogManager.getLogger(getClass());
        numberOfGenomes = this.samples.numberOfSamples() * configuration.genotypeArgs.samplePloidy;
        log10AlleleFrequencyPriorsSNPs = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.snpHeterozygosity, configuration.genotypeArgs.inputPrior);
        log10AlleleFrequencyPriorsIndels = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.indelHeterozygosity, configuration.genotypeArgs.inputPrior);

        final double refPseudocount = configuration.genotypeArgs.snpHeterozygosity / Math.pow(configuration.genotypeArgs.heterozygosityStandardDeviation,2);
        final double snpPseudocount = configuration.genotypeArgs.snpHeterozygosity * refPseudocount;
        final double indelPseudocount = configuration.genotypeArgs.indelHeterozygosity * refPseudocount;
        newAFCalculator = new AlleleFrequencyCalculator(refPseudocount, snpPseudocount, indelPseudocount, configuration.genotypeArgs.samplePloidy);
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
            if (inputPriors.size() != N) {
                throw new CommandLineException.BadArgumentValue("inputPrior", "Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");
            }

            int idx = 1;
            for (final double prior: inputPriors) {
                if (prior < 0.0) {
                    throw new CommandLineException.BadArgumentValue("Bad argument: negative values not allowed", "inputPrior");
                }
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
            throw new CommandLineException.BadArgumentValue("heterozygosity","The heterozygosity value is set too high relative to the number of samples to be processed, or invalid values specified if input priors were provided - try reducing heterozygosity value or correct input priors.");
        }
        // null frequency for AF=0 is (1 - sum(all other frequencies))
        priors[0] = Math.log10(1.0 - sum);
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
            if (inputPriors.size() != N) {
                throw new CommandLineException.BadArgumentValue("inputPrior", "Invalid length of inputPrior vector: vector length must be equal to # samples +1 ");
            }
            for (final Double prior : inputPriors) {
                if (prior <= 0 || prior >= 1) {
                    throw new CommandLineException.BadArgumentValue("inputPrior", "inputPrior vector values must be greater than 0 and less than 1");
                }
            }
            return new CustomAFPriorProvider(inputPriors);
        }
        else {
            return new HeterozygosityAFPriorProvider(heterozygosity);
        }
    }

    /**
     * Changes the logger for this genotyper engine.
     *
     * @param logger new logger.
     *
     * @throws IllegalArgumentException if {@code logger} is {@code null}.
     */
    public void setLogger(final Logger logger) {
        this.logger = Utils.nonNull(logger, "the logger cannot be null");
    }

    public Set<VCFInfoHeaderLine> getAppropriateVCFInfoHeaders() {
        final Set<VCFInfoHeaderLine> headerInfo = new LinkedHashSet<>();
        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
            headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));
        }
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
    public VariantCallContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel model, final SAMFileHeader header) {
        Utils.nonNull(vc, "vc cannot be null");
        Utils.nonNull(model, "the model cannot be null");
        return calculateGenotypes(null,null,null,null,vc,model,false,null,header);
    }

    /**
     * Main entry function to calculate genotypes of a given VC with corresponding GL's that is shared across genotypers (namely UG and HC).
     *
     * @param features                           Features
     * @param refContext                         Reference context
     * @param rawContext                         Raw context
     * @param stratifiedContexts                 Stratified alignment contexts
     * @param vc                                 Input VC
     * @param model                              GL calculation model
     * @param inheritAttributesFromInputVC       Output VC will contain attributes inherited from input vc
     * @return                                   VC with assigned genotypes
     */
    protected VariantCallContext calculateGenotypes(final FeatureContext features,
                                                    final ReferenceContext refContext,
                                                    final AlignmentContext rawContext,
                                                    Map<String, AlignmentContext> stratifiedContexts,
                                                    final VariantContext vc,
                                                    final GenotypeLikelihoodsCalculationModel model,
                                                    final boolean inheritAttributesFromInputVC,
                                                    final ReadLikelihoods<Allele> likelihoods,
                                                    final SAMFileHeader header) {
        final boolean limitedContext = features == null || refContext == null || rawContext == null || stratifiedContexts == null;
        // if input VC can't be genotyped, exit with either null VCC or, in case where we need to emit all sites, an empty call
        if (hasTooManyAlternativeAlleles(vc) || vc.getNSamples() == 0) {
            return emptyCallContext(features, refContext, rawContext, header);
        }

        final int defaultPloidy = configuration.genotypeArgs.samplePloidy;
        final int maxAltAlleles = configuration.genotypeArgs.MAX_ALTERNATE_ALLELES;


        VariantContext reducedVC = vc;
        if (maxAltAlleles < vc.getAlternateAlleles().size()) {
            final List<Allele> allelesToKeep = AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, defaultPloidy, maxAltAlleles);
            final GenotypesContext reducedGenotypes = allelesToKeep.size() == 1 ? GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy) :
                    AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), allelesToKeep, GenotypeAssignmentMethod.SET_TO_NO_CALL, vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            reducedVC = new VariantContextBuilder(vc).alleles(allelesToKeep).genotypes(reducedGenotypes).make();
        }


        final AFCalculator afCalculator = configuration.genotypeArgs.USE_NEW_AF_CALCULATOR ? newAFCalculator
                : afCalculatorProvider.getInstance(vc,defaultPloidy,maxAltAlleles);
        final AFCalculationResult AFresult = afCalculator.getLog10PNonRef(reducedVC, defaultPloidy, maxAltAlleles, getAlleleFrequencyPriors(vc,defaultPloidy,model));
        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult, vc);

        // posterior probability that at least one alt allele exists in the samples
        final double probOfAtLeastOneAltAllele = Math.pow(10, AFresult.getLog10PosteriorOfAFGT0());

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        final double log10Confidence =
                ! outputAlternativeAlleles.siteIsMonomorphic ||
                        configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES || configuration.annotateAllSitesWithPLs
                        ? AFresult.getLog10PosteriorOfAFEq0() + 0.0
                        : AFresult.getLog10PosteriorOfAFGT0() + 0.0 ;


        // Add 0.0 removes -0.0 occurrences.
        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        // skip this if we are already looking at a vc with NON_REF as the first alt allele i.e. if we are in GenotypeGVCFs
        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic)
                && !forceSiteEmission()
                && noAllelesOrFirstAlleleIsNotNonRef(outputAlternativeAlleles.alleles)) {
            // technically, at this point our confidence in a reference call isn't accurately estimated
            //  because it didn't take into account samples with no data, so let's get a better estimate
            final double[] AFpriors = getAlleleFrequencyPriors(vc, defaultPloidy, model);
            final int INDEX_FOR_AC_EQUALS_1 = 1;
            return limitedContext ? null : estimateReferenceConfidence(vc, stratifiedContexts, AFpriors[INDEX_FOR_AC_EQUALS_1], true, probOfAtLeastOneAltAllele);
        }

        // start constructing the resulting VC
        final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
        final VariantContextBuilder builder = new VariantContextBuilder(callSourceString(), vc.getContig(), vc.getStart(), vc.getEnd(), outputAlleles);

        builder.log10PError(log10Confidence);
        if ( ! passesCallThreshold(phredScaledConfidence) ) {
            builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
        }

        // create the genotypes
        //TODO: omit subsetting if output alleles is not a proper subset of vc.getAlleles
        final GenotypesContext genotypes = outputAlleles.size() == 1 ? GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy) :
                AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), outputAlleles, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));

        // calculating strand bias involves overwriting data structures, so we do it last
        final Map<String, Object> attributes = composeCallAttributes(inheritAttributesFromInputVC, vc, rawContext, stratifiedContexts, features, refContext,
                outputAlternativeAlleles.alternativeAlleleMLECounts(), outputAlternativeAlleles.siteIsMonomorphic, AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()),genotypes,model,likelihoods);

        VariantContext vcCall = builder.genotypes(genotypes).attributes(attributes).make();

        if ( annotationEngine != null && !limitedContext ) { // limitedContext callers need to handle annotations on their own by calling their own annotationEngine
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadPileup pileup = rawContext.getBasePileup();
            stratifiedContexts = AlignmentContext.splitContextBySampleName(pileup, header);

            vcCall = annotationEngine.annotateContext(vcCall, features, refContext, likelihoods, a -> true);
        }

        // if we are subsetting alleles (either because there were too many or because some were not polymorphic)
        // then we may need to trim the alleles (because the original VariantContext may have had to pad at the end).
        if ( outputAlleles.size() != vc.getAlleles().size() && !limitedContext ) // limitedContext callers need to handle allele trimming on their own to keep their alleles in sync
        {
            vcCall = GATKVariantContextUtils.reverseTrimAlleles(vcCall);
        }

        return new VariantCallContext(vcCall, confidentlyCalled(phredScaledConfidence, probOfAtLeastOneAltAllele));
    }

    @VisibleForTesting
    static boolean noAllelesOrFirstAlleleIsNotNonRef(List<Allele> altAlleles) {
        Utils.nonNull(altAlleles);
        return altAlleles.isEmpty() ||  altAlleles.get(0) != (Allele.NON_REF_ALLELE);
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
        private  final List<Allele> alleles;
        private  final boolean siteIsMonomorphic;
        private  final List<Integer> mleCounts;

        private OutputAlleleSubset(final List<Allele> alleles, final List<Integer> mleCounts, final boolean siteIsMonomorphic) {
            Utils.nonNull(alleles, "alleles");
            Utils.nonNull(mleCounts, "mleCounts");
            this.siteIsMonomorphic = siteIsMonomorphic;
            this.alleles = alleles;
            this.mleCounts = mleCounts;
        }

        private List<Allele> outputAlleles(final Allele referenceAllele) {
            return Stream.concat(Stream.of(referenceAllele), alleles.stream()).collect(Collectors.toList());
        }

        public List<Integer> alternativeAlleleMLECounts() {
            return mleCounts;
        }
    }


    /**
     * Provided the exact mode computations it returns the appropriate subset of alleles that progress to genotyping.
     * @param afCalculationResult the allele fraction calculation result.
     * @param vc the variant context
     * @return information about the alternative allele subsetting {@code null}.
     */
    private OutputAlleleSubset calculateOutputAlleleSubset(final AFCalculationResult afCalculationResult, final VariantContext vc) {
        final List<Allele> outputAlleles = new ArrayList<>();
        final List<Integer> mleCounts = new ArrayList<>();
        boolean siteIsMonomorphic = true;
        final List<Allele> alleles = afCalculationResult.getAllelesUsedInGenotyping();
        final int alternativeAlleleCount = alleles.size() - 1;
        int referenceAlleleSize = 0;
        for (final Allele allele : alleles) {
            if (allele.isReference() ) {
                referenceAlleleSize = allele.length();
            } else {
                // we want to keep the NON_REF symbolic allele but only in the absence of a non-symbolic allele, e.g.
                // if we combined a ref / NON_REF gVCF with a ref / alt gVCF
                final boolean isNonRefWhichIsLoneAltAllele = alternativeAlleleCount == 1 && allele.equals(
                        Allele.NON_REF_ALLELE);
                final boolean isPlausible = afCalculationResult.isPolymorphicPhredScaledQual(allele, configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING);

                siteIsMonomorphic &= !isPlausible;

                //it's possible that the upstream deletion that spanned this site was not emitted, mooting the symbolic spanning deletion allele
                final boolean isSpuriousSpanningDeletion = GATKVCFConstants.isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc);
                final boolean toOutput = (isPlausible || forceKeepAllele(allele) || isNonRefWhichIsLoneAltAllele) && !isSpuriousSpanningDeletion;

                if (toOutput) {
                    outputAlleles.add(allele);
                    mleCounts.add(afCalculationResult.getAlleleCountAtMLE(allele));
                    recordDeletion(referenceAlleleSize - allele.length(), vc);
                }
            }
        }

        return new OutputAlleleSubset(outputAlleles,mleCounts,siteIsMonomorphic);
    }

    void clearUpstreamDeletionsLoc() {
        upstreamDeletionsLoc.clear();
    }

    /**
     *  Record deletion to keep
     *  Add deletions to a list.
     *
     * @param deletionSize  size of deletion in bases
     * @param vc            variant context
     */
    void recordDeletion(final int deletionSize, final VariantContext vc) {

        // In a deletion
        if (deletionSize > 0) {
            final SimpleInterval genomeLoc = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart() + deletionSize);
            upstreamDeletionsLoc.add(genomeLoc);
        }
    }

    /**
     * Is the variant context covered by an upstream deletion?
     *
     * @param vc    variant context
     * @return  true if the location is covered by an upstream deletion, false otherwise
     */
    boolean isVcCoveredByDeletion(final VariantContext vc) {
        for (Iterator<SimpleInterval> it = upstreamDeletionsLoc.iterator(); it.hasNext(); ) {
            final SimpleInterval loc = it.next();
            if (!loc.getContig().equals(vc.getContig())) { // deletion is not on contig.
                it.remove();
            } else if (loc.getEnd() < vc.getStart()) { // deletion is before the start.
                it.remove();
            } else if (loc.getStart() == vc.getStart()) {
                // ignore this deletion, the symbolic one does not make reference to it.
            } else { // deletion covers.
                return true;
            }
        }

        return false;
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
     * @param conf the phred scaled quality score
     * @param PofF
     * @return {@code true} iff the variant is confidently called.
     */
    protected final boolean confidentlyCalled(final double conf, final double PofF) {
        return passesCallThreshold(conf)  ||
                (configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES
                        && passesCallThreshold(QualityUtils.phredScaleErrorRate(PofF)));
    }


    /**
     * Checks whether the variant context has too many alternative alleles for progress to genotyping the site.
     * <p>
     *     AF calculation may get into trouble with too many alternative alleles.
     * </p>
     *
     * @param vc the variant context to evaluate.
     *
     * @throws NullPointerException if {@code vc} is {@code null}.
     *
     * @return {@code true} iff there is too many alternative alleles based on
     * {@link GenotypeLikelihoods#MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED}.
     */
    protected final boolean hasTooManyAlternativeAlleles(final VariantContext vc) {
        // protect against too many alternate alleles that we can't even run AF on:
        if (vc.getNAlleles() <= GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED) {
            return false;
        }
        logger.warn("Attempting to genotype more than " + GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED +
                " alleles. Site will be skipped at location "+vc.getContig()+":"+vc.getStart());
        return true;
    }

    /**
     * Produces an empty variant-call context to output when there is no enough data provided to call anything.
     *
     * @param features feature context
     * @param ref the reference context.
     * @param rawContext the read alignment at that location.
     * @return it might be null if no enough information is provided to do even an empty call. For example when
     * we have limited-context (i.e. any of the tracker, reference or alignment is {@code null}.
     */
    protected final VariantCallContext emptyCallContext(final FeatureContext features,
                                                        final ReferenceContext ref,
                                                        final AlignmentContext rawContext,
                                                        final SAMFileHeader header) {
        if (features == null || ref == null || rawContext == null || !forceSiteEmission()) {
            return null;
        }

        VariantContext vc;

        if ( configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            final VariantContext ggaVc = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(features,
                    rawContext.getLocation(), false, configuration.genotypeFilteredAlleles, logger, configuration.alleles);
            if (ggaVc == null) {
                return null;
            }
            vc = new VariantContextBuilder(callSourceString(), ref.getInterval().getContig(), ggaVc.getStart(),
                    ggaVc.getEnd(), ggaVc.getAlleles()).make();
        } else {
            // deal with bad/non-standard reference bases
            if ( !Allele.acceptableAlleleBases(new byte[]{ref.getBase()}) ) {
                return null;
            }
            final Set<Allele> alleles = new LinkedHashSet<>(Collections.singleton(Allele.create(ref.getBase(),true)));
            vc = new VariantContextBuilder(callSourceString(), ref.getInterval().getContig(),
                    ref.getInterval().getStart(), ref.getInterval().getStart(), alleles).make();
        }

        if ( vc != null && annotationEngine != null ) {
            // Note: we want to use the *unfiltered* and *unBAQed* context for the annotations
            final ReadPileup pileup = rawContext.getBasePileup();
            vc = annotationEngine.annotateContext(vc, features, ref, null, a -> true);
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

    protected final VariantCallContext estimateReferenceConfidence(final VariantContext vc, final Map<String, AlignmentContext> contexts,
                                                                   final double log10OfTheta, final boolean ignoreCoveredSamples, final double initialPofRef) {
        if ( contexts == null ) {
            return null;
        }

        // add contribution from each sample that we haven't examined yet i.e. those with null contexts
        final double log10POfRef = Math.log10(initialPofRef) + contexts.values().stream()
                .filter(context ->  context == null || !ignoreCoveredSamples )
                .mapToInt(context -> context == null ? 0 : context.getBasePileup().size())  //get the depth
                .mapToDouble(depth -> estimateLog10ReferenceConfidenceForOneSample(depth, log10OfTheta))
                .sum();

        return new VariantCallContext(vc, passesCallThreshold(QualityUtils.phredScaleLog10CorrectRate(log10POfRef)), false);
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
    protected final double[] getAlleleFrequencyPriors( final VariantContext vc, final int defaultPloidy, final GenotypeLikelihoodsCalculationModel model ) {
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
    protected final double estimateLog10ReferenceConfidenceForOneSample(final int depth, final double log10OfTheta) {
        Utils.validateArg(depth >= 0, "depth may not be negative");
        final double log10PofNonRef = log10OfTheta + MathUtils.log10BinomialProbability(depth, 0);
        return MathUtils.log10OneMinusX(Math.pow(10.0, log10PofNonRef));
    }

    protected final boolean passesEmitThreshold(final double conf, final boolean bestGuessIsRef) {
        return (configuration.outputMode == OutputMode.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) &&
                passesCallThreshold(conf);
    }

    protected final boolean passesCallThreshold(final double conf) {
        return conf >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    protected Map<String,Object> composeCallAttributes(final boolean inheritAttributesFromInputVC, final VariantContext vc,
                                                       final AlignmentContext rawContext, final Map<String, AlignmentContext> stratifiedContexts, final FeatureContext tracker, final ReferenceContext refContext, final List<Integer> alleleCountsofMLE, final boolean bestGuessIsRef,
                                                       final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes,
                                                       final GenotypeLikelihoodsCalculationModel model, final ReadLikelihoods<Allele> likelihoods) {
        final Map<String, Object> attributes = new LinkedHashMap<>();

        final boolean limitedContext = tracker == null || refContext == null || rawContext == null || stratifiedContexts == null;

        // inherit attributes from input vc if requested
        if (inheritAttributesFromInputVC) {
            attributes.putAll(vc.getAttributes());
        }
        // if the site was down-sampled, record that fact
        if ( !limitedContext && rawContext.hasPileupBeenDownsampled() ) {
            attributes.put(GATKVCFConstants.DOWNSAMPLED_KEY, true);
        }

        // add the MLE AC and AF annotations
        if (!alleleCountsofMLE.isEmpty()) {
            attributes.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
            final List<Double> MLEfrequencies = calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes);
            attributes.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
        }

        if (doAlleleSpecificCalcs){
            List<Double> perAlleleQuals = new ArrayList<>();
            //Per-allele quals are not calculated for biallelic sites
            if (AFresult.getAllelesUsedInGenotyping().size() > 2) {
                for (final Allele a : allAllelesToUse) {
                    if (a.isNonReference()) {
                        perAlleleQuals.add(AFresult.getLog10PosteriorOfAFEq0ForAllele(a));
                    }
                }
            }
            else {
                perAlleleQuals.add(AFresult.getLog10PosteriorOfAFEq0());
            }

            attributes.put(GATKVCFConstants.AS_QUAL_KEY, perAlleleQuals);
        }

        if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
            attributes.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, vc.getAlternateAlleles().size());
        }

        return attributes;
    }

    private List<Double> calculateMLEAlleleFrequencies(final List<Integer> alleleCountsofMLE, final GenotypesContext genotypes) {
        final long AN = genotypes.stream().flatMap(g -> g.getAlleles().stream()).filter(Allele::isCalled).count();
        return alleleCountsofMLE.stream().map(AC -> Math.min(1.0, (double) AC / AN)).collect(Collectors.toList());
    }

    /**
     * Calculates the active state profile value for a single sample.
     *
     * @param log10GenotypeLikelihoods the single sample genotype likelihoods.
     * @return log10 probability from 0 to -Infinity.
     */
    public double calculateSingleSampleRefVsAnyActiveStateProfileValue(final double[] log10GenotypeLikelihoods) {
        Utils.nonNull(log10GenotypeLikelihoods, "the input likelihoods cannot be null");
        Utils.validateArg(log10GenotypeLikelihoods.length == this.configuration.genotypeArgs.samplePloidy + 1, "wrong likelihoods dimensions");

        final double[] log10Priors = log10AlleleFrequencyPriorsSNPs.forTotalPloidy(this.configuration.genotypeArgs.samplePloidy);
        final double log10ACeq0Likelihood = log10GenotypeLikelihoods[0];
        final double log10ACeq0Prior = log10Priors[0];
        final double log10ACeq0Posterior = log10ACeq0Likelihood + log10ACeq0Prior;

        // If the maximum a posteriori AC is 0 then the profile value must be 0.0 as per existing code; it does
        // not matter whether AC > 0 is at all plausible.
        boolean mapACeq0 = true;
        for (int AC = 1; AC < log10Priors.length; AC++) {
            if (log10Priors[AC] + log10GenotypeLikelihoods[AC] > log10ACeq0Posterior) {
                mapACeq0 = false;
                break;
            }
        }
        if (mapACeq0) {
            return 0.0;
        }

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
        if (normalizedLog10ACeq0Posterior >= QualityUtils.qualToErrorProbLog10(configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING)) {
            return 0.0;
        }

        return 1.0 - Math.pow(10.0, normalizedLog10ACeq0Posterior);

    }
}
