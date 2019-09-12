package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
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

    protected final AlleleFrequencyCalculator alleleFrequencyCalculator;

    protected final Config configuration;

    protected VariantAnnotatorEngine annotationEngine;

    protected Logger logger;

    protected final int numberOfGenomes;

    protected final SampleList samples;

    private final AFPriorProvider log10AlleleFrequencyPriorsSNPs;

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
                               final boolean doAlleleSpecificCalcs) {
        this.configuration = Utils.nonNull(configuration, "the configuration cannot be null");
        this.samples = Utils.nonNull(samples, "the sample list cannot be null");
        this.doAlleleSpecificCalcs = doAlleleSpecificCalcs;
        logger = LogManager.getLogger(getClass());
        numberOfGenomes = this.samples.numberOfSamples() * configuration.genotypeArgs.samplePloidy;
        log10AlleleFrequencyPriorsSNPs = composeAlleleFrequencyPriorProvider(numberOfGenomes,
                configuration.genotypeArgs.snpHeterozygosity, configuration.genotypeArgs.inputPrior);
        alleleFrequencyCalculator = AlleleFrequencyCalculator.makeCalculator(configuration.genotypeArgs);
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
     * Main entry function to calculate genotypes of a given VC with corresponding GL's that is shared across genotypers (namely UG and HC).
     *
     * Completes a variant context with genotype calls and associated annotations given the genotype likelihoods and
     * the model that need to be applied.
     *
     * @param vc                                 Input variant context to complete.
     * @param model                              GL calculation model
     * @return                                   VC with assigned genotypes
     */
    public VariantContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel model, final List<VariantContext> givenAlleles) {
        // if input VC can't be genotyped, exit with either null VCC or, in case where we need to emit all sites, an empty call
        if (hasTooManyAlternativeAlleles(vc) || vc.getNSamples() == 0) {
            return null;
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


        final AFCalculationResult AFresult = alleleFrequencyCalculator.calculate(reducedVC, defaultPloidy);
        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult, vc, givenAlleles);

        // posterior probability that at least one alt allele exists in the samples
        final double probOfAtLeastOneAltAllele = Math.pow(10, AFresult.log10ProbVariantPresent());

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        final double log10Confidence =
                ! outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                        ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0 ;


        // Add 0.0 removes -0.0 occurrences.
        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        // skip this if we are already looking at a vc with NON_REF as the first alt allele i.e. if we are in GenotypeGVCFs
        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic) && !emitAllSites()
                && noAllelesOrFirstAlleleIsNotNonRef(outputAlternativeAlleles.alleles) && givenAlleles.isEmpty()) {
            return null;
        }

        // return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
        if (! emitAllSites() && outputAlternativeAlleles.alleles.size() == 1 && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0))) {
            return null;
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
        final Map<String, Object> attributes = composeCallAttributes(vc, outputAlternativeAlleles.alternativeAlleleMLECounts(),
                AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()),genotypes);

        return builder.genotypes(genotypes).attributes(attributes).make();
    }

    public VariantContext calculateGenotypes(final VariantContext vc, final GenotypeLikelihoodsCalculationModel model) {
        return calculateGenotypes(vc, model, Collections.emptyList());
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
    private OutputAlleleSubset calculateOutputAlleleSubset(final AFCalculationResult afCalculationResult, final VariantContext vc, final List<VariantContext> givenAlleles) {
        final List<Allele> outputAlleles = new ArrayList<>();
        final List<Integer> mleCounts = new ArrayList<>();
        boolean siteIsMonomorphic = true;
        final List<Allele> alleles = afCalculationResult.getAllelesUsedInGenotyping();
        final int alternativeAlleleCount = alleles.size() - 1;
        int referenceAlleleSize = 0;

        final Set<Allele> forcedAlleles = AssemblyBasedCallerUtils.getAllelesConsistentWithGivenAlleles(givenAlleles, vc);

        for (final Allele allele : alleles) {
            if (allele.isReference() ) {
                referenceAlleleSize = allele.length();
            } else {
                // we want to keep the NON_REF symbolic allele but only in the absence of a non-symbolic allele, e.g.
                // if we combined a ref / NON_REF gVCF with a ref / alt gVCF
                final boolean isNonRefWhichIsLoneAltAllele = alternativeAlleleCount == 1 && allele.equals(Allele.NON_REF_ALLELE);
                final boolean isPlausible = afCalculationResult.passesThreshold(allele, configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING);

                //it's possible that the upstream deletion that spanned this site was not emitted, mooting the symbolic spanning deletion allele
                final boolean isSpuriousSpanningDeletion = GATKVCFConstants.isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc);

                final boolean toOutput = (isPlausible || forceKeepAllele(allele) || isNonRefWhichIsLoneAltAllele || forcedAlleles.contains(allele) ) && !isSpuriousSpanningDeletion;

                siteIsMonomorphic &= !(isPlausible && !isSpuriousSpanningDeletion);

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

    protected boolean emitAllSites() {
        return configuration.outputMode == OutputMode.EMIT_ALL_SITES;
    }

    protected final boolean passesEmitThreshold(final double conf, final boolean bestGuessIsRef) {
        return (configuration.outputMode == OutputMode.EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) &&
                passesCallThreshold(conf);
    }

    protected final boolean passesCallThreshold(final double conf) {
        return conf >= configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING;
    }

    protected Map<String,Object> composeCallAttributes(final VariantContext vc, final List<Integer> alleleCountsofMLE,
                                                       final AFCalculationResult AFresult, final List<Allele> allAllelesToUse, final GenotypesContext genotypes) {
        final Map<String, Object> attributes = new LinkedHashMap<>();

        // add the MLE AC and AF annotations
        if (!alleleCountsofMLE.isEmpty()) {
            attributes.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, alleleCountsofMLE);
            final List<Double> MLEfrequencies = calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes);
            attributes.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, MLEfrequencies);
        }

        if (doAlleleSpecificCalcs){
            List<Integer> perAlleleQuals = new ArrayList<>();
            //Per-allele quals are not calculated for biallelic sites
            if (AFresult.getAllelesUsedInGenotyping().size() > 2) {
                for (final Allele a : allAllelesToUse) {
                    if (a.isNonReference()) {
                        perAlleleQuals.add((int)Math.round(AFresult.getLog10PosteriorOfAlleleAbsent(a)*-10));
                    }
                }
            }
            else {
                perAlleleQuals.add((int)Math.round(AFresult.log10ProbOnlyRefAlleleExists()*-10));
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
        // This is another condition to return a 0.0 also present in AlleleFrequencyCalculator code as well.
        if (normalizedLog10ACeq0Posterior >= QualityUtils.qualToErrorProbLog10(configuration.genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING)) {
            return 0.0;
        }

        return 1.0 - Math.pow(10.0, normalizedLog10ACeq0Posterior);

    }
}
