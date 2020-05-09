package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.bdgenomics.adam.util.PhredUtils;
import org.broadinstitute.hellbender.tools.haplotypecaller.GenotypePriorCalculator;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.lang.reflect.Array;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Base class for genotyper engines.
 */
public abstract class GenotypingEngine<Config extends StandardCallerArgumentCollection> {

    protected final VariationalAlleleFrequencyCalculator alleleFrequencyCalculator;

    protected final Config configuration;

    protected VariantAnnotatorEngine annotationEngine;

    protected Logger logger;

    protected final int numberOfGenomes;

    protected final SampleList samples;

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
        alleleFrequencyCalculator = VariationalAlleleFrequencyCalculator.makeCalculator(configuration.genotypeArgs);
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
     * @return                                   VC with assigned genotypes
     */
    public VariantContext calculateGenotypes(final VariantContext vc, final GenotypePriorCalculator gpc, final List<VariantContext> givenAlleles) {
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
                    AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), allelesToKeep, gpc, GenotypeAssignmentMethod.SET_TO_NO_CALL, vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));
            reducedVC = new VariantContextBuilder(vc).alleles(allelesToKeep).genotypes(reducedGenotypes).make();
        }


        final AFCalculationResult AFresult = alleleFrequencyCalculator.calculate(reducedVC, defaultPloidy);
        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult, vc, givenAlleles);

        // note the math.abs is necessary because -10 * 0.0 => -0.0 which isn't nice
        final double log10Confidence =
                    !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                            ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0;

        // Add 0.0 removes -0.0 occurrences.
        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        // skip this if we are already looking at a vc with NON_REF as the first alt allele i.e. if we are in GenotypeGVCFs
        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic) && !emitAllActiveSites()
                && noAllelesOrFirstAlleleIsNotNonRef(outputAlternativeAlleles.alleles) && givenAlleles.isEmpty()) {
            return null;
        }

        // return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
        if (! emitAllActiveSites() && outputAlternativeAlleles.alleles.size() == 1 && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0))) {
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
                AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), outputAlleles, gpc, configuration.genotypeArgs.genotypeAssignmentMethod, vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));

        if (configuration.genotypeArgs.usePosteriorProbabilitiesToCalculateQual && hasPosteriors(genotypes)) {
            final double log10NoVariantPosterior = nonVariantPresentLog10PosteriorProbability(genotypes) * -.1;
            final double qualUpdate = !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                    ? log10NoVariantPosterior + 0.0 : MathUtils.log10OneMinusPow10(log10NoVariantPosterior) + 0.0;
            if (!Double.isNaN(qualUpdate)) {
                builder.log10PError(qualUpdate);
            }
        }

        // calculating strand bias involves overwriting data structures, so we do it last
        final Map<String, Object> attributes = composeCallAttributes(vc, outputAlternativeAlleles.alternativeAlleleMLECounts(),
                AFresult, outputAlternativeAlleles.outputAlleles(vc.getReference()),genotypes);

        return builder.genotypes(genotypes).attributes(attributes).make();
    }

    protected double nonVariantPresentLog10PosteriorProbability(final GenotypesContext gc) {
        return gc.stream()
                .map(gt -> gt.getExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY))
                .mapToDouble(v -> coherceToDouble(v, Double.NaN, true))
                .filter(v -> !Double.isNaN(v))
                .min().orElse(Double.NaN);
    }

    private double coherceToDouble(final Object obj, final double defaultValue, final boolean takeFirstElement) {
        if (obj == null) {
            return defaultValue;
        } else if (obj instanceof CharSequence) {
            try {
                return Double.parseDouble(obj.toString());
            } catch (final NumberFormatException ex) {
                return defaultValue;
            }
        } else if (obj instanceof Number) {
            return ((Number) obj).doubleValue();
        } else if (takeFirstElement) {
            if (obj instanceof Collection) {
                if( ((Collection)obj).isEmpty()) {
                    return defaultValue;
                } else if (obj instanceof List) {
                    final List<?> asList = (List<?>) obj;
                    return coherceToDouble(asList.get(0), defaultValue, false);
                } else {
                    final Collection<?> collection = (Collection<?>) obj;
                    return coherceToDouble(collection.iterator().next(), defaultValue, false);
                }
            } else if (obj.getClass().isArray()) {
                if (obj.getClass().getComponentType().isPrimitive()) {
                    if (Array.getLength(obj) == 0) {
                        return defaultValue;
                    } else {
                        return coherceToDouble(Array.get(obj, 1), defaultValue, false);
                    }
                } else {
                    final Object[] array = (Object[]) obj;
                    return array.length != 0 ? coherceToDouble(array[0], defaultValue, false) : defaultValue;
                }
            } else {
                return defaultValue;
            }
        } else {
            return defaultValue;
        }
    }

    private boolean hasPosteriors(final GenotypesContext gc) {
        for (final Genotype genotype : gc) {
            if (genotype.hasExtendedAttribute(VCFConstants.GENOTYPE_POSTERIORS_KEY)) {
                return true;
            }
        }
        return false;
    }

    public VariantContext calculateGenotypes(final VariantContext vc) {
        return calculateGenotypes(vc, null, Collections.emptyList());
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

    /**
     * Whether or not all calls at any region over the activity threshold should be emitted regardless of confidence.
     * This does not necessarily emit calls for all sites in a region (see {@link OutputMode}).
     */
    protected boolean emitAllActiveSites() {
        return configuration.outputMode == OutputMode.EMIT_ALL_ACTIVE_SITES;
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
                        //*-10 to convert from log10-scale to Phred-scale, as QUALs are typically represented
                        perAlleleQuals.add((int)Math.round(AFresult.getLog10PosteriorOfAlleleAbsent(a)*-10));
                    }
                }
            }
            else {
                //*-10 to convert from log10-scale to Phred-scale, as QUALs are typically represented
                perAlleleQuals.add((int)Math.round(AFresult.log10ProbOnlyRefAlleleExists()*-10));
            }

            attributes.put(GATKVCFConstants.AS_QUAL_KEY, perAlleleQuals.stream().mapToInt(q -> Math.round(q)).boxed().collect(Collectors.toList()));
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
        return alleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior(log10GenotypeLikelihoods, true);
    }
}
