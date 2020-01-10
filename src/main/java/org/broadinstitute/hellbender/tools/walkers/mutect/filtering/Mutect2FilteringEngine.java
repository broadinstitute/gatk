package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import java.nio.file.Path;

import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class Mutect2FilteringEngine {
    public static final double EPSILON = 1.0e-10;

    public static final double MIN_REPORTABLE_ERROR_PROBABILITY = 0.1;

    private final List<Mutect2VariantFilter> filters = new ArrayList<>();
    private final List<Mutect2AlleleFilter<?>> alleleFilters = new ArrayList<>();
    private final Set<String> normalSamples;

    public static final List<String> STANDARD_MUTECT_INFO_FIELDS_FOR_FILTERING = Arrays.asList(
            GATKVCFConstants.MEDIAN_MAPPING_QUALITY_KEY,
            GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY,
            GATKVCFConstants.MEDIAN_READ_POSITON_KEY,
            GATKVCFConstants.MEDIAN_FRAGMENT_LENGTH_KEY
    );

    /**
     * DATA ACCUMULATED AND LEARNED ON EACH PASS OF {@link FilterMutectCalls}
     */
    private final ThresholdCalculator thresholdCalculator;
    private final FilteringOutputStats filteringOutputStats;
    private final SomaticClusteringModel somaticClusteringModel;

    public Mutect2FilteringEngine(M2FiltersArgumentCollection MTFAC, final VCFHeader vcfHeader, final File mutectStatsTable) {
        thresholdCalculator = new ThresholdCalculator(MTFAC.thresholdStrategy, MTFAC.initialPosteriorThreshold, MTFAC.maxFalsePositiveRate, MTFAC.fScoreBeta);

        somaticClusteringModel = new SomaticClusteringModel(MTFAC, mutectStatsTable.exists() ? MutectStats.readFromFile(mutectStatsTable) : Collections.emptyList());

        normalSamples = vcfHeader.getMetaDataInInputOrder().stream()
                .filter(line -> line.getKey().equals(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER))
                .map(VCFHeaderLine::getValue)
                .collect(Collectors.toSet());

        buildFiltersList(MTFAC);
        filteringOutputStats = new FilteringOutputStats(filters);
    }

    //THE FOLLOWING ARE HELPER METHODS FOR FILTERS THAT IMPLEMENT {@link Mutect2VariantFilter}
    public boolean isNormal(final Genotype genotype) { return normalSamples.contains(genotype.getSampleName()); }

    public boolean isTumor(final Genotype genotype) { return !isNormal(genotype); }

    /**
     * Maximum probability that a potential variant is not a true somatic mutation.  Variants with error probabilities
     * at or below this threshold are called; variants with error probabilities above are filtered.
     */
    public double getThreshold() { return thresholdCalculator.getThreshold(); }

    public SomaticClusteringModel getSomaticClusteringModel() { return somaticClusteringModel; }

    public double posteriorProbabilityOfError(final VariantContext vc, final double logOddsOfRealVersusError, final int altIndex) {
        return posteriorProbabilityOfError(logOddsOfRealVersusError, getLogSomaticPrior(vc, altIndex));
    }

    public double posteriorProbabilityOfNormalArtifact(final double negativeLogOddsOfNormalArtifact) {
        return posteriorProbabilityOfError(negativeLogOddsOfNormalArtifact, somaticClusteringModel.getLogPriorOfVariantVersusArtifact());
    }

    public double getLogSomaticPrior(final VariantContext vc, final int altIndex) {
        return somaticClusteringModel.getLogPriorOfSomaticVariant(vc, altIndex);
    }

    public static double posteriorProbabilityOfError(final double logOddsOfRealVersusError, final double logPriorOfReal) {
        final double[] unweightedPosteriorOfRealAndError = new double[] {logOddsOfRealVersusError + logPriorOfReal,
                NaturalLogUtils.log1mexp(logPriorOfReal)};

        final double[] posteriorOfRealAndError = NaturalLogUtils.normalizeFromLogToLinearSpace(unweightedPosteriorOfRealAndError);

        return posteriorOfRealAndError[1];
    }

    public int[] sumADsOverSamples(final VariantContext vc, final boolean includeTumor, final boolean includeNormal) {
        final int[] ADs = new int[vc.getNAlleles()];
        vc.getGenotypes().stream().filter(g -> (includeTumor && isTumor(g)) || (includeNormal && isNormal(g)))
                .map(Genotype::getAD).forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
        return ADs;
    }

    public int[] sumStrandCountsOverSamples(final VariantContext vc, final boolean includeTumor, final boolean includeNormal) {
        final int[] result = new int[4];
        vc.getGenotypes().stream().filter(g -> (includeTumor && isTumor(g)) || (includeNormal && isNormal(g)))
                .filter(g -> g.hasExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY))
                .map(g -> StrandBiasTest.getStrandCounts(g)).forEach(sbbs -> new IndexRange(0, 4).forEach(n -> result[n] += sbbs[n]));
        return result;
    }

    public double[] weightedAverageOfTumorAFs(final VariantContext vc) {
        final MutableDouble totalWeight = new MutableDouble(0);
        final double[] AFs = new double[vc.getNAlleles() - 1];
        vc.getGenotypes().stream().filter(this::isTumor).forEach(g ->  {
            final double weight = MathUtils.sum(g.getAD());
            totalWeight.add(weight);
            final double[] sampleAFs = VariantContextGetters.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                    () -> new double[] {0.0}, 0.0);
            MathArrays.scaleInPlace(weight, sampleAFs);
            MathUtils.addToArrayInPlace(AFs, sampleAFs);
        });
        MathArrays.scaleInPlace(1/totalWeight.getValue(), AFs);
        return AFs;
    }

    /**
     * Get the (Natural) log odds of variants from the log-10 odds in a VariantContext object
     */
    public static double[] getTumorLogOdds(final VariantContext vc) {
        final double[] tumorLog10Odds = VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY);
        return tumorLog10Odds == null ? null : MathUtils.applyToArrayInPlace(tumorLog10Odds, MathUtils::log10ToLog);
    }
    // END HELPER METHODS

    /**
     * record data from a potential variant in a non-final pass of {@link FilterMutectCalls}
     */
    public void accumulateData(final VariantContext vc, final ReferenceContext referenceContext) {
        // ignore GVCF mode sites where the only alt is NON-REF
        if (vc.getAlleles().stream().noneMatch(a -> a.isNonReference() && !a.isNonRefAllele())) {
            return;
        }

        final ErrorProbabilities errorProbabilities = new ErrorProbabilities(filters, alleleFilters, vc, this, referenceContext);
        filters.forEach(f -> f.accumulateDataForLearning(vc, errorProbabilities, this));
        alleleFilters.forEach(f -> f.accumulateDataForLearning(vc, errorProbabilities, this));
        final int[] tumorADs = sumADsOverSamples(vc, true, false);
        final double[] tumorLogOdds = Mutect2FilteringEngine.getTumorLogOdds(vc);

        somaticClusteringModel.record(tumorADs, tumorLogOdds, errorProbabilities.getTechnicalArtifactProbability(),
                errorProbabilities.getNonSomaticProbability(), vc);
        thresholdCalculator.addArtifactProbability(errorProbabilities.getErrorProbability());
    }

    /**
     * Refine model parameters based on data acquired in a non-final pass of {@link FilterMutectCalls}
     */
    public void learnParameters() {
        filters.forEach(Mutect2VariantFilter::learnParametersAndClearAccumulatedData);
        somaticClusteringModel.learnAndClearAccumulatedData();
        thresholdCalculator.relearnThresholdAndClearAcumulatedProbabilities();

        filteringOutputStats.clear();
    }

    public void learnThreshold() {
        thresholdCalculator.relearnThresholdAndClearAcumulatedProbabilities();
        filteringOutputStats.clear();
    }

    /**
     * Create a filtered variant and record statistics for the final pass of {@link FilterMutectCalls}
     */
    public VariantContext applyFiltersAndAccumulateOutputStats(final VariantContext vc, final ReferenceContext referenceContext) {
        final VariantContextBuilder vcb = new VariantContextBuilder(vc).filters(new HashSet<>());

        final ErrorProbabilities errorProbabilities = new ErrorProbabilities(filters, alleleFilters, vc, this, referenceContext);
        filteringOutputStats.recordCall(errorProbabilities, getThreshold() - EPSILON);

        final boolean variantFailsFilters = errorProbabilities.getErrorProbability() > Math.min(1 - EPSILON, Math.max(EPSILON, getThreshold()));
        final double maxErrorProb = errorProbabilities.getProbabilitiesByFilter().values().stream().mapToDouble(p->p).max().orElse(1);

        for (final Map.Entry<Mutect2VariantFilter, Double> entry : errorProbabilities.getProbabilitiesByFilter().entrySet()) {
            final double errorProbability = entry.getValue();

            entry.getKey().phredScaledPosteriorAnnotationName().ifPresent(annotation -> {
                if (entry.getKey().requiredAnnotations().stream().allMatch(vc::hasAttribute)) {
                    vcb.attribute(annotation, QualityUtils.errorProbToQual(errorProbability));
                }
            });

            // error probability must exceed threshold, and just in case threshold is bad, probabilities close to 1 must be filtered
            // and probabilities close to 0 must not be filtered
            if (variantFailsFilters && errorProbability >= Math.min(maxErrorProb, MIN_REPORTABLE_ERROR_PROBABILITY)) {
                vcb.filter(entry.getKey().filterName());
            }
        }

        // apply allele specific filters
        List<String> siteFilters = new ArrayList<>();
        List<Iterator<String>> ASFilters =
                errorProbabilities.getProbabilitiesByFilterAndAllele().entrySet().stream().filter(entry -> !entry.getValue().isEmpty()).map(
                entry -> addFilterStrings(entry.getValue(), siteFilters, entry.getKey().filterName())).collect(Collectors.toList());

        siteFilters.forEach(vcb::filter);
        List<String> orderedASFilterStrings = vc.getAlleles().stream().map(allele -> allele.isReference() || allele.isSymbolic() ?
            VCFConstants.EMPTY_INFO_FIELD : getMergedFilterStringForAllele(ASFilters)).collect(Collectors.toList());
        String finalAttrString = AnnotationUtils.encodeAnyASList(orderedASFilterStrings);

        vcb.putAttributes(Collections.singletonMap(GATKVCFConstants.AS_FILTER_STATUS_KEY, finalAttrString));
        return vcb.make();
    }

    /**
     * Creates a comma separated string of all the filters that apply to the allele. This is basically
     * a pivot of the data. we have filterlist -> allele -> filterName. and we want allele -> list of filterName
     * @param alleleSpecificFilters all of the allele specific filters with the allele filter info
     * @return encoded (comma separated) list of filters that apply to the allele
     */
    private String getMergedFilterStringForAllele(List<Iterator<String>> alleleSpecificFilters) {
        // loop through each filter and pull out the filters the specified allele
        List<String> results = alleleSpecificFilters.stream().map(alleleValuesIterator -> alleleValuesIterator.next()).distinct().collect(Collectors.toList());
        if (results.size() > 1 && results.contains(VCFConstants.PASSES_FILTERS_v4)) {
            results.remove(VCFConstants.PASSES_FILTERS_v4);
        } else if (results.isEmpty()) {
            results.add(VCFConstants.PASSES_FILTERS_v4);
        }
        return AnnotationUtils.encodeStringList(results);
    }

    /**
     * For each allele, determine whether the filter should be applied. also determine if the filter should apply to the site
     * @param probabilities the probability computed by the filter for the allele
     * @param siteFilters output value - filter name is added if it should apply to the site
     * @param filterName the name of the filter used in the vcf
     * @return Iterator of filters for an allele
     */
    private Iterator<String> addFilterStrings(List<Double> probabilities, List<String> siteFilters, String filterName) {
        double thresh = getThreshold();
        double min = Math.min(1 - EPSILON, Math.max(EPSILON, thresh));
        List<String> results = probabilities.stream().map(value -> value > Math.min(1 - EPSILON, Math.max(EPSILON, getThreshold())) ?
                        filterName : VCFConstants.PASSES_FILTERS_v4).collect(Collectors.toList());
//        List<String> realFilters = results.stream().filter(x -> !x.equals(VCFConstants.EMPTY_INFO_FIELD)).collect(Collectors.toList());
        if (!results.isEmpty() && results.stream().allMatch(x -> x.equals(filterName))) {
            siteFilters.add(filterName);
        }
        return results.iterator();
    }

    public static double roundFinitePrecisionErrors(final double probability) {
        return Math.max(Math.min(probability, 1.0), 0.0);
    }

    /**
     * Write statistics collected in the final pass of {@link FilterMutectCalls}
     * @param filteringStats where to write the statistics.
     */
    public void writeFilteringStats(final Path filteringStats) {
        filteringOutputStats.writeFilteringStats(filteringStats, getThreshold(), somaticClusteringModel.clusteringMetadata());
    }

    private void buildFiltersList(final M2FiltersArgumentCollection MTFAC) {
        alleleFilters.add(new TumorEvidenceFilter());
        alleleFilters.add(new BaseQualityFilter(MTFAC.minMedianBaseQuality));
        alleleFilters.add(new MappingQualityFilter(MTFAC.minMedianMappingQuality, MTFAC.longIndelLength));
        alleleFilters.add(new DuplicatedAltReadFilter(MTFAC.uniqueAltReadCount));
        filters.add(new StrandArtifactFilter());
        alleleFilters.add(new ContaminationFilter(MTFAC.contaminationTables, MTFAC.contaminationEstimate));
        filters.add(new StrictStrandBiasFilter(MTFAC.minReadsOnEachStrand));
        alleleFilters.add(new ReadPositionFilter(MTFAC.minMedianReadPosition));
        alleleFilters.add(new MinAlleleFractionFilter(MTFAC.minAf));

        // convert to allele specific later
        // Normal Artifact Filter doesn't apply to mitochondria because we are not comparing
        // tumor and normal
        filters.add(new NormalArtifactFilter(MTFAC.normalPileupPValueThreshold));
        filters.add(new NRatioFilter(MTFAC.nRatio));

        // filters that don't apply to specific alleles
        filters.add(new PanelOfNormalsFilter());

        if (!MTFAC.readOrientationPriorTarGzs.isEmpty()) {
            final List<File> artifactTables = MTFAC.readOrientationPriorTarGzs.stream().flatMap(tarGz -> {
                final File extractDir = IOUtils.createTempDir("extract");
                IOUtils.extractTarGz(tarGz.toPath(), extractDir.toPath());
                return Arrays.stream(extractDir.listFiles());
            }).collect(Collectors.toList());

            filters.add(new ReadOrientationFilter(artifactTables));
        }

        if (MTFAC.mitochondria) {
            filters.add(new ChimericOriginalAlignmentFilter(MTFAC.maxNuMTFraction));
            alleleFilters.add(new NuMTFilter(MTFAC.medianAutosomalCoverage, MTFAC.maxNuMTAutosomalCopies));
        } else {
            filters.add(new ClusteredEventsFilter(MTFAC.maxEventsInRegion));
            filters.add(new MultiallelicFilter(MTFAC.numAltAllelesThreshold));
            filters.add(new FragmentLengthFilter(MTFAC.maxMedianFragmentLengthDifference));
            filters.add(new PolymeraseSlippageFilter(MTFAC.minSlippageLength, MTFAC.slippageRate));
            filters.add(new FilteredHaplotypeFilter(MTFAC.maxDistanceToFilteredCallOnSameHaplotype));
            filters.add(new GermlineFilter(MTFAC.tumorSegmentationTables));
        }
    }
}
