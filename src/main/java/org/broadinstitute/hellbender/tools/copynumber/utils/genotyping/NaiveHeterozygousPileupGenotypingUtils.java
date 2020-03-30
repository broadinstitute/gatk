package org.broadinstitute.hellbender.tools.copynumber.utils.genotyping;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multiset;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Naive methods for binomial genotyping of heterozygous sites from pileup allele counts.
 * Filters for total count and overlap with copy-ratio intervals are also implemented.
 */
public final class NaiveHeterozygousPileupGenotypingUtils {
    private static final Logger logger = LogManager.getLogger(NaiveHeterozygousPileupGenotypingUtils.class);

    private NaiveHeterozygousPileupGenotypingUtils() {
    }

    public static final class NaiveHeterozygousPileupGenotypingResult {
        private final ImmutableList<AllelicCountCollection> hetAllelicCountsPerSample;
        private final AllelicCountCollection hetNormalAllelicCounts;

        /**
         * @param hetNormalAllelicCounts    {@code null}, if result generated in case-only mode
         */
        private NaiveHeterozygousPileupGenotypingResult(final List<AllelicCountCollection> hetAllelicCountsPerSample,
                                                        final AllelicCountCollection hetNormalAllelicCounts) {
            Utils.nonEmpty(hetAllelicCountsPerSample);
            Utils.validateArg((int) Stream.of(hetAllelicCountsPerSample, Collections.singletonList(hetNormalAllelicCounts))
                            .flatMap(Collection::stream)
                            .filter(Objects::nonNull)
                            .map(AllelicCountCollection::getIntervals)
                            .distinct()
                            .count() == 1,
                    "Allelic-count sites must be identical across all samples.");
            CopyNumberArgumentValidationUtils.getValidatedSequenceDictionary(
                    Stream.of(hetAllelicCountsPerSample, Collections.singletonList(hetNormalAllelicCounts))
                            .flatMap(Collection::stream)
                            .toArray(AllelicCountCollection[]::new));
            this.hetAllelicCountsPerSample = ImmutableList.copyOf(hetAllelicCountsPerSample);
            this.hetNormalAllelicCounts = hetNormalAllelicCounts;
        }

        public List<AllelicCountCollection> getHetAllelicCountsPerSample() {
            return hetAllelicCountsPerSample;
        }

        /**
         * @return {@code null}, if result generated in case-only mode
         */
        public AllelicCountCollection getHetNormalAllelicCounts() {
            return hetNormalAllelicCounts;
        }
    }

    /**
     * Filters allelic counts based on total count, overlap with copy-ratio intervals, and a naive test for heterozygosity.
     * This method can be called either in matched-normal or case-only mode.
     * In the former, the test for heterozygosity is only performed on the normal and the corresponding sites are filtered
     * out of the case samples; this ensures that the set of genotyped sites is identical across all samples.
     * In the latter (activated when {@code normalAllelicCounts} is {@code null}), the test for heterozygosity is performed
     * individually on each case sample, and then the intersection of all genotyped sites is taken to ensure an identical
     * set across all samples.
     * Validation of allelic-count sites and sequence dictionaries will be performed.
     * @param allelicCountsPerSample    non-empty
     * @param normalAllelicCounts       if not {@code null}, sites that are homozygous in the normal will be filtered out;
     *                                  if {@code null} (i.e., case-only),
     * @param copyRatioIntervals        never {@code null} (but may be empty), sites not overlapping with copy-ratio intervals will be filtered out
     * @return {@link NaiveHeterozygousPileupGenotypingResult}, with hetNormalAllelicCounts set to {@code null} for case-only mode
     */
    public static NaiveHeterozygousPileupGenotypingResult genotypeHets(final List<AllelicCountCollection> allelicCountsPerSample,
                                                                       final AllelicCountCollection normalAllelicCounts,
                                                                       final SimpleIntervalCollection copyRatioIntervals,
                                                                       final int minTotalAlleleCountCase,
                                                                       final int minTotalAlleleCountNormal,
                                                                       final double genotypingHomozygousLogRatioThreshold,
                                                                       final double genotypingBaseErrorRate) {
        Utils.nonEmpty(allelicCountsPerSample);
        Utils.nonNull(copyRatioIntervals);
        ParamUtils.isPositiveOrZero(minTotalAlleleCountCase, "Minimum total allele count for the case sample must be positive or zero.");
        ParamUtils.isPositiveOrZero(minTotalAlleleCountNormal, "Minimum total allele count for the normal sample must be positive or zero.");
        ParamUtils.isFinite(genotypingHomozygousLogRatioThreshold, "Genotyping homozygous log-ratio threshold must be finite.");
        ParamUtils.inRange(genotypingBaseErrorRate, 0., 1., "Genotyping base-error rate must be in [0, 1].");
        Utils.validateArg((int) Stream.of(allelicCountsPerSample, Collections.singletonList(normalAllelicCounts))
                        .flatMap(Collection::stream)
                        .filter(Objects::nonNull)
                        .map(AllelicCountCollection::getIntervals)
                        .distinct()
                        .count() == 1,
                "Allelic-count sites must be identical across all samples.");
        CopyNumberArgumentValidationUtils.getValidatedSequenceDictionary(
                Stream.of(allelicCountsPerSample, Arrays.asList(normalAllelicCounts, copyRatioIntervals))
                        .flatMap(Collection::stream)
                        .toArray(AbstractLocatableCollection[]::new));

        logger.info("Genotyping heterozygous sites from available allelic counts...");

        final List<AllelicCountCollection> hetAllelicCountsPerSample = new ArrayList<>(allelicCountsPerSample.size());
        final AllelicCountCollection hetNormalAllelicCounts;

        if (normalAllelicCounts != null) {
            logger.info("Matched normal was provided, running in matched-normal mode...");

            final SampleLocatableMetadata normalMetadata = normalAllelicCounts.getMetadata();
            final String normalSampleName = normalMetadata.getSampleName();

            //filter on total count in matched normal
            AllelicCountCollection filteredNormalAllelicCounts = filterByTotalCount(normalAllelicCounts, minTotalAlleleCountNormal);
            logger.info(String.format("Retained %d / %d sites after filtering allelic counts with total count less than %d in matched-normal sample %s...",
                    filteredNormalAllelicCounts.size(), normalAllelicCounts.size(), minTotalAlleleCountNormal, normalSampleName));

            //filter matched normal on overlap with copy-ratio intervals, if available
            if (copyRatioIntervals.size() > 0) {
                filteredNormalAllelicCounts = filterByOverlap(filteredNormalAllelicCounts, copyRatioIntervals);
                logger.info(String.format("Retained %d / %d sites after filtering on overlap with copy-ratio intervals in matched-normal sample %s...",
                        filteredNormalAllelicCounts.size(), normalAllelicCounts.size(), normalSampleName));
            }

            //filter on heterozygosity in matched normal
            hetNormalAllelicCounts = filterByHeterozygosity(filteredNormalAllelicCounts, genotypingHomozygousLogRatioThreshold, genotypingBaseErrorRate);
            logger.info(String.format("Retained %d / %d sites after filtering on heterozygosity in matched-normal sample %s...",
                    hetNormalAllelicCounts.size(), normalAllelicCounts.size(), normalSampleName));
        } else {
            logger.info("No matched normal was provided, not running in matched-normal mode...");
            hetNormalAllelicCounts = null;
        }

        //filter and genotype all case samples
        for (final AllelicCountCollection allelicCounts : allelicCountsPerSample) {
            final String sampleName = allelicCounts.getMetadata().getSampleName();

            //filter on total count in case sample
            AllelicCountCollection filteredAllelicCounts = filterByTotalCount(allelicCounts, minTotalAlleleCountCase);
            logger.info(String.format("Retained %d / %d sites after filtering allelic counts with total count less than %d in case sample %s...",
                    filteredAllelicCounts.size(), allelicCounts.size(), minTotalAlleleCountCase, sampleName));

            //filter on overlap with copy-ratio intervals, if available
            if (copyRatioIntervals.size() > 0) {
                filteredAllelicCounts = filterByOverlap(filteredAllelicCounts, copyRatioIntervals);
                logger.info(String.format("Retained %d / %d sites after filtering on overlap with copy-ratio intervals in case sample %s...",
                        filteredAllelicCounts.size(), allelicCounts.size(), sampleName));
            }

            final AllelicCountCollection hetAllelicCounts;
            if (hetNormalAllelicCounts != null) {
                //retrieve sites that were heterozygous in normal from the case sample
                logger.info(String.format("Retaining allelic counts for case sample %s at heterozygous sites in matched-normal sample %s...",
                        hetNormalAllelicCounts.getMetadata().getSampleName(), sampleName));
                hetAllelicCounts = filterByOverlap(filteredAllelicCounts, hetNormalAllelicCounts);
            } else {
                hetAllelicCounts = filterByHeterozygosity(filteredAllelicCounts, genotypingHomozygousLogRatioThreshold, genotypingBaseErrorRate);
                logger.info(String.format("Retained %d / %d sites after filtering on heterozygosity in case sample %s...",
                        hetAllelicCounts.size(), allelicCounts.size(), sampleName));
            }
            logger.info(String.format("Retained %d / %d sites after applying all filters to case sample %s.",
                    hetAllelicCounts.size(), allelicCounts.size(), sampleName));
            hetAllelicCountsPerSample.add(hetAllelicCounts);
        }

        if (hetAllelicCountsPerSample.size() > 1) {
            //filter by intersection of heterozygous sites in all case samples
            final int maxNumHetSites = hetAllelicCountsPerSample.stream()
                    .mapToInt(AllelicCountCollection::size)
                    .max().getAsInt();
            final Multiset<SimpleInterval> hetSitesMultiset = HashMultiset.create(maxNumHetSites);
            hetAllelicCountsPerSample.stream()
                    .map(AllelicCountCollection::getIntervals)
                    .forEach(hetSitesMultiset::addAll);
            final List<SimpleInterval> hetSitesIntersectionList = hetSitesMultiset.entrySet().stream()
                    .filter(e -> e.getCount() == hetAllelicCountsPerSample.size())
                    .map(Multiset.Entry::getElement)
                    .collect(Collectors.toList());
            final SimpleIntervalCollection hetSitesIntersection = new SimpleIntervalCollection(
                    hetAllelicCountsPerSample.get(0).getMetadata(),
                    hetSitesIntersectionList);
            IntStream.range(0, hetAllelicCountsPerSample.size())
                    .forEach(i -> hetAllelicCountsPerSample.set(i, filterByOverlap(hetAllelicCountsPerSample.get(i), hetSitesIntersection)));
            logger.info(String.format("Retained %d / %d sites after taking intersection of sites in all case samples.",
                    hetSitesIntersection.size(), allelicCountsPerSample.get(0).size()));
        }

        return new NaiveHeterozygousPileupGenotypingResult(hetAllelicCountsPerSample, hetNormalAllelicCounts);
    }

    private static AllelicCountCollection filterByTotalCount(final AllelicCountCollection allelicCounts,
                                                             final int minTotalAlleleCount) {
        return minTotalAlleleCount == 0
                ? allelicCounts
                : new AllelicCountCollection(
                        allelicCounts.getMetadata(),
                        allelicCounts.getRecords().stream()
                                .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCount)
                                .collect(Collectors.toList()));
    }

    private static <T extends Locatable> AllelicCountCollection filterByOverlap(final AllelicCountCollection allelicCounts,
                                                                                final AbstractLocatableCollection<?, T> locatableCollection) {
        return locatableCollection.size() == 0
                ? new AllelicCountCollection(
                        allelicCounts.getMetadata(),
                        Collections.emptyList())
                : new AllelicCountCollection(
                        allelicCounts.getMetadata(),
                        allelicCounts.getRecords().stream()
                            .filter(ac -> locatableCollection.getOverlapDetector().overlapsAny(ac))
                                .collect(Collectors.toList()));
    }

    private static AllelicCountCollection filterByHeterozygosity(final AllelicCountCollection allelicCounts,
                                                                 final double genotypingHomozygousLogRatioThreshold,
                                                                 final double genotypingBaseErrorRate) {
        return new AllelicCountCollection(
                allelicCounts.getMetadata(),
                allelicCounts.getRecords().stream()
                        .filter(ac -> calculateHomozygousLogRatio(ac, genotypingBaseErrorRate) < genotypingHomozygousLogRatioThreshold)
                        .collect(Collectors.toList()));
    }

    /**
     * Calculates a log likelihood ratio to perform a relatively naive test for heterozygosity.
     * For a binomial likelihood, returns the log of the ratio of
     *  1) the likelihood of the success probability being in [0, {@code genotypingBaseErrorRate}] or [1 - {@code genotypingBaseErrorRate}, 1], and
     *  2) the likelihood of the success probability being in [{@code genotypingBaseErrorRate}, 1 - {@code genotypingBaseErrorRate}],
     * given the count data represented by {@code allelicCount}.
     */
    private static double calculateHomozygousLogRatio(final AllelicCount allelicCount,
                                                      final double genotypingBaseErrorRate) {
        final int r = allelicCount.getRefReadCount();
        final int n = allelicCount.getTotalReadCount();
        final double betaAll = Beta.regularizedBeta(1, r + 1, n - r + 1);
        final double betaError = Beta.regularizedBeta(genotypingBaseErrorRate, r + 1, n - r + 1);
        final double betaOneMinusError = Beta.regularizedBeta(1 - genotypingBaseErrorRate, r + 1, n - r + 1);
        final double betaHom = betaError + betaAll - betaOneMinusError;
        final double betaHet = betaOneMinusError - betaError;
        return FastMath.log(betaHom) - FastMath.log(betaHet);
    }
}
