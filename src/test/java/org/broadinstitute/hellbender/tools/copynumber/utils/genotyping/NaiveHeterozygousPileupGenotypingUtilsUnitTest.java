package org.broadinstitute.hellbender.tools.copynumber.utils.genotyping;

import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public final class NaiveHeterozygousPileupGenotypingUtilsUnitTest extends GATKBaseTest {
    private static final int RANDOM_SEED = 1;   //reset seed before each simulated test case

    private static final int NUM_SITES = 40;
    private static final int NUM_SITES_PER_BLOCK = 10;

    private static final SAMSequenceDictionary SEQUENCE_DICTIONARY = new SAMSequenceDictionary(
            Arrays.asList(
                    new SAMSequenceRecord("chr1", NUM_SITES),
                    new SAMSequenceRecord("chr2", NUM_SITES)));

    private static final SimpleIntervalCollection COPY_RATIO_INTERVALS = new SimpleIntervalCollection(
            new SimpleLocatableMetadata(SEQUENCE_DICTIONARY),
            Collections.singletonList(new SimpleInterval("chr1", 1, NUM_SITES)));
    private static final List<Boolean> COPY_RATIO_INTERVALS_IS_BLOCK_PASSING_OVERLAP_FILTER =
            Arrays.asList(true, true, true, true);
    private static final SimpleIntervalCollection EMPTY_COPY_RATIO_INTERVALS = new SimpleIntervalCollection(
            new SimpleLocatableMetadata(SEQUENCE_DICTIONARY),
            Collections.emptyList());
    private static final List<Boolean> EMPTY_COPY_RATIO_INTERVALS_IS_BLOCK_PASSING_OVERLAP_FILTER =
            Arrays.asList(true, true, true, true);
    private static final SimpleIntervalCollection NON_OVERLAPPING_COPY_RATIO_INTERVALS = new SimpleIntervalCollection(
            new SimpleLocatableMetadata(SEQUENCE_DICTIONARY),
            Collections.singletonList(new SimpleInterval("chr2", 1, NUM_SITES)));
    private static final List<Boolean> NON_OVERLAPPING_COPY_RATIO_INTERVALS_IS_BLOCK_PASSING_OVERLAP_FILTER =
            Arrays.asList(false, false, false, false);

    private static final int MIN_TOTAL_ALLELE_COUNT_CASE_IN_MATCHED_NORMAL_MODE = 0;
    private static final int MIN_TOTAL_ALLELE_COUNT_CASE_IN_CASE_ONLY_MODE = 30;
    private static final int MIN_TOTAL_ALLELE_COUNT_NORMAL = 30;
    private static final double GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD = -10.;
    private static final double GENOTYPING_BASE_ERROR_RATE = 5E-2;

    /**
     * Note that this generative model does not actually map exactly to the likelihood specified in
     * {@link NaiveHeterozygousPileupGenotypingUtils#calculateHomozygousLogRatio}.
     * In particular, {@code epsilon} does not map to {@code genotypingBaseErrorRate}.
     * We are just trying to roughly generate homozygous sites with some base errors, but we are not being
     * too careful about the generative process.
     */
    private static List<Integer> generateAltCounts(final RandomGenerator rng,
                                                   final int numSites,
                                                   final int totalDepth,
                                                   final double alternateAlleleFraction,
                                                   final double epsilon) {
        final double p = 0 < alternateAlleleFraction && alternateAlleleFraction < 1
                ? alternateAlleleFraction
                : alternateAlleleFraction + epsilon * (alternateAlleleFraction == 0 ? 1 : -1);
        final BinomialDistribution binomialDistribution = new BinomialDistribution(rng, totalDepth, p);
        return Ints.asList(binomialDistribution.sample(numSites));
    }

    /**
     * Generate a block of sites, with the first half being hom-ref and the second half being hom-alt.
     * (We assume {@link NaiveHeterozygousPileupGenotypingUtilsUnitTest#NUM_SITES_PER_BLOCK} is even.
     */
    private static List<Integer> generateHomBlockAltCounts(final RandomGenerator rng,
                                                           final int totalDepth,
                                                           final double epsilon) {
        return Stream.of(
                generateAltCounts(rng, NUM_SITES_PER_BLOCK / 2, totalDepth, 0., epsilon),
                generateAltCounts(rng, NUM_SITES_PER_BLOCK / 2, totalDepth, 1., epsilon))
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
    }

    private static List<Integer> generateHetBlockAltCounts(final RandomGenerator rng,
                                                           final int totalDepth,
                                                           final double epsilon) {
        return generateAltCounts(rng, NUM_SITES_PER_BLOCK, totalDepth, 0.5, epsilon);
    }

    private static void addBlockAllelicCounts(final List<AllelicCount> allelicCounts,
                                              final List<Integer> blockAltCounts,
                                              final int totalDepth,
                                              final String contig) {
        final int start = allelicCounts.size() + 1;
        IntStream.range(0, blockAltCounts.size()).forEach(i ->
                allelicCounts.add(new AllelicCount(
                        new SimpleInterval(contig, start + i, start + i),
                        totalDepth - blockAltCounts.get(i),
                        blockAltCounts.get(i))));
    }

    /**
     * Resolves all sample-and/or-block-level filters to a single block-level filter,
     * using logic as appropriate for matched-normal and case-only modes.
     */
    private static List<Boolean> isBlockPassing(final List<List<Boolean>> casesIsBlockPassingDepthAndHetFilters,
                                                final List<Boolean> normalIsBlockPassingDepthAndHetFilters,
                                                final List<Boolean> isBlockPassingCopyRatioOverlapFilter) {
        return IntStream.range(0, casesIsBlockPassingDepthAndHetFilters.get(0).size()).boxed()
                .map(siteIndex ->
                        (normalIsBlockPassingDepthAndHetFilters != null
                            ? normalIsBlockPassingDepthAndHetFilters.get(siteIndex)             //matched-normal
                            : IntStream.range(0, casesIsBlockPassingDepthAndHetFilters.size())  //case-only
                                .allMatch(sampleIndex -> casesIsBlockPassingDepthAndHetFilters.get(sampleIndex).get(siteIndex))) &&
                        isBlockPassingCopyRatioOverlapFilter.get(siteIndex))
                .collect(Collectors.toList());
    }

    @DataProvider(name = "dataGenotypeHets")
    @SuppressWarnings("unchecked")
    public Object[][] dataGenotypeHets() {
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        final double epsilon = 0.02;

        //normal is 3 het blocks and 1 hom block (het|het|het|hom), all counts passing normal total-depth filter
        final int normalTotalDepth = 500;
        final List<AllelicCount> normalAllelicCountsList = new ArrayList<>(NUM_SITES);
        addBlockAllelicCounts(normalAllelicCountsList, generateHetBlockAltCounts(rng, normalTotalDepth, epsilon), normalTotalDepth, "chr1");
        addBlockAllelicCounts(normalAllelicCountsList, generateHetBlockAltCounts(rng, normalTotalDepth, epsilon), normalTotalDepth, "chr1");
        addBlockAllelicCounts(normalAllelicCountsList, generateHetBlockAltCounts(rng, normalTotalDepth, epsilon), normalTotalDepth, "chr1");
        addBlockAllelicCounts(normalAllelicCountsList, generateHomBlockAltCounts(rng, normalTotalDepth, epsilon), normalTotalDepth, "chr1");
        final AllelicCountCollection normalAllelicCounts = new AllelicCountCollection(
                new SimpleSampleLocatableMetadata("normal", SEQUENCE_DICTIONARY),
                normalAllelicCountsList);
        final List<Boolean> normalIsBlockPassingDepthAndHetFilters = Arrays.asList(true, true, true, false);

        //first case is (het|het|hom|hom), all counts passing case total-depth filter
        final int caseTotalDepth = 500;
        final List<AllelicCount> firstCaseAllelicCountsList = new ArrayList<>(NUM_SITES);
        addBlockAllelicCounts(firstCaseAllelicCountsList, generateHetBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        addBlockAllelicCounts(firstCaseAllelicCountsList, generateHetBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        addBlockAllelicCounts(firstCaseAllelicCountsList, generateHomBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        addBlockAllelicCounts(firstCaseAllelicCountsList, generateHomBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        final AllelicCountCollection firstCaseAllelicCounts = new AllelicCountCollection(
                new SimpleSampleLocatableMetadata("case-1", SEQUENCE_DICTIONARY),
                firstCaseAllelicCountsList);
        final List<Boolean> firstCaseIsBlockPassingDepthAndHetFilters = Arrays.asList(true, true, false, false);

        //low coverage version of the first case (het|het|hom|hom), all counts failing case total-depth filter
        final int caseLowCoverageTotalDepth = 5;
        final List<AllelicCount> firstCaseLowCoverageAllelicCountsList = new ArrayList<>(NUM_SITES);
        addBlockAllelicCounts(firstCaseLowCoverageAllelicCountsList, generateHetBlockAltCounts(rng, caseLowCoverageTotalDepth, epsilon), caseLowCoverageTotalDepth, "chr1");
        addBlockAllelicCounts(firstCaseLowCoverageAllelicCountsList, generateHetBlockAltCounts(rng, caseLowCoverageTotalDepth, epsilon), caseLowCoverageTotalDepth, "chr1");
        addBlockAllelicCounts(firstCaseLowCoverageAllelicCountsList, generateHomBlockAltCounts(rng, caseLowCoverageTotalDepth, epsilon), caseLowCoverageTotalDepth, "chr1");
        addBlockAllelicCounts(firstCaseLowCoverageAllelicCountsList, generateHomBlockAltCounts(rng, caseLowCoverageTotalDepth, epsilon), caseLowCoverageTotalDepth, "chr1");
        final AllelicCountCollection firstCaseLowCoverageAllelicCounts = new AllelicCountCollection(
                new SimpleSampleLocatableMetadata("case-1-low", SEQUENCE_DICTIONARY),
                firstCaseLowCoverageAllelicCountsList);
        final List<Boolean> firstCaseLowCoverageIsBlockPassingDepthAndHetFilters = Arrays.asList(false, false, false, false);

        //second case is (het|hom|het|hom), all counts passing case total-depth filter
        final List<AllelicCount> secondCaseAllelicCountsList = new ArrayList<>(NUM_SITES);
        addBlockAllelicCounts(secondCaseAllelicCountsList, generateHetBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        addBlockAllelicCounts(secondCaseAllelicCountsList, generateHomBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        addBlockAllelicCounts(secondCaseAllelicCountsList, generateHetBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        addBlockAllelicCounts(secondCaseAllelicCountsList, generateHomBlockAltCounts(rng, caseTotalDepth, epsilon), caseTotalDepth, "chr1");
        final AllelicCountCollection secondCaseAllelicCounts = new AllelicCountCollection(
                new SimpleSampleLocatableMetadata("case-2", SEQUENCE_DICTIONARY),
                secondCaseAllelicCountsList);
        final List<Boolean> secondCaseIsBlockPassingDepthAndHetFilters = Arrays.asList(true, false, true, false);

        final List<List<Object>> data = new ArrayList<>();
        for (final Pair<?, ?> normalDataFilterPair : Arrays.asList(
                Pair.of(normalAllelicCounts, normalIsBlockPassingDepthAndHetFilters),
                Pair.of(null, null))) {
            for (final Pair<SimpleIntervalCollection, List<Boolean>> copyRatioDataFilterPair : Arrays.asList(
                    Pair.of(COPY_RATIO_INTERVALS, COPY_RATIO_INTERVALS_IS_BLOCK_PASSING_OVERLAP_FILTER),
                    Pair.of(EMPTY_COPY_RATIO_INTERVALS, EMPTY_COPY_RATIO_INTERVALS_IS_BLOCK_PASSING_OVERLAP_FILTER),
                    Pair.of(NON_OVERLAPPING_COPY_RATIO_INTERVALS, NON_OVERLAPPING_COPY_RATIO_INTERVALS_IS_BLOCK_PASSING_OVERLAP_FILTER))) {
                final AllelicCountCollection normalData = (AllelicCountCollection) normalDataFilterPair.getLeft();
                final List<Boolean> normalFilter = (List<Boolean>) normalDataFilterPair.getRight();
                final SimpleIntervalCollection copyRatioData = copyRatioDataFilterPair.getLeft();
                final List<Boolean> copyRatioFilter = copyRatioDataFilterPair.getRight();

                //first case + second case
                data.add(Arrays.asList(
                        Arrays.asList(
                                firstCaseAllelicCounts,
                                secondCaseAllelicCounts),
                        normalData,
                        copyRatioData,
                        isBlockPassing(
                                Arrays.asList(
                                        firstCaseIsBlockPassingDepthAndHetFilters,
                                        secondCaseIsBlockPassingDepthAndHetFilters),
                                normalFilter,
                                copyRatioFilter)));
                //first case
                data.add(Arrays.asList(
                        Collections.singletonList(
                                firstCaseAllelicCounts),
                        normalData,
                        copyRatioData,
                        isBlockPassing(
                                Collections.singletonList(
                                        firstCaseIsBlockPassingDepthAndHetFilters),
                                normalFilter,
                                copyRatioFilter)));
                //low coverage first case
                data.add(Arrays.asList(
                        Collections.singletonList(
                                firstCaseLowCoverageAllelicCounts),
                        normalData,
                        copyRatioData,
                        isBlockPassing(
                                Collections.singletonList(
                                        firstCaseLowCoverageIsBlockPassingDepthAndHetFilters),
                                normalFilter,
                                copyRatioFilter)));
            }
        }

        return data.stream().map(x -> x.toArray(new Object[0])).toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataGenotypeHets")
    public void testGenotypeHets(final List<AllelicCountCollection> allelicCountsPerSample,
                                 final AllelicCountCollection normalAllelicCounts,
                                 final SimpleIntervalCollection copyRatioIntervals,
                                 final List<Boolean> isBlockPassingExpected) {
        final NaiveHeterozygousPileupGenotypingUtils.NaiveHeterozygousPileupGenotypingResult result =
                NaiveHeterozygousPileupGenotypingUtils.genotypeHets(
                        allelicCountsPerSample, normalAllelicCounts, copyRatioIntervals,
                        normalAllelicCounts != null
                                ? MIN_TOTAL_ALLELE_COUNT_CASE_IN_MATCHED_NORMAL_MODE
                                : MIN_TOTAL_ALLELE_COUNT_CASE_IN_CASE_ONLY_MODE,
                        MIN_TOTAL_ALLELE_COUNT_NORMAL,
                        GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD,
                        GENOTYPING_BASE_ERROR_RATE);

        Assert.assertEquals(normalAllelicCounts == null, result.getHetNormalAllelicCounts() == null);

        //check that all sites across all samples in result are the same
        final List<List<SimpleInterval>> casesPassingSites = result.getHetAllelicCountsPerSample().stream()
                .map(AllelicCountCollection::getIntervals)
                .collect(Collectors.toList());
        Assert.assertEquals(casesPassingSites.stream().distinct().count(), 1);
        if (result.getHetNormalAllelicCounts() != null) {
            final List<SimpleInterval> normalPassingSites = result.getHetNormalAllelicCounts().getIntervals();
            Assert.assertEquals(casesPassingSites.get(0), normalPassingSites);
        }

        //check expected passing counts across all samples in result
        final List<Boolean> isCountPassingExpected = isBlockPassingExpected.stream()
                .map(b -> Collections.nCopies(NUM_SITES_PER_BLOCK, b))
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
        for (int caseIndex = 0; caseIndex < allelicCountsPerSample.size(); caseIndex++) {
            final Set<AllelicCount> casePassingAllelicCounts = new HashSet<>(result.getHetAllelicCountsPerSample().get(caseIndex).getRecords());
            final List<Boolean> isCountPassing = allelicCountsPerSample.get(caseIndex).getRecords().stream()
                    .map(casePassingAllelicCounts::contains)
                    .collect(Collectors.toList());
            Assert.assertEquals(isCountPassing, isCountPassingExpected);
        }
        if (normalAllelicCounts != null) {
            final Set<AllelicCount> normalPassingAllelicCounts = new HashSet<>(result.getHetNormalAllelicCounts().getRecords());
            final List<Boolean> isCountPassing = normalAllelicCounts.getRecords().stream()
                    .map(normalPassingAllelicCounts::contains)
                    .collect(Collectors.toList());
            Assert.assertEquals(isCountPassing, isCountPassingExpected);
        }
    }
}
