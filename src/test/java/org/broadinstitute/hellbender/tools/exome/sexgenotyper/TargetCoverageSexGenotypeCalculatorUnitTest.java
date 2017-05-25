package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import avro.shaded.com.google.common.collect.ImmutableMap;
import avro.shaded.com.google.common.collect.Sets;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link TargetCoverageSexGenotypeCalculator}.
 *
 * The test includes a small genotyping task on real de-identified data, as well as multiple tests on
 * random data generated according to the genotyping model (Poisson read count distribution w/ bait bias).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class TargetCoverageSexGenotypeCalculatorUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File TEST_CONTIG_PLOIDY_ANNOTS_FILE = new File(TEST_SUB_DIR, "contig_annots.tsv");
    private static final File TEST_RCC_FILE = new File(TEST_SUB_DIR, "sex_genotyper_rcc_trunc.tsv");
    private static final File TEST_SEX_GENOTYPE_FILE = new File(TEST_SUB_DIR, "sex_genotypes_agilent_trunc.tsv");
    private static final double MAPPING_ERROR_PROBABILITY = 1e-3;

    private static final Set<String> SEX_GENOTYPES = new HashSet<>(Arrays.asList("SEX_XX", "SEX_XY"));
    private static final Set<String> AUTOSOMAL_CONTIGS = new HashSet<>(Arrays.asList("1", "2", "3", "4", "5"));
    private static final Set<String> ALLOSOMAL_CONTIGS = new HashSet<>(Arrays.asList("X", "Y"));
    private static final int AUTOSOMAL_CONTIG_PLOIDY = 2;
    private static final Map<String, Map<String, Integer>> ALLOSOMAL_CONTIG_PLOIDY_MAP = ImmutableMap.of(
            "SEX_XX", ImmutableMap.of("X", 2, "Y" , 0),
            "SEX_XY", ImmutableMap.of("X", 1, "Y" , 1));

    private static final int NUMBER_OF_RANDOM_SAMPLES = 20;
    private static double MIN_BAIT_COUNT = TargetCoverageSexGenotypeCalculator.MIN_ALLOWED_BAIT_COUNT_PER_TARGET + 1.0;
    private static double MAX_BAIT_COUNT = 100.0;
    private static int MIN_TARGETS_PER_CONTIG = 10;
    private static int MAX_TARGETS_PER_CONTIG = 20;
    private static int MAX_TARGET_LENGTH = 1000;
    private static int MAX_TARGET_SPACING = 1000;
    private static double MIN_READ_DEPTH = 200.0;
    private static double MAX_READ_DEPTH = 500.0;
    private static long RANDOM_SEED = 1984;
    private static final double EPSILON = 1e-12;

    private static final List<ContigGermlinePloidyAnnotation> HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS;
    static {
        HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS = ContigGermlinePloidyAnnotationTableReader
                .readContigGermlinePloidyAnnotationsFromFile(TEST_CONTIG_PLOIDY_ANNOTS_FILE);
    }

    private static final List<Target> RANDOM_TARGETS_WITH_BAIT_COUNTS = annotateTargetsWithRandomBaitCounts(
            generateRandomTargets(new Random(RANDOM_SEED)), new Random(RANDOM_SEED));
    private static final List<Target> RANDOM_TARGETS_WITHOUT_BAIT_COUNTS = generateRandomTargets(new Random(RANDOM_SEED));

    private static final ReadCountCollection RANDOM_READ_COUNTS_WITH_BAIT_BIAS;
    private static final List<String> TRUTH_SEX_GENOTYPES_WITH_BAIT_BIAS;

    private static final ReadCountCollection RANDOM_READ_COUNTS_WITHOUT_BAIT_BIAS;
    private static final List<String> TRUTH_SEX_GENOTYPES_WITHOUT_BAIT_BIAS;

    static {
        final ImmutablePair<List<String>, ReadCountCollection> randomData = generateRandomReadCountCollection(
                RANDOM_TARGETS_WITH_BAIT_COUNTS, NUMBER_OF_RANDOM_SAMPLES, new Well19937c(RANDOM_SEED), true);
        RANDOM_READ_COUNTS_WITH_BAIT_BIAS = randomData.getRight();
        TRUTH_SEX_GENOTYPES_WITH_BAIT_BIAS = randomData.getLeft();
    }

    static {
        final ImmutablePair<List<String>, ReadCountCollection> randomData = generateRandomReadCountCollection(
                RANDOM_TARGETS_WITHOUT_BAIT_COUNTS, NUMBER_OF_RANDOM_SAMPLES, new Well19937c(RANDOM_SEED), false);
        RANDOM_READ_COUNTS_WITHOUT_BAIT_BIAS = randomData.getRight();
        TRUTH_SEX_GENOTYPES_WITHOUT_BAIT_BIAS = randomData.getLeft();
    }

    /**
     * Note: this test relies on the content of the test resource files
     */
    @Test
    public static void testOnRealData() {
        final ReadCountCollection readCounts;
        try {
            readCounts = ReadCountCollectionUtils.parse(TEST_RCC_FILE);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource files");
        }
        final TargetCoverageSexGenotypeCalculator genotypeCalculator = new TargetCoverageSexGenotypeCalculator(readCounts,
                readCounts.targets(), HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);

        /* assert all ploidy tags are present */
        final Set<String> loadedSexGenotypes = genotypeCalculator.getSexGenotypeIdentifiers();
        Assert.assertTrue(loadedSexGenotypes.equals(SEX_GENOTYPES));

        /* assert autosomal and allosomal target lists */
        Assert.assertEquals(genotypeCalculator.getAutosomalTargetList(), readCounts.targets().stream()
            .filter(target -> AUTOSOMAL_CONTIGS.contains(target.getContig()))
            .collect(Collectors.toList()));
        Assert.assertEquals(genotypeCalculator.getAllosomalTargetList(), readCounts.targets().stream()
                .filter(target -> ALLOSOMAL_CONTIGS.contains(target.getContig()))
                .collect(Collectors.toList()));

        /* assert ploidy on autosomal targets */
        final int[] autosomalTargetPloidies = genotypeCalculator.getAutosomalTargetGermlinePloidies();
        ArrayAsserts.assertArrayEquals(
                IntStream.range(0, autosomalTargetPloidies.length).map(i -> AUTOSOMAL_CONTIG_PLOIDY).toArray(),
                autosomalTargetPloidies);

        /* assert ploidy on allosomal targets */
        final Map<String, int[]> allosomalTargetPloidies = genotypeCalculator.getAllosomalTargetGermlinePloidiesMap();
        int[] expected;
        final List<Target> allosomalTargetList = genotypeCalculator.getAllosomalTargetList();
        expected = allosomalTargetList.stream()
                .mapToInt(t -> ALLOSOMAL_CONTIG_PLOIDY_MAP.get("SEX_XX").get(t.getContig()))
                .toArray();
        ArrayAsserts.assertArrayEquals(expected, allosomalTargetPloidies.get("SEX_XX"));
        expected = allosomalTargetList.stream()
                .mapToInt(t -> ALLOSOMAL_CONTIG_PLOIDY_MAP.get("SEX_XY").get(t.getContig()))
                .toArray();
        ArrayAsserts.assertArrayEquals(expected, allosomalTargetPloidies.get("SEX_XY"));

        /* assert the correctness of inferred sex genotypes */
        final SexGenotypeDataCollection result = genotypeCalculator.inferSexGenotypes();
        try {
            Assert.assertEquals(result, new SexGenotypeDataCollection(TEST_SEX_GENOTYPE_FILE));
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource files");
        }
    }

    /**
     * Generates a random list of autosomal and allosomal targets
     */
    private static List<Target> generateRandomTargets(final Random rng) {
        final List<Target> targetList = new ArrayList<>();
        for (final String contig : Sets.union(AUTOSOMAL_CONTIGS, ALLOSOMAL_CONTIGS)) {
            final int numTargets = MIN_TARGETS_PER_CONTIG + rng.nextInt(MAX_TARGETS_PER_CONTIG - MIN_TARGETS_PER_CONTIG);
            int lastEndpoint = 1;
            for (int ti = 0; ti < numTargets; ti++) {
                final int start = lastEndpoint + rng.nextInt(MAX_TARGET_SPACING) + 1;
                final int end = start + rng.nextInt(MAX_TARGET_LENGTH) + 1;
                final String name = String.format("%s_%d_%d", contig, start, end);
                targetList.add(new Target(name, new SimpleInterval(contig, start, end)));
                lastEndpoint = end;
            }
        }
        return targetList;
    }

    private static List<Target> annotateTargetsWithRandomBaitCounts(final List<Target> targets, final Random rng) {
        return targets.stream()
                .map(target -> new Target(target.getName(), target.getInterval(),
                        new TargetAnnotationCollection(ImmutableMap.of(TargetAnnotation.BAIT_COUNT,
                                String.valueOf(MIN_BAIT_COUNT + rng.nextDouble() * (MAX_BAIT_COUNT - MIN_BAIT_COUNT))))))
                .collect(Collectors.toList());
    }

    /**
     * Generates a random read count collection with a given number of samples with random sex genotypes.
     * Bait count bias is taken into account.
     *
     * @param targets a valid target list
     * @param numSamples number of samples to generate
     * @return pair of sex genotypes list and read counts
     */
    private static ImmutablePair<List<String>, ReadCountCollection> generateRandomReadCountCollection(
            final List<Target> targets, final int numSamples, final RandomGenerator rng, final boolean withBaitCountBias) {
        /* generate sample names and sex genotypes */
        final List<String> sampleNames = IntStream.range(0, numSamples)
                .mapToObj(si -> String.format("SAMPLE_%d", si))
                .collect(Collectors.toList());
        final int numTargets = targets.size();
        final String[] sexGenotypesArray = SEX_GENOTYPES.toArray(new String[0]);
        final List<String> sampleSexGenotypes = IntStream.range(0, numSamples)
                .mapToObj(si -> sexGenotypesArray[rng.nextInt(sexGenotypesArray.length)])
                .collect(Collectors.toList());

        /* draw random read counts from a Poisson distribution */
        final RealMatrix readCounts = new BlockRealMatrix(numTargets, numSamples);
        final double[] baitCounts;
        if (withBaitCountBias) {
            baitCounts = targets.stream()
                    .map(Target::getAnnotations)
                    .mapToDouble(annotations -> annotations.getDouble(TargetAnnotation.BAIT_COUNT))
                    .toArray();
        } else {
            baitCounts = new double[0];
        }
        for (int si = 0; si < numSamples; si++) {
            final String sampleSexGenotype = sampleSexGenotypes.get(si);
            final double readDepth = MIN_READ_DEPTH + rng.nextDouble() * (MAX_READ_DEPTH - MIN_READ_DEPTH);
            for (int ti = 0; ti < numTargets; ti++) {
                final String contig = targets.get(ti).getContig();
                final int ploidy = AUTOSOMAL_CONTIGS.contains(contig)
                        ? AUTOSOMAL_CONTIG_PLOIDY
                        : ALLOSOMAL_CONTIG_PLOIDY_MAP.get(sampleSexGenotype).get(contig);
                final double lambda = (1 - MAPPING_ERROR_PROBABILITY) * readDepth *
                        (withBaitCountBias ? baitCounts[ti] : 1.0) * ploidy +
                        MAPPING_ERROR_PROBABILITY * readDepth * AUTOSOMAL_CONTIG_PLOIDY;
                final int readCount = new PoissonDistribution(rng, lambda,
                        PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS).sample();
                readCounts.setEntry(ti, si, readCount);
            }
        }
        return ImmutablePair.of(sampleSexGenotypes, new ReadCountCollection(targets, sampleNames, readCounts));
    }

    @Test
    public void testOnRandomPoissonReadCountsWithBaitCountBias() {
        final TargetCoverageSexGenotypeCalculator genotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                RANDOM_READ_COUNTS_WITH_BAIT_BIAS, RANDOM_TARGETS_WITH_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);
        final List<String> inferredSexGenotypes = genotypeCalculator.inferSexGenotypes().getSexGenotypeDataList()
                .stream().map(SexGenotypeData::getSexGenotype).collect(Collectors.toList());
        Assert.assertEquals(TRUTH_SEX_GENOTYPES_WITH_BAIT_BIAS, inferredSexGenotypes);
    }

    @Test
    public void testOnRandomPoissonReadCountsWithoutBaitCountBias() {
        final TargetCoverageSexGenotypeCalculator genotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                RANDOM_READ_COUNTS_WITHOUT_BAIT_BIAS, RANDOM_TARGETS_WITHOUT_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);
        final List<String> inferredSexGenotypes = genotypeCalculator.inferSexGenotypes().getSexGenotypeDataList()
                .stream().map(SexGenotypeData::getSexGenotype).collect(Collectors.toList());
        Assert.assertEquals(TRUTH_SEX_GENOTYPES_WITHOUT_BAIT_BIAS, inferredSexGenotypes);
    }

    @Test
    public void testOnRandomPoissonReadCountsWithExcludedIntervals() {
        final int maxContigLength = (MAX_TARGET_LENGTH + MAX_TARGET_SPACING) * MAX_TARGETS_PER_CONTIG * 2;
        final List<SimpleInterval> genotypingIntervals = new ArrayList<>();
        genotypingIntervals.add(new SimpleInterval("1", 1, maxContigLength));
        genotypingIntervals.add(new SimpleInterval("X", 1, maxContigLength));
        final TargetCoverageSexGenotypeCalculator genotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                RANDOM_READ_COUNTS_WITH_BAIT_BIAS, RANDOM_TARGETS_WITH_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, genotypingIntervals, MAPPING_ERROR_PROBABILITY);
        final List<String> inferredSexGenotypes = genotypeCalculator.inferSexGenotypes().getSexGenotypeDataList()
                .stream().map(SexGenotypeData::getSexGenotype).collect(Collectors.toList());
        Assert.assertEquals(TRUTH_SEX_GENOTYPES_WITH_BAIT_BIAS, inferredSexGenotypes);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoAllosomalTargetsException() {
        final int maxContigLength = (MAX_TARGET_LENGTH + MAX_TARGET_SPACING) * MAX_TARGETS_PER_CONTIG * 2;
        final List<SimpleInterval> genotypingIntervals = new ArrayList<>();
        genotypingIntervals.add(new SimpleInterval("1", 1, maxContigLength));
        new TargetCoverageSexGenotypeCalculator(
                RANDOM_READ_COUNTS_WITH_BAIT_BIAS, RANDOM_TARGETS_WITH_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, genotypingIntervals, MAPPING_ERROR_PROBABILITY);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoAutosomalTargetsException() {
        final int maxContigLength = (MAX_TARGET_LENGTH + MAX_TARGET_SPACING) * MAX_TARGETS_PER_CONTIG * 2;
        final List<SimpleInterval> genotypingIntervals = new ArrayList<>();
        genotypingIntervals.add(new SimpleInterval("X", 1, maxContigLength));
        genotypingIntervals.add(new SimpleInterval("Y", 1, maxContigLength));
        new TargetCoverageSexGenotypeCalculator(
                RANDOM_READ_COUNTS_WITH_BAIT_BIAS, RANDOM_TARGETS_WITH_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, genotypingIntervals, MAPPING_ERROR_PROBABILITY);
    }

    /**
     * If there is just a single sample, totally uncovered targets must not be masked
     */
    @Test
    public void testTotallyUncoveredTargetsMaskedSingleSample() {
        final int numTargets = RANDOM_TARGETS_WITH_BAIT_COUNTS.size();
        final RealMatrix readCountsMatrix = new BlockRealMatrix(numTargets, 1);
        final double[] data = new double[numTargets];
        Arrays.fill(data, 1.0);
        readCountsMatrix.setColumn(0, data);
        final ReadCountCollection singleSampleReadCountCollection = new ReadCountCollection(RANDOM_TARGETS_WITH_BAIT_COUNTS,
                Collections.singletonList("SOME_SAMPLE_NAME"), readCountsMatrix);
        final TargetCoverageSexGenotypeCalculator sexGenotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                singleSampleReadCountCollection, RANDOM_TARGETS_WITH_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);
        final int[] autosomalTargetMasks = sexGenotypeCalculator.getAutosomalTargetMasks();
        final int[] allosomalTargetMasks = sexGenotypeCalculator.getAllosomalTargetMasks();

        /* generate the expected mask */
        final int numAutosomalTargets = (int)RANDOM_TARGETS_WITH_BAIT_COUNTS.stream()
                .filter(target -> AUTOSOMAL_CONTIGS.contains(target.getContig()))
                .count();
        final int numAllosomalTargets = (int)RANDOM_TARGETS_WITH_BAIT_COUNTS.stream()
                .filter(target -> ALLOSOMAL_CONTIGS.contains(target.getContig()))
                .count();
        final int[] expectedAutosomalTargetMasks = new int[numAutosomalTargets];
        Arrays.fill(expectedAutosomalTargetMasks, 1);
        final int[] expectedAllosomalTargetMasks = new int[numAllosomalTargets];
        Arrays.fill(expectedAllosomalTargetMasks, 1);

        ArrayAsserts.assertArrayEquals(autosomalTargetMasks, expectedAutosomalTargetMasks);
        ArrayAsserts.assertArrayEquals(allosomalTargetMasks, expectedAllosomalTargetMasks);
    }

    /**
     * Must be able to detect and mask totally uncovered targets
     */
    @Test
    public void testTotallyUncoveredTargetsMaskedMultipleSample() {
        final int numTargets = RANDOM_TARGETS_WITH_BAIT_COUNTS.size();
        /* vanishing counts on every other target */
        final RealMatrix readCountsMatrix = new BlockRealMatrix(numTargets, NUMBER_OF_RANDOM_SAMPLES);
        final double[] data = new double[numTargets];
        Arrays.fill(data, 1.0);
        IntStream.range(0, numTargets/2).forEach(ti -> data[2*ti] = 0);
        for (int si = 0; si < NUMBER_OF_RANDOM_SAMPLES; si++) {
            readCountsMatrix.setColumn(si, data);
        }
        final Set<Target> badTargets = IntStream.range(0, numTargets/2)
                .mapToObj(ti -> RANDOM_TARGETS_WITH_BAIT_COUNTS.get(2*ti))
                .collect(Collectors.toSet());
        final ReadCountCollection readCountCollection = new ReadCountCollection(RANDOM_TARGETS_WITH_BAIT_COUNTS,
                RANDOM_READ_COUNTS_WITH_BAIT_BIAS.columnNames(), readCountsMatrix);

        final TargetCoverageSexGenotypeCalculator sexGenotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                readCountCollection, RANDOM_TARGETS_WITH_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);
        final int[] autosomalTargetMasks = sexGenotypeCalculator.getAutosomalTargetMasks();
        final int[] allosomalTargetMasks = sexGenotypeCalculator.getAllosomalTargetMasks();

        /* generate the expected mask */
        final int numAutosomalTargets = (int)RANDOM_TARGETS_WITH_BAIT_COUNTS.stream()
                .filter(target -> AUTOSOMAL_CONTIGS.contains(target.getContig()))
                .count();
        final int numAllosomalTargets = (int)RANDOM_TARGETS_WITH_BAIT_COUNTS.stream()
                .filter(target -> ALLOSOMAL_CONTIGS.contains(target.getContig()))
                .count();
        final int[] expectedAutosomalTargetMasks = new int[numAutosomalTargets];
        Arrays.fill(expectedAutosomalTargetMasks, 1);
        final int[] expectedAllosomalTargetMasks = new int[numAllosomalTargets];
        Arrays.fill(expectedAllosomalTargetMasks, 1);
        IntStream.range(0, numAutosomalTargets)
                .filter(ti -> badTargets.contains(sexGenotypeCalculator.getAutosomalTargetList().get(ti)))
                .forEach(ti -> expectedAutosomalTargetMasks[ti] = 0);
        IntStream.range(0, numAllosomalTargets)
                .filter(ti -> badTargets.contains(sexGenotypeCalculator.getAllosomalTargetList().get(ti)))
                .forEach(ti -> expectedAllosomalTargetMasks[ti] = 0);

        ArrayAsserts.assertArrayEquals(autosomalTargetMasks, expectedAutosomalTargetMasks);
        ArrayAsserts.assertArrayEquals(allosomalTargetMasks, expectedAllosomalTargetMasks);
    }

    /**
     * Targets with low bait counts must be masked
     */
    @Test
    public void testLowBaitCountMasked() {
        final int numTargets = RANDOM_TARGETS_WITHOUT_BAIT_COUNTS.size();
        final RealMatrix readCountsMatrix = new BlockRealMatrix(numTargets, NUMBER_OF_RANDOM_SAMPLES);
        final double[] data = new double[numTargets];
        Arrays.fill(data, 1.0);
        for (int si = 0; si < NUMBER_OF_RANDOM_SAMPLES; si++) {
            readCountsMatrix.setColumn(si, data);
        }
        final ReadCountCollection readCountCollection = new ReadCountCollection(RANDOM_TARGETS_WITHOUT_BAIT_COUNTS,
                RANDOM_READ_COUNTS_WITHOUT_BAIT_BIAS.columnNames(), readCountsMatrix);

        /* low bait counts on every other target */
        final List<Target> targetsWithBadBaits = IntStream.range(0, numTargets)
                .mapToObj(ti -> {
                    final Target target = RANDOM_TARGETS_WITHOUT_BAIT_COUNTS.get(ti);
                    final String baitCount = ti % 2 == 0
                            ? String.valueOf(0.5 * TargetCoverageSexGenotypeCalculator.MIN_ALLOWED_BAIT_COUNT_PER_TARGET)
                            : String.valueOf(2.0 * TargetCoverageSexGenotypeCalculator.MIN_ALLOWED_BAIT_COUNT_PER_TARGET);
                    return new Target(target.getName(), target.getInterval(), new TargetAnnotationCollection(
                            ImmutableMap.of(TargetAnnotation.BAIT_COUNT, baitCount)));

                })
                .collect(Collectors.toList());
        final Set<Target> badTargets = IntStream.range(0, numTargets)
                .filter(ti -> ti % 2 == 0)
                .mapToObj(targetsWithBadBaits::get)
                .collect(Collectors.toSet());

        final TargetCoverageSexGenotypeCalculator sexGenotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                readCountCollection, targetsWithBadBaits,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);
        final int[] autosomalTargetMasks = sexGenotypeCalculator.getAutosomalTargetMasks();
        final int[] allosomalTargetMasks = sexGenotypeCalculator.getAllosomalTargetMasks();

        /* generate the expected mask */
        final int numAutosomalTargets = (int)RANDOM_TARGETS_WITH_BAIT_COUNTS.stream()
                .filter(target -> AUTOSOMAL_CONTIGS.contains(target.getContig()))
                .count();
        final int numAllosomalTargets = (int)RANDOM_TARGETS_WITH_BAIT_COUNTS.stream()
                .filter(target -> ALLOSOMAL_CONTIGS.contains(target.getContig()))
                .count();
        final int[] expectedAutosomalTargetMasks = new int[numAutosomalTargets];
        Arrays.fill(expectedAutosomalTargetMasks, 1);
        final int[] expectedAllosomalTargetMasks = new int[numAllosomalTargets];
        Arrays.fill(expectedAllosomalTargetMasks, 1);
        IntStream.range(0, numAutosomalTargets)
                .filter(ti -> badTargets.contains(sexGenotypeCalculator.getAutosomalTargetList().get(ti)))
                .forEach(ti -> expectedAutosomalTargetMasks[ti] = 0);
        IntStream.range(0, numAllosomalTargets)
                .filter(ti -> badTargets.contains(sexGenotypeCalculator.getAllosomalTargetList().get(ti)))
                .forEach(ti -> expectedAllosomalTargetMasks[ti] = 0);

        ArrayAsserts.assertArrayEquals(autosomalTargetMasks, expectedAutosomalTargetMasks);
        ArrayAsserts.assertArrayEquals(allosomalTargetMasks, expectedAllosomalTargetMasks);
    }

    /**
     * Without bait counts, all biases must be 1.0
     */
    @Test
    public void testNoBiasWithoutBaitCounts() {
        final TargetCoverageSexGenotypeCalculator sexGenotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                RANDOM_READ_COUNTS_WITHOUT_BAIT_BIAS, RANDOM_TARGETS_WITHOUT_BAIT_COUNTS,
                HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS, null, MAPPING_ERROR_PROBABILITY);
        final double[] expectedAutosomalBiases = new double[sexGenotypeCalculator.getAutosomalTargetList().size()];
        Arrays.fill(expectedAutosomalBiases, 1.0);
        final double[] expectedAllosomalBiases = new double[sexGenotypeCalculator.getAllosomalTargetList().size()];
        Arrays.fill(expectedAllosomalBiases, 1.0);
        ArrayAsserts.assertArrayEquals(expectedAutosomalBiases, sexGenotypeCalculator.getAutosomalTargetBiases(), EPSILON);
        ArrayAsserts.assertArrayEquals(expectedAllosomalBiases, sexGenotypeCalculator.getAllosomalTargetBiases(), EPSILON);
    }

    /**
     * With bait counts, biases must be proportional to bait count
     */
    @Test
    public void testBiasWithBaitCounts() {
        /* generate a read count collection with no fully uncovered targets */
        final int numTargets = RANDOM_TARGETS_WITH_BAIT_COUNTS.size();
        final RealMatrix readCountsMatrix = new BlockRealMatrix(numTargets, NUMBER_OF_RANDOM_SAMPLES);
        final double[] data = new double[numTargets];
        Arrays.fill(data, 1.0);
        for (int si = 0; si < NUMBER_OF_RANDOM_SAMPLES; si++) {
            readCountsMatrix.setColumn(si, data);
        }
        final ReadCountCollection readCountCollection = new ReadCountCollection(RANDOM_TARGETS_WITH_BAIT_COUNTS,
                RANDOM_READ_COUNTS_WITH_BAIT_BIAS.columnNames(), readCountsMatrix);
        final TargetCoverageSexGenotypeCalculator sexGenotypeCalculator = new TargetCoverageSexGenotypeCalculator(
                readCountCollection, RANDOM_TARGETS_WITH_BAIT_COUNTS, HOMO_SAPIENS_GERMLINE_CONTIG_PLOIDY_ANNOTATIONS,
                null, MAPPING_ERROR_PROBABILITY);

        /* generate expected biases */
        final double[] expectedAutosomalBiases = sexGenotypeCalculator.getAutosomalTargetList().stream()
                .mapToDouble(target -> target.getAnnotations().getDouble(TargetAnnotation.BAIT_COUNT))
                .toArray();
        final double[] expectedAllosomalBiases = sexGenotypeCalculator.getAllosomalTargetList().stream()
                .mapToDouble(target -> target.getAnnotations().getDouble(TargetAnnotation.BAIT_COUNT))
                .toArray();
        final double meanBias = (Arrays.stream(expectedAutosomalBiases).sum() + Arrays.stream(expectedAllosomalBiases).sum()) /
                (expectedAutosomalBiases.length + expectedAllosomalBiases.length);

        IntStream.range(0, expectedAutosomalBiases.length).forEach(ti -> expectedAutosomalBiases[ti] /= meanBias);
        IntStream.range(0, expectedAllosomalBiases.length).forEach(ti -> expectedAllosomalBiases[ti] /= meanBias);

        ArrayAsserts.assertArrayEquals(expectedAutosomalBiases, sexGenotypeCalculator.getAutosomalTargetBiases(), EPSILON);
        ArrayAsserts.assertArrayEquals(expectedAllosomalBiases, sexGenotypeCalculator.getAllosomalTargetBiases(), EPSILON);
    }
}
