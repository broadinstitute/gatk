package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.TargetCodec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Integration tests for {@link AnnotateTargets}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AnnotateTargetsIntegrationTest extends CommandLineProgramTest {

    private static final File REFERENCE = new File(BaseTest.hg19MiniReference);

    // Test meta-parameters:
    private static final int MIN_TARGET_SIZE = 10;
    private static final int MAX_TARGET_SIZE = 1000;
    private static final int MEAN_TARGET_SIZE = 50;
    private static final double TARGET_SIZE_STDEV = MEAN_TARGET_SIZE;
    private static final int NUMBER_OF_TARGETS = 100;

    private static final int MIN_BAIT_SIZE = 5;
    private static final int MAX_BAIT_SIZE = 500;
    private static final int MEAN_BAIT_SIZE = 20;
    private static final double BAIT_SIZE_STDEV = MEAN_BAIT_SIZE;
    private static final int NUMBER_OF_BAITS = 1000; // must be greater than NUMBER_OF_TARGETS
    // End of meta-parameters

    private static final Random RANDOM = new Random(1313);
    private static final File TARGET_FILE = createTempFile("ati-target-file", ".tsv");
    private static final File TARGET_FILE_IDX = Tribble.indexFile(TARGET_FILE);

    private static final File BAIT_FILE = createTempFile("ati-bait-file", ".tsv");
    private static final File BAIT_FILE_IDX = Tribble.indexFile(BAIT_FILE);

    private static Map<Target, Double> EXPECTED_BAIT_COUNTS;
    private static final double EPSILON = 1e-12;

    @BeforeClass
    public void createRandomTargetAndBaitFiles() throws IOException  {
        final SAMSequenceDictionary referenceDictionary = resolveReferenceDictionary();
        final List<SimpleInterval> targetIntervals = createRandomIntervals(referenceDictionary, NUMBER_OF_TARGETS,
                MIN_TARGET_SIZE, MAX_TARGET_SIZE, MEAN_TARGET_SIZE, TARGET_SIZE_STDEV);
        final List<SimpleInterval> baitIntervals = Stream.concat(
                createRandomIntervals(referenceDictionary, NUMBER_OF_BAITS - NUMBER_OF_TARGETS, MIN_BAIT_SIZE,
                        MAX_BAIT_SIZE, MEAN_BAIT_SIZE, BAIT_SIZE_STDEV).stream(),
                createOverlappingIntervals(referenceDictionary, targetIntervals, MIN_BAIT_SIZE, MAX_BAIT_SIZE,
                        MEAN_BAIT_SIZE, BAIT_SIZE_STDEV).stream())
                .sorted(Comparator.comparing(SimpleInterval::getContig,
                        Comparator.comparingInt(referenceDictionary::getSequenceIndex))
                        .thenComparingInt(SimpleInterval::getStart)
                        .thenComparingInt(SimpleInterval::getEnd))
                .collect(Collectors.toList());

        /* index and write targets */
        final List<Target> targets = targetIntervals.stream().map(Target::new).collect(Collectors.toList());
        TargetWriter.writeTargetsToFile(TARGET_FILE, targets);
        final Index targetIndex = IndexFactory.createIndex(TARGET_FILE, new TargetCodec(), IndexFactory.IndexType.LINEAR);
        final LittleEndianOutputStream targetIndexStream = new LittleEndianOutputStream(new FileOutputStream(TARGET_FILE_IDX));
        targetIndex.write(targetIndexStream);
        targetIndexStream.close();

        /* index and write baits */
        final List<Target> baits = baitIntervals.stream().map(Target::new).collect(Collectors.toList());
        TargetWriter.writeTargetsToFile(BAIT_FILE, baits);
        final Index baitIndex = IndexFactory.createIndex(BAIT_FILE, new TargetCodec(), IndexFactory.IndexType.LINEAR);
        final LittleEndianOutputStream baitIndexStream = new LittleEndianOutputStream(new FileOutputStream(BAIT_FILE_IDX));
        baitIndex.write(baitIndexStream);
        baitIndexStream.close();

        /* calculate bait counts */
        final List<Set<SimpleInterval>> overlappingBaitsPerTargets = targets.stream()
                .map(target -> baits.stream()
                        .filter(bait -> bait.getInterval().overlaps(target.getInterval()))
                        .map(Target::getInterval)
                        .collect(Collectors.toSet()))
                .collect(Collectors.toList());
        final double[] baitCountsArray = AnnotateTargets.calculateBaitCounts(targets, overlappingBaitsPerTargets);
        EXPECTED_BAIT_COUNTS = IntStream.range(0, targets.size())
                .mapToObj(ti -> ImmutablePair.of(targets.get(ti), baitCountsArray[ti]))
                .collect(Collectors.toMap(Pair<Target, Double>::getKey, Pair<Target, Double>::getValue));
    }

    @AfterClass
    public void disposeTargetFile() {
        TARGET_FILE.delete();
        TARGET_FILE_IDX.delete();
        BAIT_FILE.delete();
        BAIT_FILE_IDX.delete();
    }

    private List<SimpleInterval> createRandomIntervals(final SAMSequenceDictionary referenceDictionary,
                                                       final int numberOfIntervals, final int minIntervalSize,
                                                       final int maxIntervalSize, final int meanIntervalSize,
                                                       final double intervalSizeStdev) {
        final List<SimpleInterval> result = new ArrayList<>(numberOfIntervals);
        final int numberOfSequences = referenceDictionary.getSequences().size();
        for (int i = 0; i < numberOfIntervals; i++) {
            final SAMSequenceRecord contig = referenceDictionary.getSequence(RANDOM.nextInt(numberOfSequences));
            final String contigName = contig.getSequenceName();
            final int intervalSize = Math.min(maxIntervalSize, (int) Math.max(minIntervalSize,
                    Math.round(RANDOM.nextDouble() * intervalSizeStdev + meanIntervalSize)));
            final int intervalStart = 1 + RANDOM.nextInt(contig.getSequenceLength() - intervalSize);
            final int intervalEnd = intervalStart + intervalSize - 1;
            final SimpleInterval interval = new SimpleInterval(contigName, intervalStart, intervalEnd);
            result.add(interval);
        }

        final Comparator<SimpleInterval> comparator = Comparator.comparing(SimpleInterval::getContig,
                Comparator.comparingInt(referenceDictionary::getSequenceIndex))
                .thenComparingInt(SimpleInterval::getStart)
                .thenComparingInt(SimpleInterval::getEnd);
        result.sort(comparator);
        return result;
    }

    private List<SimpleInterval> createOverlappingIntervals(final SAMSequenceDictionary referenceDictionary,
                                                            final List<SimpleInterval> originalIntervals,
                                                            final int minIntervalSize, final int maxIntervalSize,
                                                            final int meanIntervalSize, final double intervalSizeStdev) {
        final int numberOfIntervals = originalIntervals.size();
        final List<SimpleInterval> result = new ArrayList<>(numberOfIntervals);
        for (int i = 0; i < numberOfIntervals; i++) {
            final SimpleInterval originalInterval = originalIntervals.get(i);
            final String contigName = originalInterval.getContig();
            final int midPoint = (int)(0.5 * (originalInterval.getEnd() + originalInterval.getEnd()));
            final int intervalSize = Math.min(maxIntervalSize, (int) Math.max(minIntervalSize,
                    Math.round(RANDOM.nextDouble() * intervalSizeStdev + meanIntervalSize)));
            final int intervalStart = Math.max(1, midPoint - (int)(0.5 * intervalSize));
            final int intervalEnd = Math.max(intervalStart, midPoint + (int)(0.5 * intervalSize));
            final SimpleInterval interval = new SimpleInterval(contigName, intervalStart, intervalEnd);
            result.add(interval);
        }

        final Comparator<SimpleInterval> comparator = Comparator.comparing(SimpleInterval::getContig,
                Comparator.comparingInt(referenceDictionary::getSequenceIndex))
                .thenComparingInt(SimpleInterval::getStart)
                .thenComparingInt(SimpleInterval::getEnd);
        result.sort(comparator);
        return result;
    }

    private static SAMSequenceDictionary resolveReferenceDictionary() {
        return new ReferenceFileSource(REFERENCE).getSequenceDictionary();
    }

    @Override
    public String getTestedClassName() {
        return AnnotateTargets.class.getSimpleName();
    }

    @Test()
    public void testAllAnnotations() throws IOException {
        final File outputFile = createTempFile("ati-test-out-",".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.addAll(baseArguments(outputFile));
        arguments.addAll(gcContentDependenciesArguments());
        arguments.addAll(baitCountDependenciesArguments());
        runCommandLine(arguments);
        checkOutputFileContent(outputFile, true, true, true);
    }

    @Test()
    public void testGCContentAnnotationOnly() throws IOException {
        final File outputFile = createTempFile("ati-test-out-",".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.addAll(baseArguments(outputFile));
        arguments.addAll(gcContentDependenciesArguments());
        runCommandLine(arguments);
        checkOutputFileContent(outputFile, true, false, false);
    }

    @Test()
    public void testBaitCountAnnotationOnly() throws IOException {
        final File outputFile = createTempFile("ati-test-out-",".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.addAll(baseArguments(outputFile));
        arguments.addAll(baitCountDependenciesArguments());
        runCommandLine(arguments);
        checkOutputFileContent(outputFile, false, true, false);
    }

    private void checkOutputFileContent(final File outputFile,
                                        final boolean mustHaveGCContent,
                                        final boolean mustHaveBaitCount,
                                        final boolean mustHaveRepeats) throws IOException {
        try (final TargetTableReader outputReader = new TargetTableReader(outputFile);
             final TargetTableReader inputReader = new TargetTableReader(TARGET_FILE);
             final ReferenceFileSource reference = new ReferenceFileSource(REFERENCE)) {
            final SAMSequenceDictionary dictionary = reference.getSequenceDictionary();
            Target outputTarget;
            Target inputTarget;
            do {
                outputTarget = outputReader.readRecord();
                inputTarget = inputReader.readRecord();
                if (outputTarget == inputTarget) {
                    Assert.assertNull(outputTarget);
                    break;
                }
                Assert.assertNotNull(outputTarget, "too few targets in output");
                Assert.assertNotNull(inputTarget, "too many targets in output");
                Assert.assertEquals(outputTarget.getName(), inputTarget.getName());
                Assert.assertEquals(outputTarget.getInterval(), inputTarget.getInterval());
                final TargetAnnotationCollection annotations = outputTarget.getAnnotations();

                if (mustHaveGCContent) {
                    Assert.assertTrue(annotations.hasAnnotation(TargetAnnotation.GC_CONTENT));
                    checkOutputGCContent(reference, outputTarget.getInterval(), annotations.getDouble(TargetAnnotation.GC_CONTENT));
                } else {
                    Assert.assertFalse(annotations.hasAnnotation(TargetAnnotation.GC_CONTENT));
                }

                if (mustHaveBaitCount) {
                    Assert.assertTrue(annotations.hasAnnotation(TargetAnnotation.BAIT_COUNT));
                    checkOutputBaitCount(inputTarget, annotations.getDouble(TargetAnnotation.BAIT_COUNT));
                } else {
                    Assert.assertFalse(annotations.hasAnnotation(TargetAnnotation.BAIT_COUNT));
                }
            } while (true);
        }
    }

    private void checkOutputGCContent(final ReferenceFileSource reference, final SimpleInterval interval, final double observed) {
        final ReferenceSequence sequence = reference.queryAndPrefetch(interval);
        int total = 0;
        int gc = 0;
        for (final byte base : sequence.getBases()) {
            switch (Character.toLowerCase(base)) {
                case 'g':
                case 'c':
                    gc++; total++; break;
                case 'a':
                case 't':
                case 'u':
                    total++; break;
                default:
            }
        }
        final double expectedGCContent = total == 0 ? Double.NaN : ((double) gc / (double) total);
        if (Double.isNaN(expectedGCContent)) {
            Assert.assertTrue(Double.isNaN(observed));
        } else {
            Assert.assertEquals(observed, expectedGCContent, 0.0001);
        }
    }

    private void checkOutputBaitCount(final Target target, final double observed) {
        Assert.assertEquals(observed, EXPECTED_BAIT_COUNTS.get(target), EPSILON);
    }

    private Collection<? extends String> gcContentDependenciesArguments() {
        return Arrays.asList("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE.getPath());
    }

    private Collection<? extends String> baitCountDependenciesArguments() {
        return Arrays.asList("-" + ExomeStandardArgumentDefinitions.BAITS_FILE_SHORT_NAME, BAIT_FILE.getPath());
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoAnnotations() {
        final File outputFile = createTempFile("ati-test-out-",".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.addAll(baseArguments(outputFile));
        runCommandLine(arguments);
    }

    private List<String> baseArguments(final File outputFile) {
        return Arrays.asList("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                outputFile.getAbsolutePath(),
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME,
                TARGET_FILE.getAbsolutePath());
    }
}
