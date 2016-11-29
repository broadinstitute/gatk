package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Integration tests for {@link FilterTargets}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class FilterTargetsIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return FilterTargets.class.getSimpleName();
    }

    @Test(dataProvider = "targetTargetSizeFilterData")
    public void testTargetSizeFilter(final int min, final int max, final List<Target> targets)
        throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getInterval().size() >= min && t.getInterval().size() <= max)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> t.getInterval().size() < min || t.getInterval().size() > max)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MINIMUM_TARGET_SIZE_SHORT_NAME, String.valueOf(min),
                "-" + FilterTargets.MAXIMUM_TARGET_SIZE_SHORT_NAME, String.valueOf(max),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExtremeTargetSize);
    }

    @Test(dataProvider = "targetExcludedIntervalsData")
    public void testExcludedIntervalsThroughFile(final List<String> excludeSpecs,
                                      final List<SimpleInterval> excludedIntervals,
                                      final List<Target> targets)
            throws IOException
    {
        final File excludedIntervalsFile = createExcludedIntervalsFile(excludeSpecs);
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> excludedIntervals.stream().noneMatch(e -> e.overlaps(t.getInterval())))
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> excludedIntervals.stream().anyMatch(e -> e.overlaps(t.getInterval())))
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        final List<String> commonArgument = new ArrayList<>(Arrays.asList("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath()));
        final List<String> arguments = new ArrayList<>(commonArgument);
        arguments.add("-" + FilterTargets.EXCLUDED_INTERVALS_SHORT_NAME);
        arguments.add(excludedIntervalsFile.getPath());
        runCommandLine(arguments);

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExcludedInterval);
    }

    private File createExcludedIntervalsFile(final List<String> excludeSpecs) throws IOException {
        final File result = createTempFile("excluded-interval",".list");
        try (final PrintWriter writer = new PrintWriter(new FileWriter(result))) {
            writer.println("# comment1");
            writer.println();
            for (final String s : excludeSpecs) {
                writer.println(s);
            }
        }
        return result;
    }

    @Test(dataProvider = "targetExcludedIntervalsData")
    public void testExcludedIntervals(final List<String> excludeSpecs,
                                      final List<SimpleInterval> excludedIntervals,
                                      final List<Target> targets)
            throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> excludedIntervals.stream().noneMatch(e -> e.overlaps(t.getInterval())))
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> excludedIntervals.stream().anyMatch(e -> e.overlaps(t.getInterval())))
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        final List<String> commonArgument = new ArrayList<>(Arrays.asList("-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath()));
        final List<String> arguments = new ArrayList<>(commonArgument);
        for (final String spec : excludeSpecs) {
            arguments.add("-" + FilterTargets.EXCLUDED_INTERVALS_SHORT_NAME);
            arguments.add(spec);
        }
        runCommandLine(arguments);

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExcludedInterval);
    }

    private void checkLeftInTargets(final List<Target> leftInTargets, final File outputFile) {
        final TargetCollection<Target> outputTargetCollection = TargetArgumentCollection.readTargetCollection(outputFile);
        final List<Target> outputTargets = outputTargetCollection.targets();
        Assert.assertEquals(outputTargets.size(), leftInTargets.size());
        Assert.assertEquals(new HashSet<>(outputTargets), new HashSet<>(leftInTargets));
        for (int i = 0; i < outputTargets.size(); i++) {
            Assert.assertEquals(outputTargets.get(i), leftInTargets.get(i));
            Assert.assertEquals(outputTargets.get(i).getInterval(),
                    leftInTargets.get(i).getInterval());
        }
    }

    @Test(dataProvider = "targetTargetSizeFilterData")
    public void testTargetSizeFilterWithoutRejectionFile(final int min, final int max, final List<Target> targets) {
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getInterval().size() >= min && t.getInterval().size() <= max)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_TARGET_SIZE_SHORT_NAME, String.valueOf(max),
                "-" + FilterTargets.MINIMUM_TARGET_SIZE_SHORT_NAME, String.valueOf(min),
        });
        checkLeftInTargets(leftInTargets, outputFile);
    }

    @Test(dataProvider = "targetExtremeGCFilterData")
    public void testTargetExtremeGCFilter(final List<Target> targets, final double extremeGC) throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.singleton(TargetAnnotation.GC_CONTENT));
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) >= extremeGC
                && t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) <= 1 - extremeGC)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) < extremeGC
                        && t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) > 1.0 - extremeGC)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_GC_CONTENT_SHORT_NAME, String.valueOf(1.0 - extremeGC),
                "-" + FilterTargets.MINIMUM_GC_CONTENT_SHORT_NAME, String.valueOf(extremeGC),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExtremeGCContent);
    }

    private void checkRejectedTargets(final List<Target> leftOutTargets, final File rejectionFile,
                                      final FilterTargets.TargetFilter filter) throws IOException {
        try (final TableReader<String[]> rejectionReader = new TableReader<String[]>(rejectionFile) {
            @Override
            protected String[] createRecord(DataLine dataLine) {
                return dataLine.toArray();
            }
        }) {
            Assert.assertTrue(rejectionReader.columns().containsAll(
                    FilterTargets.REJECT_FILE_REASON_COLUMN_NAME,
                    FilterTargets.REJECT_FILE_FILTER_COLUMN_NAME,
                    FilterTargets.REJECT_FILE_TARGET_COLUMN_NAME
            ));
            for (final Target target : leftOutTargets) {
                final String[] values = rejectionReader.readRecord();
                Assert.assertEquals(values[rejectionReader.columns().indexOf(FilterTargets.REJECT_FILE_TARGET_COLUMN_NAME)],
                        target.getName());
                Assert.assertEquals(values[rejectionReader.columns().indexOf(FilterTargets.REJECT_FILE_FILTER_COLUMN_NAME)],
                        filter.toString());
            }
        }
    }

    @Test(dataProvider = "targetExtremeVarianceData")
    public void testTargetExtremeCoverageVarianceFilterUsingAnnotationFile(final List<Target> targets, final double min, final double max)
            throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final File annotationFile = createTargetFile(targets, Collections.singleton(TargetAnnotation.COVERAGE_VARIANCE));
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.COVERAGE_VARIANCE) >= min
                        && t.getAnnotations().getDouble(TargetAnnotation.COVERAGE_VARIANCE) <= max)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.COVERAGE_VARIANCE) < min
                        || t.getAnnotations().getDouble(TargetAnnotation.COVERAGE_VARIANCE) > max)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + TargetArgumentCollection.TARGET_ANNOTATION_FILES_SHORT_NAME, annotationFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_COVERAGE_VARIANCE_SHORT_NAME, String.valueOf(max),
                "-" + FilterTargets.MINIMUM_COVERAGE_VARIANCE_SHORT_NAME, String.valueOf(min),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExtremeCoverageVariance);
    }

    @Test(dataProvider = "targetExtremeInterquartileRangeData")
    public void testTargetExtremeCoverageInterquartileRangeFilterUsingAnnotationFile(final List<Target> targets, final double max)
            throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final File annotationFile = createTargetFile(targets, Collections.singleton(TargetAnnotation.COVERAGE_INTERQUARTILE_RANGE));
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.COVERAGE_INTERQUARTILE_RANGE) <= max)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.COVERAGE_INTERQUARTILE_RANGE) > max)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + TargetArgumentCollection.TARGET_ANNOTATION_FILES_SHORT_NAME, annotationFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_COVERAGE_INTERQUARTILE_RANGE_SHORT_NAME, String.valueOf(max),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExtremeCoverageInterquartileRange);
    }

    @Test(dataProvider = "targetExtremeGCFilterData")
    public void testTargetExtremeGCFilterUsingAnnotationFile(final List<Target> targets, final double extremeGC)
            throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final File annotationFile = createTargetFile(targets, Collections.singleton(TargetAnnotation.GC_CONTENT));
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) >= extremeGC
                        && t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) <= 1 - extremeGC)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) < extremeGC
                        && t.getAnnotations().getDouble(TargetAnnotation.GC_CONTENT) > 1.0 - extremeGC)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + TargetArgumentCollection.TARGET_ANNOTATION_FILES_SHORT_NAME, annotationFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_GC_CONTENT_SHORT_NAME, String.valueOf(1.0 - extremeGC),
                "-" + FilterTargets.MINIMUM_GC_CONTENT_SHORT_NAME, String.valueOf(extremeGC),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExtremeGCContent);
    }

    @Test(dataProvider = "targetExtremeGCFilterData", expectedExceptions = UserException.BadInput.class)
    public void testTargetExtremeGCFilterMissingAnnotation(final List<Target> targets, final double extremeGC) {
        if (extremeGC <= 0.0) {
            throw new UserException.BadInput("when minGC is 0.0 or less we don't throw an exception");
        }
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[]{
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_GC_CONTENT_SHORT_NAME, String.valueOf(1.0 - extremeGC),
                "-" + FilterTargets.MINIMUM_GC_CONTENT_SHORT_NAME, String.valueOf(extremeGC),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });
    }

    @Test(dataProvider = "targetExtremeGCFilterData", expectedExceptions = UserException.BadInput.class)
    public void testTargetExtremeRepeatFilterMissingAnnotation(final List<Target> targets, final double maxRepeatContent) {
        if (maxRepeatContent >= 1.0) {
            throw new UserException.BadInput("when maxRepeat is 1.0 or greater we don't throw an exception");
        }
        final File targetFile = createTargetFile(targets, Collections.emptySet());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[]{
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_REPEAT_CONTENT_SHORT_NAME, String.valueOf(maxRepeatContent),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });
    }

    @Test(dataProvider = "targetExtremeRepeatFilterData")
    public void testTargetExtremeRepeatFilter(final List<Target> targets, final double maxRepeatContent)
            throws IOException
    {
        final File targetFile = createTargetFile(targets, Collections.singleton(TargetAnnotation.REPEAT_FRACTION));
        final List<Target> leftInTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.REPEAT_FRACTION) <= maxRepeatContent)
                .collect(Collectors.toList());
        final List<Target> leftOutTargets = targets.stream()
                .filter(t -> t.getAnnotations().getDouble(TargetAnnotation.REPEAT_FRACTION) > maxRepeatContent)
                .collect(Collectors.toList());
        final File outputFile = createTempFile("output", ".tab");
        final File rejectionFile = createTempFile("reject", ".tab");
        runCommandLine(new String[] {
                "-" + TargetArgumentCollection.TARGET_FILE_SHORT_NAME, targetFile.getPath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getPath(),
                "-" + FilterTargets.MAXIMUM_REPEAT_CONTENT_SHORT_NAME, String.valueOf(maxRepeatContent),
                "-" + FilterTargets.REJECT_OUTPUT_FILE_SHORT_NAME, rejectionFile.getPath(),
        });

        checkLeftInTargets(leftInTargets, outputFile);
        checkRejectedTargets(leftOutTargets, rejectionFile, FilterTargets.TargetFilter.ExtremeRepeatContent);
    }

    private File createTargetFile(final List<Target> targets, final Set<TargetAnnotation> annotations) {
        final File result = createTempFile("filter-test-target", ".tab");
        try (final TargetWriter writer = new TargetWriter(result, annotations)) {
            for (final Target target : targets.stream()
                    .sorted((a, b) -> IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR.compare(a.getInterval(), b.getInterval()))
                    .collect(Collectors.toList())) {
                writer.writeRecord(target);
            }
            return result;
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    @DataProvider(name = "targetExcludedIntervalsData")
    public Object[][] targetExcludedIntervalsData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval((i & 1) == 0 ?"chr1" : "chr2", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))))
                .collect(Collectors.toList());

        return new Object[][] {
                { Collections.singletonList("chr2:1-2000"),
                  Collections.singletonList(new SimpleInterval("chr2", 1, 2000)), targets},
                { Collections.singletonList("chr1"),
                  Collections.singletonList(new SimpleInterval("chr1")), targets},
                { Collections.singletonList("chr1:10000+"),
                  Collections.singletonList(new SimpleInterval("chr1:10000+")), targets},
                { Arrays.asList("chr2", "chr1:1-10000"),
                  Arrays.asList(new SimpleInterval("chr1:1-10000"), new SimpleInterval("chr2")), targets},
                { Arrays.asList("chr2:1003", "chr1:2010", "chr2:200000+", "chr1:500000-700000"),
                  Arrays.asList(new SimpleInterval("chr2:1003"), new SimpleInterval("chr1:2010"),
                                new SimpleInterval("chr2:200000+"),
                                new SimpleInterval("chr1:500000-700000")) , targets},
                { Collections.emptyList(), Collections.emptyList(), targets},
        };
    }

    @DataProvider(name = "targetTargetSizeFilterData")
    public Object[][] targetTargetSizeFilterData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval("chr1", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))))
                .collect(Collectors.toList());

        return new Object[][] {
                { -1, 10000, targets},
                { 0, Integer.MAX_VALUE, targets},
                { 1, 100, targets},
                { 10, 500, targets},
                { 100, 200, targets},
                { 500, 1001, targets},
        };
    }

    @DataProvider(name = "targetExtremeGCFilterData")
    public Object[][] targetExtremeGCFilterData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval("chr1", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))
                , new TargetAnnotationCollection(Collections.singletonMap(TargetAnnotation.GC_CONTENT, String.valueOf(rdn.nextDouble())))))
                .collect(Collectors.toList());

        return new Object[][] {
                { targets,  0},
                { targets,  0.1},
                { targets,  0.2},
                { targets,  0.49},
                { targets,  0.5},
        };
    }

    @DataProvider(name = "targetExtremeRepeatFilterData")
    public Object[][] targetExtremeRepeatFilterData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval("chr1", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))
                        , new TargetAnnotationCollection(Collections.singletonMap(TargetAnnotation.REPEAT_FRACTION, String.valueOf(rdn.nextDouble())))))
                .collect(Collectors.toList());

        return new Object[][] {
                { targets,  0},
                { targets,  0.1},
                { targets,  0.2},
                { targets,  0.49},
                { targets,  0.5},
                { targets,  1.0},
                { targets,  0.9},
                { targets,  -1.0},
        };
    }

    @DataProvider(name = "targetExtremeCoverageMeanData")
    public Object[][] targetExtremeCoverageMeanData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval("chr1", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))
                        , new TargetAnnotationCollection(Collections.singletonMap(TargetAnnotation.MEAN_COVERAGE, String.valueOf(rdn.nextDouble() * 1000)))))
                .collect(Collectors.toList());

        return new Object[][] {
                { targets,  0, 10},
                { targets,  10, 50},
                { targets,  500, 600},
                { targets,  -1.0, Double.POSITIVE_INFINITY},
                { targets, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY},
        };
    }

    @DataProvider(name = "targetExtremeVarianceData")
    public Object[][] targetExtremeCoverageVarianceData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval("chr1", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))
                        , new TargetAnnotationCollection(Collections.singletonMap(TargetAnnotation.COVERAGE_VARIANCE, String.valueOf(rdn.nextDouble() * 100)))))
                .collect(Collectors.toList());

        return new Object[][] {
                { targets,  0, 10},
                { targets,  10, 50},
                { targets,  90, 110},
                { targets,  -1.0, Double.POSITIVE_INFINITY},
                { targets, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY},
        };
    }

    @DataProvider(name = "targetExtremeInterquartileRangeData")
    public Object[][] targetExtremeCoverageInterquartileRangeData() {
        final Random rdn = new Random(1313);
        final List<Target> targets = IntStream.range(1, 1001)
                .mapToObj(i -> new Target("target_" + i,
                        new SimpleInterval("chr1", i * 1002, i * 1002 + 1 + rdn.nextInt(1000))
                        , new TargetAnnotationCollection(Collections.singletonMap(TargetAnnotation.COVERAGE_INTERQUARTILE_RANGE, String.valueOf(rdn.nextDouble() * 100)))))
                .collect(Collectors.toList());

        return new Object[][] {
                { targets,  10},
                { targets,  50},
                { targets,  110},
                { targets,  Double.POSITIVE_INFINITY},
        };
    }
}
