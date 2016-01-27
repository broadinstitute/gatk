package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModel;
import org.broadinstitute.hellbender.tools.exome.hmm.CopyNumberTriStateHiddenMarkovModelArgumentCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Integration tests for {@link DiscoverCopyNumberTriStateSegments}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class DiscoverCopyNumberTriStateIntegrationTest extends CommandLineProgramTest {

    private static final File REALISTIC_TARGETS_FILE = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/hmm/realistic-targets.tab");

    private static final TargetCollection<Target> REALISTIC_TARGETS;

    static {
        try (final TargetTableReader reader = new TargetTableReader(REALISTIC_TARGETS_FILE)) {
            REALISTIC_TARGETS = new HashedListTargetCollection<>(reader.stream().collect(Collectors.toList()));
        } catch (final IOException ex) {
            throw new GATKException("could not read the realistic-targets file", ex);
        }
    }

    private static final CopyNumberTriStateHiddenMarkovModel[] TEST_MODELS = new CopyNumberTriStateHiddenMarkovModel[] {
        new CopyNumberTriStateHiddenMarkovModel(1e-4, 6, -3, 3, 70_000)
    };

    /**
     * Checks that some called segment fits exactly between one or two XHMM called events.
     * <p>
     *     XHMM does not output diploid segments, so when we find out in the discovery output we
     *     make sure that it fits exactly between two events declared by XHMM.
     * </p>
     * @param before the XHMM event before (null if {@code between} is not preceeded by any event in the sample.
     * @param between the diploid call to assess.
     * @param after the XHMM event that follows {@code between}. If null indicates that {@code between} is the last segment of the sample.
     * @param targets the target collection.
     */
    public static void assertsItFitsInNicely(final XHMMOutputRecord before, final CopyNumberTriStateSegment between, final XHMMOutputRecord after,
                                             final TargetCollection<Target> targets, final String errorInfo) {
        Utils.nonNull(between);
        Assert.assertEquals(between.getCall(), CopyNumberTriState.NEUTRAL);
        final IndexRange betweenIndexRange = targets.indexRange(between.getInterval());
        if (before != null) {
            final IndexRange beforeIndexRange = targets.indexRange(before.interval);
            Assert.assertTrue(betweenIndexRange.from >= beforeIndexRange.to);
        } else {
            Assert.assertTrue(betweenIndexRange.from == 0 || !targets.target(betweenIndexRange.from - 1).getContig().equals(between.getInterval().getContig()));
        }
        if (after != null) {
            final IndexRange afterIndexRange = targets.indexRange(after.interval);
            Assert.assertTrue(betweenIndexRange.to <= afterIndexRange.from, " " + between.getInterval() + " " + after.interval + " error-info: " + errorInfo);
        } else {
            Assert.assertTrue(betweenIndexRange.to == targets.targetCount() || !targets.target(betweenIndexRange.to + 1).getContig().equals(between.getContig()));
        }
    }

    @Override
    public String getTestedClassName() {
        return DiscoverCopyNumberTriStateSegments.class.getSimpleName();
    }

    @Test(dataProvider = "testBadModelArgumentsData", expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadModelArguments(final String argumentShortName, final double value) {
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + argumentShortName);
        arguments.add(String.valueOf(value));
        final File inputFile = createTempFile("input", ".tab");
        Assert.assertTrue(inputFile.delete());
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.getAbsolutePath());
        final File outputFile = createTempFile("output", ".tab");
        Assert.assertTrue(outputFile.delete());
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments.toArray(new String[arguments.size()]));
    }

    private List<CopyNumberTriStateSegmentRecord> readOutputRecords(final File outputFile) {
        try (final CopyNumberTriStateSegmentRecordReader reader = new CopyNumberTriStateSegmentRecordReader(outputFile)) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException ex) {
            Assert.fail("problems reading the output file " + outputFile);
            throw new RuntimeException(ex);
        }
    }

    @Test(dataProvider = "simulatedChainData")
    public void testXhmmConcordance(final HiddenMarkovModelChain chain, final File xhmmOutputFile) {
        final File inputFile = writeChainInTempFile(chain);
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(inputFile.getAbsolutePath());
        //final File outputFile = new File("/tmp/333");
        final File outputFile = createTempFile("output", ".tab");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        // The model arguments:
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getDeletionMean()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getDuplicationMean()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.EVENT_RATE_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getEventRate()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getMeanEventTargetDistance()));
        arguments.add("-" + CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_TARGETS_PER_SEGMENT_SHORT_NAME);
        arguments.add(String.valueOf(chain.model.getNumberOfTargetsInEvent()));
        arguments.add("-" + DiscoverCopyNumberTriStateSegments.INPUT_ARE_ZSCORES_SHORT_NAME);
        runCommandLine(arguments.toArray(new String[arguments.size()]));
        Assert.assertTrue(outputFile.exists());
        final TargetCollection<Target> targets = TargetArgumentCollection.readTargetCollection(REALISTIC_TARGETS_FILE);
        final List<CopyNumberTriStateSegmentRecord> outputRecords = readOutputRecords(outputFile);
        assertOutputIsInOrder(outputRecords, targets);
        assertOutputHasConsistentNumberOfTargets(outputRecords, targets);
        final Map<String, List<CopyNumberTriStateSegmentRecord>> outputRecordsBySample = splitOutputRecordBySample(outputRecords);
        assertSampleNames(outputRecordsBySample.keySet(), chain);
        for (final List<CopyNumberTriStateSegmentRecord> sampleRecords : outputRecordsBySample.values()) {
            assertSampleSegmentsCoverAllTargets(sampleRecords, targets);
            assertSampleSegmentsCoordinates(sampleRecords, targets);
        }
        final List<XHMMOutputRecord> xhmmOutputRecords = readXhmmOutputRecords(xhmmOutputFile);
        final Map<String, List<XHMMOutputRecord>> xhmmOutputRecordsBySample = xhmmOutputRecords.stream()
                .collect(Collectors.groupingBy(r -> r.sample));
        for (final String sample : outputRecordsBySample.keySet()) {
            assertConcordantCalls(outputRecordsBySample.get(sample), xhmmOutputRecordsBySample.get(sample), targets);
        }

    }

    private void assertConcordantCalls(final List<CopyNumberTriStateSegmentRecord> ownRecords,
                                       final List<XHMMOutputRecord> xhmmOutputRecords, final TargetCollection<Target> targets) {
        Collections.sort(xhmmOutputRecords, Comparator.comparingInt(x -> targets.indexRange(x.interval).from));
        Collections.sort(ownRecords, Comparator.comparingInt(x -> targets.indexRange(x.getSegment().getInterval()).from));

        int nextOwn = 0; int nextXhmm = 0;
        XHMMOutputRecord xhmmPrevious = null;
        while (nextOwn < ownRecords.size()) {
            if (ownRecords.get(nextOwn).getSegment().getCall() == CopyNumberTriState.NEUTRAL) {
                assertsItFitsInNicely(xhmmPrevious, ownRecords.get(nextOwn++).getSegment(),
                        nextXhmm < xhmmOutputRecords.size() ? xhmmOutputRecords.get(nextXhmm) : null, targets, xhmmPrevious != null ?"" + xhmmPrevious.interval : "null");
            } else {
                xhmmOutputRecords.get(nextXhmm).assertComparable(ownRecords.get(nextOwn++));
                xhmmPrevious = xhmmOutputRecords.get(nextXhmm++);
            }
        }

    }

    private void assertSampleSegmentsCoordinates(List<CopyNumberTriStateSegmentRecord> sampleRecords, TargetCollection<Target> targets) {
        for (final CopyNumberTriStateSegmentRecord record : sampleRecords) {
            final IndexRange range = targets.indexRange(record.getSegment());
            Assert.assertTrue(range.size() > 0);
            Assert.assertEquals(record.getSegment().getContig(),targets.location(range.from).getContig());
            Assert.assertEquals(record.getSegment().getStart(), targets.location(range.from).getStart());
            Assert.assertEquals(record.getSegment().getEnd(), targets.location(range.to - 1).getEnd());
        }
    }

    private void assertSampleSegmentsCoverAllTargets(final List<CopyNumberTriStateSegmentRecord> sampleRecords, final TargetCollection<Target> targets) {
        int next = 0;
        for (final CopyNumberTriStateSegmentRecord record : sampleRecords) {
            final IndexRange range = targets.indexRange(record.getSegment());
            Assert.assertEquals(range.from, next);
            next = range.to;
        }
    }

    private void assertSampleNames(final Set<String> samples, final HiddenMarkovModelChain chain) {
        final int numberOfSamples = chain.data.size();
        final List<String> sampleNames = IntStream.range(0, chain.data.size()).mapToObj(a -> "SAMPLE_" + a).collect(Collectors.toList());
        Assert.assertEquals(samples.size(), numberOfSamples);
        for (final String sample : sampleNames) {
            Assert.assertTrue(samples.contains(sample));
        }
    }

    private Map<String,List<CopyNumberTriStateSegmentRecord>> splitOutputRecordBySample(final List<CopyNumberTriStateSegmentRecord> outputRecords) {
            return outputRecords.stream().collect(Collectors.groupingBy(CopyNumberTriStateSegmentRecord::getSampleName));
    }

    private void assertOutputIsInOrder(final List<CopyNumberTriStateSegmentRecord> outputRecords, final TargetCollection<Target> targets) {
        for (int i = 1; i < outputRecords.size(); i++) {
            final CopyNumberTriStateSegmentRecord nextRecord = outputRecords.get(i);
            final CopyNumberTriStateSegmentRecord previousRecord = outputRecords.get(i - 1);
            final IndexRange nextRange = targets.indexRange(nextRecord.getSegment());
            final IndexRange previousRange = targets.indexRange(previousRecord.getSegment());
            Assert.assertTrue(nextRange.from >= previousRange.from);
        }
    }

    private void assertOutputHasConsistentNumberOfTargets(final List<CopyNumberTriStateSegmentRecord> outputRecords, final TargetCollection<Target> targets) {
        for (final CopyNumberTriStateSegmentRecord nextRecord : outputRecords) {
            final IndexRange indexRange = targets.indexRange(nextRecord.getSegment());
            Assert.assertEquals(indexRange.to - indexRange.from, nextRecord.getSegment().getTargetCount());
        }
    }

    private List<XHMMOutputRecord> readXhmmOutputRecords(final File xhmmOutputFile) {
        try (final XHMMOutputReader reader = new XHMMOutputReader(xhmmOutputFile)) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException ex) {
            Assert.fail("could not xhmm output file: " + xhmmOutputFile);
            throw new RuntimeException(ex);
        }
    }

    private File writeChainInTempFile(final HiddenMarkovModelChain chain) {
        final File result = createTempFile("chain-",".tab");
        //final File result = new File("/tmp/input");
        final List<String> sampleNames = IntStream.range(0, chain.data.size()).mapToObj(a -> "SAMPLE_" + a).collect(Collectors.toList());
        final List<String> columnNames = new ArrayList<>(sampleNames.size() + 4);
        columnNames.addAll(Arrays.asList("CONTIG", "START", "END", "NAME"));
        columnNames.addAll(sampleNames);
        try (final TableWriter<Integer> writer = new TableWriter<Integer>(result,
                new TableColumnCollection(columnNames)) {
            @Override
            protected void composeLine(final Integer record, final DataLine dataLine) {
                dataLine.append(chain.targets.get(record).getContig())
                        .append(chain.targets.get(record).getStart())
                        .append(chain.targets.get(record).getEnd())
                        .append(chain.targets.get(record).getName());
                for (final List<Double> sampleData : chain.data) {
                    dataLine.append(sampleData.get(record).doubleValue());
                }
            }}) {

            writer.writeAllRecords(IntStream.range(0, chain.targets.size()).boxed().collect(Collectors.toList()));
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
        return result;
    }

    @DataProvider(name = "testBadModelArgumentsData")
    public Object[][] testBadModelArgumentsData() {
        return new Object[][] {
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.EVENT_RATE_SHORT_NAME, -1.0D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DELETION_COVERAGE_SHIFT_SHORT_NAME, 1.1D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DUPLICATION_COVERAGE_SHIFT_SHORT_NAME, -1.1D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_DISTANCE_BETWEEN_TARGETS_IN_SEGMENT_SHORT_NAME, -1.0D},
                {CopyNumberTriStateHiddenMarkovModelArgumentCollection.MEAN_TARGETS_PER_SEGMENT_SHORT_NAME, -2.0D},
        };
    }

    @Test
    public void testSimulateChainData() {
        final Random rdn = new Random(131313);
        final List<Object[]> result = Stream.of(TEST_MODELS)
                .map(m -> new Object[] { simulateChain(m, REALISTIC_TARGETS, 3, rdn) })
                .collect(Collectors.toList());
    }

    @DataProvider(name = "simulatedChainData")
    public Object[][] simulateChainData() {
        final Random rdn = new Random(131313);
        final List<Object[]> result = Stream.of(TEST_MODELS)
                .map(m -> new Object[] { simulateChain(m, REALISTIC_TARGETS, 3, rdn), new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/discover-germline-xhmm-output-4-6-70-3-3.tab") })
                .collect(Collectors.toList());
        return result.toArray(new Object[result.size()][]);
    }

    private HiddenMarkovModelChain simulateChain(final CopyNumberTriStateHiddenMarkovModel model, final TargetCollection<Target> targetCollection, final int sampleCount, final Random rdn) {
        final List<Target> targets = targetCollection.targets();
        final List<CopyNumberTriState> truth = new ArrayList<>(targets.size());
        final List<List<Double>> data = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++) {
            final List<Double> sampleData = randomSampleData(model, rdn, targets, truth);
            data.add(sampleData);
        }
        return new HiddenMarkovModelChain(model, targets, truth, data);
    }

    private List<Double> randomSampleData(CopyNumberTriStateHiddenMarkovModel model, Random rdn, List<Target> targets, List<CopyNumberTriState> truth) {
        final List<Double> sampleData = new ArrayList<>(targets.size());
        final Target firstPosition = targets.get(0);
        final CopyNumberTriState firstState = generateFirstState(model, firstPosition, rdn);
        truth.add(firstState);
        sampleData.add(model.randomDatum(firstState, rdn));
        CopyNumberTriState previousState = firstState;
        Target previousPosition = firstPosition;
        for (int i = 1; i < targets.size(); i++) {
            final Target position = targets.get(i);
            final CopyNumberTriState nextState = generateNextState(model, previousState, previousPosition, position, rdn);
            sampleData.add(model.randomDatum(nextState, rdn));
            truth.add(previousState = nextState);
            previousPosition = position;
        }
        return sampleData;
    }


    private class HiddenMarkovModelChain {

        public final List<Target> targets;
        public final List<CopyNumberTriState> truth;
        public final List<List<Double>> data;
        public final CopyNumberTriStateHiddenMarkovModel model;

        private HiddenMarkovModelChain(final CopyNumberTriStateHiddenMarkovModel model, final List<Target> targets, final List<CopyNumberTriState> truth, final List<List<Double>> data) {
            this.targets = targets;
            this.truth = truth;
            this.data = data;
            this.model = model;
        }
    }

    private CopyNumberTriState generateNextState(final CopyNumberTriStateHiddenMarkovModel model, final CopyNumberTriState previousState, final Target previousPosition, final Target nextPosition, final Random rdn) {
        final double[] logTransitionProbs = new double[CopyNumberTriState.values().length];
        for (int i = 0; i < logTransitionProbs.length; i++) {
            logTransitionProbs[i] = model.logTransitionProbability(previousState, previousPosition, CopyNumberTriState.values()[i], nextPosition);
        }
        final double total = DoubleStream.of(logTransitionProbs).map(Math::exp).sum();
        final double uniform = rdn.nextDouble() * total;
        double accumulative = 0;
        for (int i = 0; i < logTransitionProbs.length; i++) {
            accumulative += Math.exp(logTransitionProbs[i]);
            if (uniform < accumulative) {
                return CopyNumberTriState.values()[i];
            }
        }
        throw new IllegalStateException("may be some of the trans probs are negative");
    }

    private CopyNumberTriState generateFirstState(final CopyNumberTriStateHiddenMarkovModel model, final Target firstPosition, final Random rdn) {

        final double[] priors = new double[CopyNumberTriState.values().length];
        for (int i = 0; i < priors.length; i++) {
            priors[i] = model.logPriorProbability(CopyNumberTriState.values()[i], firstPosition);
        }
        final double uniform = rdn.nextDouble() * DoubleStream.of(priors).map(Math::exp).sum();

        double accumulative = 0;
        for (int i = 0; i < priors.length; i++) {
            accumulative += Math.exp(model.logPriorProbability(CopyNumberTriState.values()[i], firstPosition));
            if (uniform < accumulative) {
                return CopyNumberTriState.values()[i];
            }
        }
        throw new IllegalStateException("bug in testing code: seems that the test model priors sum to more than 1: " + Arrays.toString(priors));
    }
}
