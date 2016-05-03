package org.broadinstitute.hellbender.utils.hmm;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Unit test for HMM api classes including {@link ViterbiAlgorithm}, {@link ForwardBackwardAlgorithm} and
 * {@link ForwardBackwardAlgorithm.Result}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HMMUnitTest extends BaseTest {

    private final static Random RANDOM = new Random(1313);

    private static List<TestHMModel> TEST_MODELS = new ArrayList<>(Arrays.asList(
            TestHMModel.fromPhredProbabilities(
                    20, 0, 20, // priors
                    0, 20, 20, // transition
                    20, 0, 20,
                    20, 20, 0,
                    0, 20, 20, // emission
                    20, 0, 20,
                    20, 20, 0),
            TestHMModel.fromPhredProbabilities(
                    0, 0, 0, // priors
                    0, .2, .2, // transition
                    .2, 0, .2,
                    .2, .2, 0,
                    0, .2, .2, // emission
                    .2, 0, .2,
                    .2, .2, 0),
            TestHMModel.fromPhredProbabilities(
                    0, 0, 0, // priors
                    0, 0, 0, // transition
                    0, 0, 0,
                    0, 0, 0,
                    0, 0, 0, // emission
                    0, 0, 0,
                    0, 0, 0),
            TestHMModel.fromPhredProbabilities(
                    30.1, 0.6, 10.2, // priors
                    0.7, 30.3, 10.4, // transition
                    10.5, 0.4, 30.6,
                    30.7, 10.8, 0.5,
                    0.3, 10.9, 20.1, // emission
                    40.2, 0.2, 10.3,
                    5.5, 5.4, 0.1)
    ));

    private static final int TEST_PATH_LENGTH = 200;

    private static List<Pair<List<TestHMModel.State>, List<TestHMModel.Datum>>> TEST_SEQUENCES;

    private static File TEST_SEQUENCE_FILE;

    private static File TEST_R_RESULTS_FILE;

    private static List<ExpectedResult> TEST_EXPECTED_RESULTS;

    @Test
    public void testSetup() throws IOException {
        setUp();
    }

    /**
     * This method runs the Test models and sequences against R's HMM package to
     * obtain independent estimates for the expected values for the Viterbi and Forward-Backward
     * algorithm.
     * <p>
     *     Typically this would be annotated as a <code>@BeforeClass</> method but unfortunately the data-provider
     *     code may be run before it and so we need to implement lazy initialization.
     * </p>
     * <p>
     *     Therefore the data-providers will call it once and only once before any of their data consumer test are run.
     * </p>
     *
     * @throws IOException
     */
    private void setUp() throws IOException {
        final List<Integer> positions = IntStream.range(0, TEST_PATH_LENGTH)
                .boxed()
                .collect(Collectors.toList());
        TEST_SEQUENCES = TEST_MODELS.stream().map(m -> m.generate(positions, RANDOM))
                .collect(Collectors.toList());

        TEST_SEQUENCE_FILE = createTempFile("sequences", ".seq");
        final PrintWriter writer = new PrintWriter(new FileWriter(TEST_SEQUENCE_FILE));
        TEST_SEQUENCES.forEach(s -> {
            final StringBuilder builder = new StringBuilder(TEST_PATH_LENGTH);
            s.getSecond().forEach(c -> {
                builder.append(c.toString());
                builder.append(',');
            });
            builder.setLength(builder.length() - 1);
            writer.println(builder.toString());
        });
        writer.close();
        TEST_R_RESULTS_FILE = createTempFile("r-output", ".tab");
        RScriptExecutor rExecutor = new RScriptExecutor();
        final File script = composeRScriptFile();
        rExecutor.addScript(script);
        if (!rExecutor.exec()) {
            Assert.fail("could not obtain expected results from R");
        }
        @SuppressWarnings({"rawtypes", "unchecked"})
        final List<RResultRecord>[][] recordsByModelAndSequence = (List<RResultRecord>[][]) new List[TEST_MODELS.size()][TEST_SEQUENCES.size()];
        try (final RResultReader reader = new RResultReader(TEST_R_RESULTS_FILE)) {
            reader.stream().forEach(r -> {
                if (recordsByModelAndSequence[r.modelIndex][r.sequenceIndex] == null) {
                    recordsByModelAndSequence[r.modelIndex][r.sequenceIndex] = new ArrayList<>(TEST_PATH_LENGTH);
                }
                recordsByModelAndSequence[r.modelIndex][r.sequenceIndex].add(r);
            });
        }
        TEST_EXPECTED_RESULTS = new ArrayList<>(TEST_MODELS.size() * TEST_SEQUENCES.size());
        for (int i = 0; i < recordsByModelAndSequence.length; i++) {
            for (int j = 0; j < recordsByModelAndSequence[i].length; j++) {
                TEST_EXPECTED_RESULTS.add(composeExpectedResult(recordsByModelAndSequence[i][j], i, j));
            }
        }
    }

    private ExpectedResult composeExpectedResult(final List<RResultRecord> rResultRecords, final int modelIndex, final int sequenceIndex) {
        final List<TestHMModel.State> bestPath = rResultRecords.stream()
                .sorted(Comparator.comparing(r -> r.dataIndex))
                .map(r -> r.bestPathState)
                .collect(Collectors.toList());
        final List<double[]> logForwardProbs = rResultRecords.stream()
                .sorted(Comparator.comparing(r -> r.dataIndex))
                .map(r -> r.logForwardProbs)
                .collect(Collectors.toList());
        final List<double[]> logBackwardProbs = rResultRecords.stream()
                .sorted(Comparator.comparing(r -> r.dataIndex))
                .map(r -> r.logBackwardProbs)
                .collect(Collectors.toList());
        final List<double[]> logProbabilities = rResultRecords.stream()
                .sorted(Comparator.comparing(r -> r.dataIndex))
                .map(r -> r.logProbability)
                .collect(Collectors.toList());

        return new ExpectedResult(TEST_MODELS.get(modelIndex), TEST_SEQUENCES.get(sequenceIndex).getSecond(), bestPath,
                logForwardProbs, logBackwardProbs, logProbabilities);
    }

    private File composeRScriptFile() throws IOException {
        final File script = createTempFile("r-script", ".R");
        try (final PrintWriter scriptWriter = new PrintWriter(new FileWriter(script))) {
            scriptWriter.println(TestHMModel.toHMMInstallRString());
            scriptWriter.println("outfile = \"" + TEST_R_RESULTS_FILE.getAbsolutePath() + '"');
            scriptWriter.println("sequences = strsplit(readLines(\"" + TEST_SEQUENCE_FILE.getPath() + "\"), ',')");
            scriptWriter.println("sequences.n = length(sequences)");
            scriptWriter.println();
            scriptWriter.println("models = " + TestHMModel.toHMMModelDeclarationRString(TEST_MODELS));
            scriptWriter.println("models.n = length(models)");
            scriptWriter.println("write(x = paste('SEQUENCE', 'MODEL', 'POSITION', 'BEST_PATH'," +
                    " 'FW_A', 'FW_B', 'FW_C'," +
                    " 'BW_A', 'BW_B', 'BW_C'," +
                    " 'PP_A', 'PP_B', 'PP_C', sep ='\\t'), file = outfile, append = F)");
            scriptWriter.println("for (i in 1:models.n) {");
            scriptWriter.println("    for (j in 1:sequences.n) {");
            scriptWriter.println("        best_path = viterbi(models[[i]], sequences[[j]])");
            scriptWriter.println("        fw = forward(models[[i]], sequences[[j]])");
            scriptWriter.println("        bw = backward(models[[i]], sequences[[j]])");
            scriptWriter.println("        pp = log(posterior(models[[i]], sequences[[j]]))");
            scriptWriter.println("        for (k in 1:length(best_path)) {");
            scriptWriter.println("            write(x = paste(j - 1, i - 1, k - 1, best_path[k]," +
                    " fw['A', k], fw['B', k], fw['C', k]," +
                    " bw['A', k], bw['B', k], bw['C', k]," +
                    " pp['A', k], pp['B', k], pp['C', k], sep = '\\t'), file = outfile, append = T)");
            scriptWriter.println("}}}");
        }
        return script;
    }

    @SuppressWarnings("all")
    @AfterClass
    private void tearDown() {
        if (TEST_SEQUENCE_FILE != null) {
            TEST_SEQUENCE_FILE.delete();
        }
        if (TEST_R_RESULTS_FILE != null) {
            TEST_R_RESULTS_FILE.delete();
        }
    }

    @Test(dataProvider = "testViterbiData")
    public void testViterbi(final TestHMModel model, final List<TestHMModel.Datum> sequence, final List<TestHMModel.State> expected) {
        final List<Integer> positions = IntStream.range(0, sequence.size()).boxed().collect(Collectors.toList());
        final List<TestHMModel.State> observed = ViterbiAlgorithm.apply(sequence, positions, model);
        try {
            Assert.assertEquals(observed, expected);
        } catch (final AssertionError ex) {
            // If not the same sequence of states, we need to show that they have the same probability.
            final ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State>
                    fb = ForwardBackwardAlgorithm.apply(sequence, positions, model);
            final double observedProb = fb.logProbability(observed);
            final double expectedProb = fb.logProbability(expected);
            Assert.assertEquals(observedProb, expectedProb, 0.0000001);
        }
    }

    /**
     * Test {@link ForwardBackwardAlgorithm.Result} methods not covered by other tests.
     */
    @Test(dataProvider = "testFBResultData")
    public void testFBResult(final TestHMModel model, final List<TestHMModel.Datum> sequence) {
        final List<Integer> positions = IntStream.range(0, sequence.size()).boxed().collect(Collectors.toList());
        final ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> forwardBackwardResult
                = ForwardBackwardAlgorithm.apply(sequence, positions, model);
        Assert.assertNotNull(forwardBackwardResult);
        Assert.assertEquals(forwardBackwardResult.data(), sequence);
        Assert.assertEquals(forwardBackwardResult.positions(), positions);
        Assert.assertSame(forwardBackwardResult.model(), model);
    }

    @Test(dataProvider = "testForwardData")
    public void testForward(final TestHMModel model, final List<TestHMModel.Datum> sequence, final List<double[]> expected) {
        final List<Integer> positions = IntStream.range(0, sequence.size()).boxed().collect(Collectors.toList());
        final ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> forwardBackwardResult
                = ForwardBackwardAlgorithm.apply(sequence, positions, model);
        Assert.assertNotNull(forwardBackwardResult);
        for (final Integer position : positions) {
            for (final TestHMModel.State state : TestHMModel.State.values()) {
                final double expectedValue = expected.get(position)[state.ordinal()];
                final double observedValue = forwardBackwardResult.logForwardProbability(position, state);
                final double observedValueUsingIndex = forwardBackwardResult.logForwardProbability(position.intValue(), state);
                Assert.assertEquals(observedValueUsingIndex, observedValue);
                final double epsilon = 0.001 * Math.min(Math.abs(expectedValue),
                        Math.abs(observedValue));
                Assert.assertEquals(observedValue, expectedValue, epsilon);
            }
        }
    }

    @Test(dataProvider = "testBackwardData")
    public void testBackward(final TestHMModel model, final List<TestHMModel.Datum> sequence, final List<double[]> expected) {
        final List<Integer> positions = IntStream.range(0, sequence.size()).boxed().collect(Collectors.toList());
        final ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> forwardBackwardResult
                = ForwardBackwardAlgorithm.apply(sequence, positions, model);
        Assert.assertNotNull(forwardBackwardResult);
        for (int position = positions.size() - 1; position >= 0; position--) {
            for (final TestHMModel.State state : TestHMModel.State.values()) {
                final double expectedValue = expected.get(position)[state.ordinal()];
                final double observedValue = forwardBackwardResult.logBackwardProbability(new Integer(position), state);
                final double observedValueUsingIndex = forwardBackwardResult.logBackwardProbability(position, state);
                Assert.assertEquals(observedValueUsingIndex, observedValue);
                final double epsilon = 0.001 * Math.min(Math.abs(expectedValue),
                        Math.abs(observedValue));
                Assert.assertEquals(observedValue, expectedValue, epsilon, " " + state + " " + position);
            }
        }
    }

    @Test(dataProvider = "testPosteriorData")
    public void testPosterior(final TestHMModel model, final List<TestHMModel.Datum> sequence, final List<double[]> expected) {
        final List<Integer> positions = IntStream.range(0, sequence.size()).boxed().collect(Collectors.toList());
        final ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> forwardBackwardResult
                = ForwardBackwardAlgorithm.apply(sequence, positions, model);
        Assert.assertNotNull(forwardBackwardResult);
        for (final Integer position : positions) {
            for (final TestHMModel.State state : TestHMModel.State.values()) {
                final double expectedValue = expected.get(position)[state.ordinal()];
                final double observedValue = forwardBackwardResult.logProbability(position, state);
                final double observedValueUsingIndex = forwardBackwardResult.logProbability(position.intValue(), state);
                Assert.assertEquals(observedValue, observedValueUsingIndex);
                // When probs are close to 1, is better to compare them after exp. them.
                if (Math.exp(observedValue) > .5) {
                    final double epsilon = 0.001 * Math.min(Math.exp(observedValue),
                            Math.exp(expectedValue));
                    Assert.assertEquals(Math.exp(observedValue), Math.exp(expectedValue), epsilon);
                } else {
                    final double epsilon = 0.001 * Math.min(Math.abs(observedValue),
                            Math.abs(expectedValue));
                    Assert.assertEquals(observedValue, expectedValue, epsilon);
                }
            }
        }
    }

    @Test(dataProvider = "testModelData")
    public void testViterbiOnEmptyData(final TestHMModel model) {
        final List<TestHMModel.State> bestPath = ViterbiAlgorithm.apply(new ArrayList<>(), new ArrayList<>(), model);
        Assert.assertNotNull(bestPath);
        Assert.assertEquals(bestPath, new ArrayList<>());
    }

    @Test(dataProvider = "testModelData")
    public void testForwardBackwardOnEmptyData(final TestHMModel model) {
        final ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> fbResult = ForwardBackwardAlgorithm.apply(new ArrayList<>(), new ArrayList<>(), model);
        Assert.assertNotNull(fbResult);
        Assert.assertEquals(fbResult.data(), new ArrayList<>());
        Assert.assertEquals(fbResult.positions(), new ArrayList<>());
        Assert.assertEquals(fbResult.logDataLikelihood(), 0.0);
        Assert.assertEquals(fbResult.logProbability(new ArrayList<>()), 0.0);
    }

    // Checks that if all transmissions, priors and emission have the same prob,
    // all possible paths have the same
    // prob as well.
    @Test
    public void testAllPathEqualProbabilityProperty() {
       final int numStates = 3;
       final int length = 5; // dont make this too big as the number of paths to test is numStates ^ length.
       final HiddenMarkovModel<Integer, Integer, Integer> model = new UninformativeTestHMModel(numStates);
       final List<Integer> positions = IntStream.range(0, length).boxed().collect(Collectors.toList());
       final List<Integer> data = Collections.nCopies(length, -1); // data is irrelevant for this model.
       final ForwardBackwardAlgorithm.Result<Integer, Integer, Integer> fbResult =
               ForwardBackwardAlgorithm.apply(data, positions, model);
       final List<Integer> states = model.hiddenStates();
       Assert.assertEquals(states.size(), numStates);
       final List<Integer> nextPath = new ArrayList<>(Collections.nCopies(length, states.get(0)));
       final int numPaths = (int) Math.pow(numStates, length);
       final double expectedPathProb = - Math.log(numPaths);
       final double epsilon = Math.abs(expectedPathProb) * 0.0001;
       for (int i = 0; i < numPaths; i++) {
           increasePath(nextPath, numStates);
           final double logProb = fbResult.logProbability(nextPath);
           Assert.assertEquals(logProb, expectedPathProb, epsilon);
       }
       for (final Integer position : positions) {
           Assert.assertEquals(fbResult.logDataLikelihood(position), 0.0, 0.000001);
       }
    }

    // Checks that when a state is "heavier" than the rest
    // (priors and transitions to it are greater than any other state) and emission
    // priors are uniform, that the path that stays in that state the whole chain
    // is the most likely one.
    @Test
    public void testHeavyStateProperty() {
        final int numStates = 3;
        final int heavyState = numStates >> 1; // the (low) median is the heavy state.
        final double heavyStateWeight = (1.0 / numStates) + 0.0001; // just slightly heavier.
        final int length = 5; // dont make this too big as the number of paths to test is numStates ^ length.
        final HiddenMarkovModel<Integer, Integer, Integer> model = new HeavyStateTestHMModel(numStates, heavyState, heavyStateWeight);
        final List<Integer> positions = IntStream.range(0, length).boxed().collect(Collectors.toList());
        final List<Integer> data = Collections.nCopies(length, -1); // data is irrelevant for this model.
        final List<Integer> states = model.hiddenStates();
        Assert.assertEquals(states.size(), numStates);

        // Checking the Viterbi algorithm behavior:
        final List<Integer> bestPath = ViterbiAlgorithm.apply(data, positions, model);
        final List<Integer> expectedBestPath = Collections.nCopies(length, heavyState);
        Assert.assertEquals(bestPath, expectedBestPath);

        // Checking th FW algorithm behavior:
        final ForwardBackwardAlgorithm.Result<Integer, Integer, Integer> fbResult =
                ForwardBackwardAlgorithm.apply(data, positions, model);

        final int numPaths = (int) Math.pow(numStates, length);
        final double uniformPathProb = - Math.log(numPaths);

        final double heavyPathProb = fbResult.logProbability(expectedBestPath);
        @SuppressWarnings("all")
        final double lightPathProb = fbResult.logProbability(Collections.nCopies(length, heavyState == 0 ? 1 : 0));

        Assert.assertTrue(heavyPathProb > uniformPathProb);
        Assert.assertTrue(lightPathProb < uniformPathProb);
    }

    // Checks that the FW algorithm produce posterior probabilities for the
    // hidden states at each position (i) in the chain accordingly with what its
    // expected given a random prior (P) and transition matrix (T): P x T^i .
    @Test
    public void testMultiplicativeTransitionsProperty() {
        final int numStates = 3;
        final int length = 5;
        final Random rdn = new Random(313123);
        final RealVector priors = randomPriors(rdn, numStates);
        final RealMatrix transitions = randomTransitions(rdn, numStates);
        final FlatRealTestHMModel model = new FlatRealTestHMModel(priors, transitions);
        final List<Integer> positions = IntStream.range(0, length).boxed().collect(Collectors.toList());
        final List<Integer> data = Collections.nCopies(length, -1); // data is irrelevant for this model.

        final ForwardBackwardAlgorithm.Result<Integer, Integer, Integer> fbResult =
                ForwardBackwardAlgorithm.apply(data, positions, model);

        for (int p = 0; p < positions.size(); p++) {
            final RealVector expectedProbs = transitions.power(p).preMultiply(priors);
            for (int i = 0; i < numStates; i++) {
                Assert.assertEquals(fbResult.logProbability(positions.get(p), (Integer) i),
                        fbResult.logProbability(p, (Integer) i));
                Assert.assertEquals(fbResult.logProbability(positions.get(p), (Integer) i),
                        Math.log(expectedProbs.getEntry(i)), 0.0001, "p = " + p + " i = " + i);
            }
        }

        transitions.preMultiply(priors);

    }

    private RealVector randomPriors(final Random rdn, final int numStates) {
        final int maxPhred = 30;
        final double[] values = IntStream.range(0, numStates)
                .map(i -> rdn.nextInt(maxPhred + 1) - 1).mapToDouble(i -> Math.pow(10, -.1 *i)).toArray();
        final double sum = DoubleStream.of(values).sum();

        final double[] normalizedValues = DoubleStream.of(values).map(v -> v / sum).toArray();
        return new ArrayRealVector(normalizedValues);
    }

    private RealMatrix randomTransitions(final Random rdn, final int numStates) {
        final int maxPhred = 30;
        final RealMatrix result = new Array2DRowRealMatrix(numStates, numStates);
        final double[] rowTotals = new double[numStates];
        result.walkInRowOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                final double phred = rdn.nextInt(maxPhred + 1) - 1;
                final double prob =  Math.pow(10, phred * -.1);
                rowTotals[row] += prob;
                return prob;
            }
        });
        result.walkInColumnOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return value / rowTotals[row];
            }
        });
        return result;
    }

    private void increasePath(final List<Integer> path, final int numStates) {
        for (int i = 0; i < path.size(); i++) {
            final int state = path.get(i);
            final int newState = (state + 1) % numStates;
            path.set(i, newState);
            if (newState != 0) {
                break;
            }
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testForwardBackwardWithMismatchedDataPositionsLength() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Y)),
                new ArrayList<>(Collections.singletonList(1)), model);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testViterbiWithMismatchedDataPositionsLength() {
        final TestHMModel model = TEST_MODELS.get(0);
        ViterbiAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Y)),
                new ArrayList<>(Collections.singletonList(1)), model);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testForwardBackwardWithNullDataPoints() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, null, TestHMModel.Datum.Y)),
                new ArrayList<>(Arrays.asList(1, 2, 3)), model);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithNullState() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(2, (TestHMModel.State)null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithNullStateUsingPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(new Integer(2), (TestHMModel.State)null);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogDataLikelihoodWithWrongTarget() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        result.logDataLikelihood(14);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogDataLikelihoodWithWrongTargetUsingPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        result.logDataLikelihood(new Integer(14));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBackProbabilityWithNullState() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logBackwardProbability(2, null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBackProbabilityWithNullStateUsingPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logBackwardProbability(new Integer(2), null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogForwardProbabilityWithNullState() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logForwardProbability(2, null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogForwardProbabilityWithNullStateUsingPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logForwardProbability(new Integer(2), null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithNullPosition() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(null, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithNegativePositionIndex() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(-13, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithMadeUpPosition() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(14, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithMadeUpPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(new Integer(14), TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithASequenceThatIsTwoLong() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(2, Arrays.asList(TestHMModel.State.A, TestHMModel.State.A, TestHMModel.State.B));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithASequenceThatIsTwoLongUsingPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(new Integer(2), Arrays.asList(TestHMModel.State.A, TestHMModel.State.A, TestHMModel.State.B));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBackwardProbabilityWithNullPosition() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logBackwardProbability(null, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBackwardProbabilityWithNegativePositionIndex() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logBackwardProbability(-13, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogForwardProbabilityWithNullPosition() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logForwardProbability(null, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogForwardProbabilityWithNegativePositionIndex() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logForwardProbability(-13, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithNonExistentPosition() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(4, TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithNonExistentPositionObject() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(new Integer(4), TestHMModel.State.A);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithANull() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogProbabilityWithSequenceWrongLength() {
        final TestHMModel model = TEST_MODELS.get(0);
        ForwardBackwardAlgorithm.Result<TestHMModel.Datum, Integer, TestHMModel.State> result = null;
        try {
            result = ForwardBackwardAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, TestHMModel.Datum.Z, TestHMModel.Datum.Y)),
                    new ArrayList<>(Arrays.asList(1, 2, 3)), model);
        } catch (final Exception ex) {
            Assert.fail("this part of the test must not fail");

        }
        Assert.assertNotNull(result);
        result.logProbability(Arrays.asList(TestHMModel.State.A, TestHMModel.State.A));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testViterbiWithNullDataPoints() {
        final TestHMModel model = TEST_MODELS.get(0);
        ViterbiAlgorithm.apply(new ArrayList<>(Arrays.asList(TestHMModel.Datum.X, null, TestHMModel.Datum.Y)),
                new ArrayList<>(Arrays.asList(1, 2, 3)), model);
    }

    @DataProvider(name = "testViterbiData")
    public Object[][] testViterbiData() throws IOException {
        if (TEST_EXPECTED_RESULTS == null) {
            setUp();
        }
        return TEST_EXPECTED_RESULTS.stream()
                .map(er -> new Object[]{er.model, er.data, er.bestPath})
                .toArray(Object[][]::new);
    }

    @DataProvider(name = "testForwardData")
    public Object[][] testForwardData() throws IOException {
        if (TEST_EXPECTED_RESULTS == null) {
            setUp();
        }
        return TEST_EXPECTED_RESULTS.stream()
                .map(er -> new Object[]{er.model, er.data, er.logForwardProbs})
                .toArray(Object[][]::new);
    }

    @DataProvider(name = "testFBResultData")
    public Object[][] testFBResultData() throws IOException {
        if (TEST_EXPECTED_RESULTS == null) {
            setUp();
        }
        return TEST_EXPECTED_RESULTS.stream()
                .map(er -> new Object[]{er.model, er.data})
                .toArray(Object[][]::new);
    }

    @DataProvider(name = "testModelData")
    public Object[][] testModelData() throws IOException {
        return TEST_MODELS.stream()
                .map(model -> new Object[]{ model })
                .toArray(Object[][]::new);
    }

    @DataProvider(name = "testBackwardData")
    public Object[][] testBackwardData() throws IOException {
        if (TEST_EXPECTED_RESULTS == null) {
            setUp();
        }
        return TEST_EXPECTED_RESULTS.stream()
                .map(er -> new Object[]{er.model, er.data, er.logBackwardProbs})
                .toArray(Object[][]::new);
    }

    @DataProvider(name = "testPosteriorData")
    public Object[][] testPosteriorData() throws IOException {
        if (TEST_EXPECTED_RESULTS == null) {
            setUp();
        }
        return TEST_EXPECTED_RESULTS.stream()
                .map(er -> new Object[]{er.model, er.data, er.logProbabilities})
                .toArray(Object[][]::new);
    }

    /**
     * Expected result.
     *
     * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
     */
    public static class ExpectedResult {

        public final TestHMModel model;
        public final List<TestHMModel.Datum> data;
        public final List<TestHMModel.State> bestPath;
        public final List<double[]> logBackwardProbs;
        public final List<double[]> logForwardProbs;
        public final List<double[]> logProbabilities;

        public ExpectedResult(final TestHMModel model, final List<TestHMModel.Datum> data, final List<TestHMModel.State> bestPath,
                              final List<double[]> logForwardProbs, final List<double[]> logBackwardProbs,
                              final List<double[]> logProbabilities) {
            this.model = model;
            this.data = data;
            this.bestPath = bestPath;
            this.logForwardProbs = logForwardProbs;
            this.logBackwardProbs = logBackwardProbs;
            this.logProbabilities = logProbabilities;
        }
    }

    private static class RResultReader extends TableReader<RResultRecord> {
        public RResultReader(File file) throws IOException {
            super(file);
        }

        @Override
        protected RResultRecord createRecord(final DataLine dataLine) {
            return new RResultRecord(dataLine);
        }
    }

    private static class RResultRecord {
        private final int modelIndex;
        private final int sequenceIndex;
        private final int dataIndex;
        private final TestHMModel.State bestPathState;
        private final double[] logForwardProbs;
        private final double[] logBackwardProbs;
        private final double[] logProbability;

        RResultRecord(final DataLine dataLine) {
            modelIndex = dataLine.getInt("MODEL");
            sequenceIndex = dataLine.getInt("SEQUENCE");
            dataIndex = dataLine.getInt("POSITION");
            bestPathState = TestHMModel.State.valueOf(dataLine.get("BEST_PATH"));
            logForwardProbs = extractLogForwardProbs(dataLine);
            logBackwardProbs = extractLogBackwardProbs(dataLine);
            logProbability = extractLogPosteriorProbability(dataLine);
        }

        private double[] extractLogForwardProbs(final DataLine dataLine) {
            return extractLogProbs(dataLine, "FW_");
        }

        private double[] extractLogBackwardProbs(final DataLine dataLine) {
            return extractLogProbs(dataLine, "BW_");
        }

        private double[] extractLogPosteriorProbability(final DataLine dataLine) {
            return extractLogProbs(dataLine, "PP_");
        }

        private double[] extractLogProbs(final DataLine dataLine, final String prefix) {
            final double[] result = new double[TestHMModel.State.values().length];
            for (final TestHMModel.State state : TestHMModel.State.values()) {
                result[state.ordinal()] = dataLine.getDouble(prefix + state.name());
            }
            return result;
        }
    }
}