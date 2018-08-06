package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class MachineLearningUtilsUnitTest extends GATKBaseTest {
    private static final Random random = new Random();
    // arbitrary number to test some of the array-building / manipulating routines that should not expected to be
    // brittle to differences between positive numbers. Choosing a large one to verify that it scales to reasonable getNumRows
    // data sets.
    private static final int NUM_ELEMENTS_TEST = 1001001;
    private static final int NUM_COLUMNS_TEST = 4;
    private static final int NUM_CLASSES_TEST = 3;
    private static final double TRAINING_FRACTION = 0.3;
    private static final int NUM_CROSSVALIDATION_FOLDS = 5;
    // choose an arbitrary large value in between 1 and NUM_ELEMENTS_TEST
    private static final int NUM_STRATIFY_CLASSES = (int)Math.round(Math.sqrt(NUM_ELEMENTS_TEST));
    private static final int MIN_COUNTS_PER_STRATIFY_BIN = 1000;
    // choose arbitrary number of ClassifierParamRange samples to test
    private static final int NUM_PARAM_RANGE_SAMPLES = 60;
    // choose arbitrary ranges to test over
    private static final double LINEAR_LOW = -1.0;
    private static final double LINEAR_HIGH = 1.0;
    private static final double LOG_LOW = 0.1;
    private static final double LOG_HIGH = 10.0;
    private static final int LINEAR_INT_LOW = -100;
    private static final int LINEAR_INT_HIGH = 100;
    private static final int LOG_INT_LOW = 1;
    private static final int LOG_INT_HIGH = 100;
    private static final double[] CHECK_RANGE_BALANCE_PROPORTION = new double[] {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    @Test(groups = "sv")
    protected void testGetRange() {
        try {
            MachineLearningUtils.getRange(-1);
            Assert.fail("getRange(negativeInteger) should throw IllegalArgumentException");
        } catch (IllegalArgumentException ignored) {
        }
        Integer negativeInteger = -1;
        try {
            MachineLearningUtils.getRange(negativeInteger);
            Assert.fail("getRange(negativeInteger) should throw IllegalArgumentException");
        } catch (IllegalArgumentException ignored) {
        }
        assertEmpty(MachineLearningUtils.getRange(0), "getRange(0) should be empty");
        assertEmpty(MachineLearningUtils.getRange((Integer) 0), "getRange((Integer)0) should be empty");

        // if these work for any given positive integer they'll work for the rest. Just test one.
        final int[] intRange = MachineLearningUtils.getRange(NUM_ELEMENTS_TEST);
        final Integer[] integerRange = MachineLearningUtils.getRange((Integer) NUM_ELEMENTS_TEST);
        Assert.assertEquals(intRange.length, NUM_ELEMENTS_TEST, "getRange(int) returned wrong number of elements");
        Assert.assertEquals(integerRange.length, NUM_ELEMENTS_TEST, "getRange(Integer) returned wrong number of elements");
        for (int j = 0; j < NUM_ELEMENTS_TEST; ++j) {
            Assert.assertEquals(intRange[j], j, "getRange(int) has wrong value");
            Assert.assertEquals(integerRange[j], (Integer) j, "getRange(Integer) has wrong value");
        }
    }

    @Test(groups = "sv")
    protected void testArrayOrderManipulation() {
        try {
            MachineLearningUtils.getRandomPermutation(random, -1);
            Assert.fail("getRandomPermutation(random, negativeInteger) should throw IllegalArgumentException");
        } catch (IllegalArgumentException ignored) {
        }
        assertEmpty(MachineLearningUtils.getRandomPermutation(random, 0), "getRandomPermutation(random, 0) should be empty");

        // if these work for any given positive integer they'll work for the rest. Just test one.
        final int[] permutation = MachineLearningUtils.getRandomPermutation(random, NUM_ELEMENTS_TEST);
        Assert.assertEquals(permutation.length, NUM_ELEMENTS_TEST, "getRandomPermutation() returned wrong number of elements");
        // assert the permutation is non-trivial. **Technically** it's possible this will fail, but it's so unlikely it
        // should never happen during the epoch where the universe is capable of hosting life...
        final int[] range = MachineLearningUtils.getRange(NUM_ELEMENTS_TEST);
        assertArraysNotEqual(permutation, range,"getPermutation returned trivial permutaiton.");

        // use argsort to find indices that invert the permutation.
        final int[] sortedInds = MachineLearningUtils.argsort(permutation);
        // it's easy to use sliceAssign to get the correct inverse permutation
        final int[] inverseInds = new int[NUM_ELEMENTS_TEST];
        MachineLearningUtils.sliceAssign(inverseInds, permutation, range);
        assertArrayEquals(sortedInds, inverseInds, "sortedInds don't equal inverseInds");
        // using slice should put the permutation back into order
        Assert.assertEquals(MachineLearningUtils.slice(permutation, sortedInds), range,
                "argsort did not correctly unscramble the permutation.");
        // check corner-cases
        assertEmpty(MachineLearningUtils.slice(permutation, new int[0]), "empty slice should return empty result");
        assertEmpty(MachineLearningUtils.argsort(new int[0]), "argsort of empty array should return empty result");
        try {
            MachineLearningUtils.slice(permutation, new int[] {NUM_ELEMENTS_TEST});
            Assert.fail("slice should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }
        try {
            MachineLearningUtils.slice(permutation, new int[] {-1});
            Assert.fail("slice should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }

        // get two differently-sized arbitrary arrays of indices that are a subset of valid indices into permutation.
        final int[] sliceArr1 = MachineLearningUtils.slice(
                permutation, MachineLearningUtils.getRange(NUM_ELEMENTS_TEST / 2)
        );
        final int[] sliceArr2 = MachineLearningUtils.slice(
                permutation, MachineLearningUtils.getRange(NUM_ELEMENTS_TEST / 4)
        );
        try {
            MachineLearningUtils.sliceAssign(permutation, sliceArr1, sliceArr2);
            Assert.fail("sliceAssign should throw IllegalArgumentException when slice indices don't have same getNumRows as slice values");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.sliceAssign(permutation, sliceArr2, sliceArr1);
            Assert.fail("sliceAssign should throw IllegalArgumentException when slice indices don't have same getNumRows as slice values");
        } catch (IllegalArgumentException ignored) {
        }
        // sliceAssign should throw IndexOutOfBoundsException when asked for a bad index
        try {
            MachineLearningUtils.sliceAssign(permutation, new int[] {NUM_ELEMENTS_TEST}, new int[] {42});
            Assert.fail("sliceAssign should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }
        try {
            MachineLearningUtils.sliceAssign(permutation, new int[] {-1}, new int[] {42});
            Assert.fail("sliceAssign should throw ArrayIndexOutOfBoundsException when passed a bad index");
        } catch (ArrayIndexOutOfBoundsException ignored) {
        }
    }

    @Test(groups = "sv")
    /**
     * Test forming stratify array from TruthSet. Don't check ability to balance data set, that is checked in
     * testTrainTestSplit
     */
    protected void testGetStratifyArray() {
        final int numRows = NUM_ELEMENTS_TEST;
        // pick something arbitrary that makes the binning do something
        final int numDataValues = MachineLearningUtils.DEFAULT_NUM_STRATIFY_BINS * 2;

        final RealMatrix features = new Array2DRowRealMatrix(
                IntStream.range(0, numRows).mapToObj(
                        r -> IntStream.range(0, NUM_COLUMNS_TEST).mapToDouble(
                                c -> (double) random.nextInt(numDataValues)
                        ).toArray()
                ).toArray(double[][]::new)
        );
        final int[] classLabels = IntStream.range(0, numRows).map(r->random.nextInt(NUM_CLASSES_TEST)).toArray();

        final MachineLearningUtils.TruthSet truthSet = new MachineLearningUtils.TruthSet(features, classLabels);

        final int[] stratify = truthSet.getStratifyArray(MachineLearningUtils.DEFAULT_NUM_STRATIFY_BINS,
                                                         MIN_COUNTS_PER_STRATIFY_BIN);
        // check length is consistent with num rows
        Assert.assertEquals(stratify.length, truthSet.getNumRows(), "stratify has incorrect length");
        final Map<Integer, List<Integer>> stratifyInds = IntStream.range(0, stratify.length).boxed().collect(
                Collectors.groupingBy(index -> (Integer)stratify[index])
        );
        // Check values are contiguous [0, 1, ..., maxStratifyValue]
        final int[] stratifyValues = stratifyInds.keySet().stream().sorted().mapToInt(Integer::intValue).toArray();
        assertArrayEquals(stratifyValues, MachineLearningUtils.getRange(stratifyValues.length),
                "stratify values are not successive integers starting at 0");
        // Check no bins have too few counts
        for(final List<Integer> inds : stratifyInds.values()) {
            final int stratifyCount = inds.size();
            Assert.assertTrue(stratifyCount >= MIN_COUNTS_PER_STRATIFY_BIN,
                    "Found stratify value with " + stratifyCount + " counts, which is below minimum allowed ("
                            + MIN_COUNTS_PER_STRATIFY_BIN + ")");
        }
        // Check that there are more stratify values than NUM_CLASSES_TEST or numDataValues, i.e. data was melded
        Assert.assertTrue(stratifyInds.size() > Math.max(NUM_CLASSES_TEST, numDataValues),
                "stratify did not meld data from multiple columns. Produced " + stratifyInds.size()
                        + " stratify values from " + NUM_CLASSES_TEST + " classes and " + numDataValues + " bins");
        // Check that stratify values respect class labels (i.e. for each stratify value, all its elements have exactly
        // one class label
        for(final List<Integer> inds: stratifyInds.values()) {
            final long numClassLabels = Arrays.stream(
                    MachineLearningUtils.slice(
                        truthSet.classLabels,
                        inds.stream().mapToInt(Integer::intValue).toArray()
                )
            ).distinct().count();
            Assert.assertEquals(numClassLabels, 1, "stratify value encoded rows with more than one class label");
        }

    }

    @Test(groups = "sv")
    /**
     * Test expected behavior of TrainTestSplit factory functions.
     * Basic functionality is tested in actual use of tuning / cross-validating. Here check that the details match
     * expectations (data set is divided evenly according to specification).
     */
    protected void testTrainTestSplit() {
        // Create stratify array, random array of class membership indexes used to stratify splits. Make it class
        // populations uneven to further stress-test algorithms.
        final int[] stratify = Arrays.stream(
                new MachineLearningUtils.ClassifierIntegerLogParamRange(1, NUM_STRATIFY_CLASSES)
                .getRandomSamples(random, NUM_ELEMENTS_TEST)
        ).mapToInt(i->i).toArray();
        Assert.assertEquals(stratify.length, NUM_ELEMENTS_TEST, "stratify was generated improperly");

        // Test stratified and un-stratified train-test split
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        TRAINING_FRACTION, NUM_ELEMENTS_TEST, random, stratify
                ),
                NUM_ELEMENTS_TEST, TRAINING_FRACTION, stratify
        );

        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        TRAINING_FRACTION, NUM_ELEMENTS_TEST, random, null
                ),
                NUM_ELEMENTS_TEST, TRAINING_FRACTION, null
        );


        // Test stratified and un-stratified cross-validation splits
        final double cvTrainingFraction = 1.0 - 1.0 / (double)NUM_CROSSVALIDATION_FOLDS;

        final Iterator<MachineLearningUtils.TrainTestSplit> stratifiedCvSplits
                = MachineLearningUtils.TrainTestSplit.getCrossvalidationSplits(
                        NUM_CROSSVALIDATION_FOLDS, NUM_ELEMENTS_TEST, random, stratify
            );
        while(stratifiedCvSplits.hasNext()) {
            final MachineLearningUtils.TrainTestSplit stratifiedCvSplit = stratifiedCvSplits.next();
            assertGoodSplit(stratifiedCvSplit, NUM_ELEMENTS_TEST, cvTrainingFraction, stratify);
        }

        final Iterator<MachineLearningUtils.TrainTestSplit> flatCvSplits
                = MachineLearningUtils.TrainTestSplit.getCrossvalidationSplits(
                NUM_CROSSVALIDATION_FOLDS, NUM_ELEMENTS_TEST, random, null
        );
        while(flatCvSplits.hasNext()) {
            final MachineLearningUtils.TrainTestSplit flatCvSplit = flatCvSplits.next();
            assertGoodSplit(flatCvSplit, NUM_ELEMENTS_TEST, cvTrainingFraction, null);
        }

        // check corner-cases
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        0.0, NUM_ELEMENTS_TEST, random, stratify
                ),
                NUM_ELEMENTS_TEST, 0.0, stratify
        );
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        1.0, NUM_ELEMENTS_TEST, random, stratify
                ),
                NUM_ELEMENTS_TEST, 1.0, stratify
        );
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        0.0, NUM_ELEMENTS_TEST, random, null
                ),
                NUM_ELEMENTS_TEST, 0.0, null
        );
        assertGoodSplit(
                MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                        1.0, NUM_ELEMENTS_TEST, random, null
                ),
                NUM_ELEMENTS_TEST, 1.0, null
        );

        // ensure errors are thrown when crazy values are passed
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    -0.01, NUM_ELEMENTS_TEST, random, stratify
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    1.01, NUM_ELEMENTS_TEST, random, stratify
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    -0.01, NUM_ELEMENTS_TEST, random, null
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }
        try {
            MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                    1.01, NUM_ELEMENTS_TEST, random, null
            );
            Assert.fail("getTrainTestSplit should throw IllegalArgmentException when trainingFraction is not in range [0, 1]");
        } catch (IllegalArgumentException ignored) {
        }

        try {
            MachineLearningUtils.TrainTestSplit.getCrossvalidationSplits(
                    1, NUM_ELEMENTS_TEST, random, stratify
            );
            Assert.fail("getCrossvalidationSplits should throw IllegalArgmentException when numCrossvalidationFolds < 2");
        } catch (IllegalArgumentException ignored) {
        }
    }

    @Test(groups = "sv")
    /**
     * Test expected behavior of ParamRange subclasses.
     * Basic functionality is tested in actual use of tuning / cross-validating. Here check that the details match
     * expectations (getRandomSamples returns randomly-ordered samples from range, evenly spaced linearly or logarithmically,
     * with different param ranges scrambled independently, and of the appropriate data type).
     */
    protected void testParamRanges() {
        final Double[] linearSamples = new MachineLearningUtils.ClassifierLinearParamRange(
                LINEAR_LOW, LINEAR_HIGH
        ).getRandomSamples(random, NUM_PARAM_RANGE_SAMPLES);
        final Double[] logSamples = new MachineLearningUtils.ClassifierLogParamRange(
                LOG_LOW, LOG_HIGH
        ).getRandomSamples(random, NUM_PARAM_RANGE_SAMPLES);
        final Integer[] linearIntSamples = new MachineLearningUtils.ClassifierIntegerLinearParamRange(
                LINEAR_INT_LOW, LINEAR_INT_HIGH
        ).getRandomSamples(random, NUM_PARAM_RANGE_SAMPLES);
        final Integer[] logIntSamples = new MachineLearningUtils.ClassifierIntegerLogParamRange(
                LOG_INT_LOW, LOG_INT_HIGH
        ).getRandomSamples(random, NUM_PARAM_RANGE_SAMPLES);

        // check ranges have appropriate min and max
        assertArrayEquals(
                new double[] {Arrays.stream(linearSamples).min(Double::compare).get(), Arrays.stream(linearSamples).max(Double::compare).get()},
                new double[] {LINEAR_LOW, LINEAR_HIGH},
                "linearSamples have incorrect range"
        );
        assertArrayEquals(
                new double[] {Arrays.stream(logSamples).min(Double::compare).get(), Arrays.stream(logSamples).max(Double::compare).get()},
                new double[] {LOG_LOW, LOG_HIGH},
                "logSamples have incorrect range"
        );
        assertArrayEquals(
                new int[] {Arrays.stream(linearIntSamples).min(Integer::compare).get(), Arrays.stream(linearIntSamples).max(Integer::compare).get()},
                new int[] {LINEAR_INT_LOW, LINEAR_INT_HIGH},
                "linearIntSamples have incorrect range"
        );
        assertArrayEquals(
                new int[] {Arrays.stream(logIntSamples).min(Integer::compare).get(), Arrays.stream(logIntSamples).max(Integer::compare).get()},
                new int[] {LOG_INT_LOW, LOG_INT_HIGH},
                "logIntSamples have incorrect range"
        );

        // check ranges are sampled uniformly
        for(final double checkProportion : CHECK_RANGE_BALANCE_PROPORTION) {
            final double linearDivider = LINEAR_LOW + plusEpsilon(LINEAR_HIGH - LINEAR_LOW) * checkProportion;
            assertEvenParamRange(linearSamples, linearDivider, checkProportion);
            final double logDivider = LOG_LOW * Math.pow(plusEpsilon(LOG_HIGH / LOG_LOW), checkProportion);
            assertEvenParamRange(logSamples, logDivider, checkProportion);
            final double linearIntDivider = LINEAR_INT_LOW + plusEpsilon(LINEAR_INT_HIGH - LINEAR_INT_LOW) * checkProportion;
            assertEvenParamRange(linearIntSamples, linearIntDivider, checkProportion);
            final double logIntDivider = LOG_INT_LOW * Math.pow(plusEpsilon(LOG_INT_HIGH / (double)LOG_INT_LOW), checkProportion);
            assertEvenParamRange(logIntSamples, logIntDivider, checkProportion);
        }

        // check ranges do not co-vary, permutations should be distinct
        final int[] linearPermutation = MachineLearningUtils.argsort(linearSamples);
        final int[] logPermutation = MachineLearningUtils.argsort(logSamples);
        final int[] linearIntPermutation = MachineLearningUtils.argsort(linearIntSamples);
        final int[] logIntPermutation = MachineLearningUtils.argsort(logIntSamples);
        assertArraysNotEqual(linearPermutation, logPermutation, "Permuations repeated, parameter space will not be evenly sampled");
        assertArraysNotEqual(linearPermutation, linearIntPermutation, "Permuations repeated, parameter space will not be evenly sampled");
        assertArraysNotEqual(linearPermutation, logIntPermutation, "Permuations repeated, parameter space will not be evenly sampled");
        assertArraysNotEqual(logPermutation, linearIntPermutation, "Permuations repeated, parameter space will not be evenly sampled");
        assertArraysNotEqual(logPermutation, logIntPermutation, "Permuations repeated, parameter space will not be evenly sampled");
        assertArraysNotEqual(linearIntPermutation, logIntPermutation, "Permuations repeated, parameter space will not be evenly sampled");

        // check InvalidArguments
        assertBadParamRangeThrowsException(LINEAR_HIGH, LINEAR_LOW, MachineLearningUtils.ClassifierLinearParamRange.class,
                IllegalArgumentException.class, "ClassifierParamRange with low > high should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(LOG_HIGH, LOG_LOW, MachineLearningUtils.ClassifierLogParamRange.class,
                IllegalArgumentException.class, "ClassifierParamRange with low > high should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(LINEAR_INT_HIGH, LINEAR_INT_LOW, MachineLearningUtils.ClassifierIntegerLinearParamRange.class,
                IllegalArgumentException.class, "ClassifierParamRange with low > high should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(LOG_INT_HIGH, LOG_INT_LOW, MachineLearningUtils.ClassifierIntegerLogParamRange.class,
                IllegalArgumentException.class, "ClassifierParamRange with low > high should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(-LOG_HIGH, LOG_LOW, MachineLearningUtils.ClassifierLogParamRange.class,
                IllegalArgumentException.class, "Log-distributed ClassifierParamRange with low and high of different signs should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(0.0, LOG_LOW, MachineLearningUtils.ClassifierLogParamRange.class,
                IllegalArgumentException.class, "Log-distributed ClassifierParamRange with low and high of different signs should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(-LOG_INT_HIGH, LOG_INT_LOW, MachineLearningUtils.ClassifierIntegerLogParamRange.class,
                IllegalArgumentException.class, "Log-distributed ClassifierParamRange with low and high of different signs should throw IllegalArgumentException");
        assertBadParamRangeThrowsException(0, LOG_INT_HIGH, MachineLearningUtils.ClassifierIntegerLogParamRange.class,
                IllegalArgumentException.class, "Log-distributed ClassifierParamRange with low and high of different signs should throw IllegalArgumentException");
    }

    private static <R extends MachineLearningUtils.ClassifierParamRange<Double>, E extends Exception>
    void assertBadParamRangeThrowsException(final Double low, final Double high, Class<R> clazz, Class<E> exceptionClazz,
                                            final String message) {
        try {
            final R foo = clazz.getConstructor(double.class, double.class).newInstance(low, high);
            Assert.fail(message);
        } catch(Exception exception) {
            if(!exceptionClazz.isInstance(exception.getCause())) {
                // Wrong type, this is an unexpected error
                throw new GATKException(exception.getClass().getName() + ": " + exception.getCause());
            }
        }
    }
    private static <R extends MachineLearningUtils.ClassifierParamRange<Integer>, E extends Exception>
    void assertBadParamRangeThrowsException(final Integer low, final Integer high, Class<R> clazz, Class<E> exceptionClazz,
                                            final String message) {
        try {
            final R foo = clazz.getConstructor(int.class, int.class).newInstance(low, high);
            Assert.fail(message);
        } catch(Exception exception) {
            if(!exceptionClazz.isInstance(exception.getCause())) {
                // Wrong type, this is an unexpected error
                throw new GATKException(exception.getClass().getName() + ": " + exception.getCause());
            }
        }
    }

    private static double plusEpsilon(final double x) {
        return x + 2 * Math.ulp(x);
    }

    private static void assertEvenParamRange(final Double[] randomSamples, final double divider, final double checkProportion) {
        final long numLessDivider = Arrays.stream(randomSamples).filter(s -> s < divider).count();
        final double actualProportion = numLessDivider / (double)randomSamples.length;
        final double tol = 1.0 / randomSamples.length;
        Assert.assertEquals(actualProportion, checkProportion, tol, "Proportion of samples below divider deviated from expected");
    }

    private static void assertEvenParamRange(final Integer[] linearSamples, final double divider, final double checkProportion) {
        // because of quantization, actual proportion less than divider may vary by quite a bit. Bracket divider and ensure
        // that proportions bracket checkProportion
        final long numLessFloorDivider = Arrays.stream(linearSamples).filter(s -> s < (int)Math.floor(divider)).count();
        final long numLessCeilDivider = Arrays.stream(linearSamples).filter(s -> s <= (int)Math.ceil(divider)).count();
        final long numExpectedLessDivider = Math.round(checkProportion * linearSamples.length);
        Assert.assertTrue(numLessFloorDivider <= numExpectedLessDivider && numExpectedLessDivider <= numLessCeilDivider,
                "Proportion of samples below divider deviated from expected");
    }

    private static void assertEmpty(final int[] arr, final String message) {
        if(arr.length != 0) {
            Assert.fail(message);
        }
    }

    private static <T> void assertEmpty(final T[] arr, final String message) {
        if(arr.length != 0) {
            Assert.fail(message);
        }
    }

    private static void assertArraysNotEqual(final int[] actuals, final int[] expecteds, final String message) {
        if(actuals.length != expecteds.length) {
            return;
        }
        for (int index = 0; index < expecteds.length; ++index) {
            if(actuals[index] != expecteds[index]) {
                return;
            }
        }
        Assert.fail(message);
    }

    private static void assertArrayEquals(final int[] actuals, final int[] expecteds, final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for (int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], "at index=" + index + ": " + message);
        }
    }

    private static void assertArrayEquals(final double[] actuals, final double[] expecteds, final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for (int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], "at index=" + index + ": " + message);
        }
    }

    private static void assertArrayEquals(final int[] actuals, final double[] expecteds, final double tol,
                                          final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for (int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], tol, "at index=" + index + ": " + message);
        }
    }

    private static void assertArraysDisjoint(final int[] arr1, final int[] arr2, final String message) {
        final Set<Integer> arr1Values = Arrays.stream(arr1).boxed().collect(Collectors.toSet());
        for(final int val2 : arr2) {
            if(arr1Values.contains(val2)) {
                Assert.fail(message + ": " + val2 + " in both arrays");
            }
        }
    }

    private static void assertGoodSplit(final MachineLearningUtils.TrainTestSplit split, final int numElements,
                                        final double trainingFraction, final int[] stratify) {
        // first check if the split is trivial
        if(trainingFraction == 0.0) {
            assertEmpty(split.trainRows, "trainRows should be empty");
            assertArrayEquals(split.testRows, MachineLearningUtils.getRange(numElements),
                    "testRows should be all elements");
            return;
        } else if(trainingFraction == 1.0) {
            assertEmpty(split.testRows, "testRows should be empty");
            assertArrayEquals(split.trainRows, MachineLearningUtils.getRange(numElements),
                    "trainRows should be all elements");
            return;
        }

        // next check that it *is* a split
        Assert.assertEquals(split.testRows.length + split.trainRows.length, numElements,
                "number of trainRows + number of testRows != number of elements");
        assertArraysDisjoint(split.testRows, split.trainRows, "trainRows and testRows have overlap");
        final IntSummaryStatistics trainStats = Arrays.stream(split.trainRows).summaryStatistics();
        final IntSummaryStatistics testStats = Arrays.stream(split.testRows).summaryStatistics();
        final int minIndex = Math.min(trainStats.getMin(), testStats.getMin());
        final int maxIndex = Math.max(trainStats.getMax(), testStats.getMax());
        Assert.assertEquals(minIndex, 0, "split contains indices less than 0");
        Assert.assertEquals(maxIndex, numElements - 1, "split contains indices past end of data");

        // finally check balance of division
        Assert.assertEquals(split.trainRows.length, trainingFraction * numElements, 2,
                "Number of training rows differs from expected by more than 2");
        if(stratify != null) {
            // check balance of division for each stratify value
            final double[] expectedStratifyCounts = Arrays.stream(
                    getStratifyCounts(stratify)
            ).mapToDouble(c -> c * trainingFraction).toArray();
            final int[] trainStratifyCounts = getStratifyCounts(MachineLearningUtils.slice(stratify, split.trainRows));
            assertArrayEquals(trainStratifyCounts, expectedStratifyCounts, 2,
                    "split wasn't properly balanced for stratify");
        }
    }

    private static int[] getStratifyCounts(final int[] stratify) {
        final Map<Integer, Integer> stratifyCountsMap = new HashMap<>();
        for(final Integer value: stratify) {
            stratifyCountsMap.put(value, stratifyCountsMap.getOrDefault(value, 0) + 1);
        }
        if(stratifyCountsMap.isEmpty()) {
            return new int[0];
        }
        final int maxValue = stratifyCountsMap.keySet().stream().mapToInt(Integer::intValue).max().getAsInt();
        final int[] stratifyCounts = new int[maxValue + 1];
        for(final Map.Entry<Integer, Integer> entry : stratifyCountsMap.entrySet()) {
            stratifyCounts[entry.getKey()] = entry.getValue();
        }
        return stratifyCounts;
    }
}
