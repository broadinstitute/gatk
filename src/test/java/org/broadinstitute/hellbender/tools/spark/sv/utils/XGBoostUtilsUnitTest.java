package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class XGBoostUtilsUnitTest  extends GATKBaseTest {
    private static final String SV_UTILS_TEST_DIR = toolsTestDir + "spark/sv/utils/";
    private static final String TEST_MATRIX_DATA_FILE = SV_UTILS_TEST_DIR + "agaricus-integers.csv.gz";
    private static final MachineLearningUtils.TruthSet TEST_SET = MachineLearningUtils.loadCsvFile(TEST_MATRIX_DATA_FILE);
    private static final XGBoostUtils.GATKXGBooster classifier = new XGBoostUtils.GATKXGBooster();
    private static final int NUM_TRAINING_ROUNDS = XGBoostUtils.DEFAULT_NUM_TRAINING_ROUNDS;
    private static final int NUM_EVAL_METRIC_TRAINING_ROUNDS = 100; // for tests where accurate results are not needed
    private static final int EARLY_STOPPING_ROUNDS = XGBoostUtils.DEFAULT_EARLY_STOPPING_ROUNDS;
    private static final int NUM_CROSSVALIDATION_FOLDS = MachineLearningUtils.DEFAULT_NUM_CROSSVALIDATION_FOLDS;
    private static final int NUM_TUNING_ROUNDS = XGBoostUtils.DEFAULT_NUM_TUNING_ROUNDS; // keep tests quick
    private static final double TUNING_FRACTION = 1.0 / (1.0 + NUM_TUNING_ROUNDS);
    private static final Map<String, Object> CLASSIFIER_PARAMS = XGBoostUtils.DEFAULT_CLASSIFIER_PARAMETERS;
    private static final String[] TEST_EVAL_METRICS = {
            "rmse", "mae", "logloss", "error", "error@0.75", "auc", "ndcg", "map", "map@7000", "poisson-nloglik"
    };
    private static final Random random = new Random(XGBoostUtils.DEFAULT_SEED);

    private static final double MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY = 0.99;

    private static void assertLabelsEqual(final int[] actuals, final int[] expecteds, final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], "at index=" + index + ": " + message);
        }
    }

    @Test(groups = "sv")
    protected void testBasicTrain() throws IOException {
        // check that a classifier can be trained from data
        final Map<String, Object> classifierParams = new HashMap<>(CLASSIFIER_PARAMS);
        classifierParams.put(XGBoostUtils.NUM_TRAINING_ROUNDS_KEY, NUM_EVAL_METRIC_TRAINING_ROUNDS);
        classifier.train(classifierParams, TEST_SET);

        // check that the classifier can make predictions on the same data;
        final int[] predictedLabels = classifier.predictClassLabels(TEST_SET.features);

        // check that the predictions match the data (this data is easy to predict, they should)
        assertLabelsEqual(predictedLabels, TEST_SET.classLabels,
                "Predicted labels not identical to actual labels");

        // predict probabilities of whole matrix
        final double[][] probabilities = classifier.predictProbability(TEST_SET.features);
        // check that you get indentical results when predicting row-by-row
        for(int row = 0; row < TEST_SET.getNumRows(); ++row) {
            final double[] rowProbability = classifier.predictProbability(TEST_SET.features.getRow(row));
            assertArrayEquals(rowProbability, probabilities[row], 0,
                    "Row " + row + ": different probabilities predicted for matrix and row-by-row");
        }

        // save classifier to temporary file
        File tempFile = File.createTempFile("gatk-xgboost-classifier", "kryo");
        classifier.save(tempFile.getAbsolutePath());
        // load classifier from temporary file

        final MachineLearningUtils.GATKClassifier loadedClassifier = MachineLearningUtils.GATKClassifier.load(
                tempFile.getAbsolutePath()
        );
        // check that you get identical results to first probability predictions
        final double[][] loadedProbabilities = loadedClassifier.predictProbability(TEST_SET.features);
        assertMatrixEquals(loadedProbabilities, probabilities, 0.0,
                "Probabilities predicted by loaded classifier not equal to original");
    }

    @Test(groups = "sv")
    protected void testGetTrainingTrace() {
        final int[] stratify = TEST_SET.classLabels;

        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_SET.getNumRows(), random, stratify
        );
        final MachineLearningUtils.TruthSet trainSet = TEST_SET.sliceRows(hyperSplit.trainRows);
        final MachineLearningUtils.TruthSet validateSet = TEST_SET.sliceRows(hyperSplit.testRows);

        final double[] trainingTrace = classifier.trainAndReturnQualityTrace(CLASSIFIER_PARAMS, trainSet, validateSet,
                NUM_TRAINING_ROUNDS, EARLY_STOPPING_ROUNDS);
        Assert.assertEquals(
                trainingTrace.length, NUM_TRAINING_ROUNDS,
                "Training trace did not have requested number of rounds (" + trainingTrace.length
                        + " instead of " + NUM_TRAINING_ROUNDS
        );
    }

    @Test(groups = "sv")
    protected void testGetMaximizeEvalMetric() {
        final int[] stratify = TEST_SET.classLabels;

        // Just get a small subset of data, and measure what the training algorithm is trying to do on that data.
        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_SET.getNumRows(), random, stratify
        );
        final MachineLearningUtils.TruthSet trainSet = TEST_SET.sliceRows(hyperSplit.trainRows);
        final MachineLearningUtils.TruthSet validateSet = trainSet;

        final XGBoostUtils.GATKXGBooster classifier = new XGBoostUtils.GATKXGBooster();
        final Map<String, Object> classifierParams = new HashMap<>(XGBoostUtils.DEFAULT_CLASSIFIER_PARAMETERS);
        for(final String evalMatric : TEST_EVAL_METRICS) {
            classifierParams.put(XGBoostUtils.EVAL_METRIC_KEY, evalMatric);
            final boolean maximizeEvalMetric = classifier.getMaximizeEvalMetric(classifierParams);
            final double[] trainingTrace = classifier.trainAndReturnQualityTrace(
                    classifierParams, trainSet, validateSet, NUM_EVAL_METRIC_TRAINING_ROUNDS, Integer.MAX_VALUE);
            final double startVal = trainingTrace[0];
            final DoubleSummaryStatistics traceStats = Arrays.stream(trainingTrace).summaryStatistics();
            final double maxVal = traceStats.getMax();
            final double minVal = traceStats.getMin();
            final boolean traceIncreasedMoreThanDecreased = maxVal - startVal > startVal - minVal;
            Assert.assertEquals(traceIncreasedMoreThanDecreased, maximizeEvalMetric,
                    "For " + XGBoostUtils.EVAL_METRIC_KEY + "=" + evalMatric + ": maximizeEvalMetric="
                            + maximizeEvalMetric + ", but traceIncreasedMoreThanDecreased=" + traceIncreasedMoreThanDecreased);
        }
    }

    @Test(groups = "sv")
    protected void testCrossvalidatedTuneAndTrain() {
        final int minCountsPerStratifyValue = (int)Math.round(NUM_CROSSVALIDATION_FOLDS / TUNING_FRACTION);
        final int[] stratify = TEST_SET.getStratifyArray(MachineLearningUtils.DEFAULT_NUM_STRATIFY_BINS,
                                                         minCountsPerStratifyValue);

        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                TUNING_FRACTION, TEST_SET.getNumRows(), random, stratify
        );
        final MachineLearningUtils.TruthSet tuneSet = TEST_SET.sliceRows(hyperSplit.trainRows);
        final int[] tuneStratify = MachineLearningUtils.slice(stratify, hyperSplit.trainRows);
        final MachineLearningUtils.TruthSet validateSet =TEST_SET.sliceRows(hyperSplit.testRows);
        final int[] validateStratify = MachineLearningUtils.slice(stratify, hyperSplit.testRows);

        // set the number of threads
        classifier.chooseNumThreads(CLASSIFIER_PARAMS, XGBoostUtils.NUM_THREADS_KEY, tuneSet);

        final Map<String, Object> bestClassifierParameters = classifier.tuneClassifierParameters(
                CLASSIFIER_PARAMS, XGBoostUtils.DEFAULT_TUNING_PARAMETERS,
                MachineLearningUtils.ClassifierTuningStrategy.RANDOM,
                tuneSet, random, tuneStratify, NUM_CROSSVALIDATION_FOLDS, NUM_TRAINING_ROUNDS,
                EARLY_STOPPING_ROUNDS, NUM_TUNING_ROUNDS
        );

        final int[] predictedTestLabels = classifier.crossvalidatePredict(
                validateSet, bestClassifierParameters, random, validateStratify, NUM_CROSSVALIDATION_FOLDS
        );

        final double accuracy = MachineLearningUtils.getPredictionAccuracy(
                predictedTestLabels, validateSet.classLabels
        );

        Assert.assertTrue(accuracy >= MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY,
                "Crossvalidated prediction accuracy (" + accuracy + ") less than passing ("
                        + MINIMUM_ALLOWED_CROSSVALIDATED_ACCURACY + ")");
    }

    private static void assertArrayEquals(final double[] actuals, final double[] expecteds, final double tol,
                                          final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Lengths not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            Assert.assertEquals(actuals[index], expecteds[index], tol, "at index=" + index + ": " + message);
        }
    }

    private static void assertMatrixEquals(final double[][] actuals, final double[][] expecteds, final double tol,
                                           final String message) {
        Assert.assertEquals(actuals.length, expecteds.length, "Number of rows not equal: " + message);
        for(int index = 0; index < expecteds.length; ++index) {
            assertArrayEquals(actuals[index], expecteds[index], tol, "at row=" + index + ": " + message);
        }
    }
}
