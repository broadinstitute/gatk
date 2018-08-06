package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.KryoSerializable;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * This is a base utils class that provides an API for packages that need to use classifiers. Specific classifiers can
 * be provided by extending GATKClassifier. Additional hyperparameter tuning strategies can be provided by extending
 * ClassifierTuner. These abstract base classes implement some common routines that should make adding new classifiers
 * or tuning strategies a little less tedious.
 */
public class MachineLearningUtils {
    public static final String NUM_TRAINING_ROUNDS_KEY = "num_training_rounds";
    public static final int DEFAULT_NUM_CROSSVALIDATION_FOLDS = 5;
    public static final int DEFAULT_NUM_TUNING_ROUNDS = 100;
    public static final int DEFAULT_NUM_STRATIFY_BINS = 5;

    private static final int NUM_CALIBRATION_TRAINING_ROUNDS = 500;
    private static final int NUM_CALIBRATION_CLASS_ROWS = 500;
    public static final int CLASS_LABEL_COLUMN = 0;
    private static final Logger localLogger = LogManager.getLogger(MachineLearningUtils.class);

    /**
     * Simple class to hold truth sets: matrix of features paired with array of class labels
     * Note: if classLabels is passed in as null, it will be converted to an empty array so that checks for presence of
     * classLabels can just check .length
     * */
    public static class TruthSet {
        public final RealMatrix features;
        public final int[] classLabels;
        TruthSet(final RealMatrix features, final int[] classLabels) {
            if(classLabels != null && classLabels.length > 0 && features.getRowDimension() != classLabels.length) {
                throw new IllegalArgumentException(
                        "classLabels must either be null, empty, or have same length as the number of rows in features."
                );
            }
            this.features = features;
            this.classLabels = classLabels == null ? new int[0] : classLabels;
        }

        public int getNumRows() {
            return features.getRowDimension();
        }

        public int getNumColumns() {
            return features.getColumnDimension();
        }

        public TruthSet sliceRows(final int[] rowIndices) {
            return new TruthSet(
                    MachineLearningUtils.sliceRows(features, rowIndices),
                    slice(classLabels, rowIndices)
            );
        }

        /**
         * Extract a stratify array from a TruthSet. Columns are binned and stratify array values are obtained by finding
         * unique rows.
         * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
         * @param numBins Number of bins into which to group values in each column
         * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
         *                                  the number of columns in consideration is reduced and the under-count rows are
         *                                  searched again for unique values.
         * @return stratifyArray, int array with length = number of rows in features
         */
        public int[] getStratifyArray(final int numBins, final int minCountsPerStratifyValue) {
            return getStratifyArray(numBins, minCountsPerStratifyValue, new HashSet<>());
        }

        /**
         * Extract a stratify array from a TruthSet. Columns are binned and stratify array values are obtained by finding
         * unique rows. This version can specify some columns as categorical (unbinnable).
         * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
         * @param numBins Number of bins into which to group values in each column
         * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
         *                                  the number of columns in consideration is reduced and the under-count rows are
         *                                  searched again for unique values.
         * @param categoricalColumns Columns in categoricalColumns are not binned.
         * @return stratifyArray, int array with length = number of rows in features
         */
        public int[] getStratifyArray(final int numBins, final int minCountsPerStratifyValue,
                                      final Set<Integer> categoricalColumns) {
            // class labels are the most important stratify column, so make them first. Also make them categorical
            final RealMatrix stratifyMatrix = concatenateColumns(
                    new Array2DRowRealMatrix(IntStream.of(classLabels).mapToDouble(Double::valueOf).toArray()),
                    features
            );
            final Set<Integer> concatenatedCategoricalColumns = categoricalColumns == null ? new HashSet<>()
                    : categoricalColumns.stream().map(x -> x + 1).collect(Collectors.toSet());
            concatenatedCategoricalColumns.add(0); // column 0 (classLabel) is categorical
            return stratifyMatrixToStratifyArray(stratifyMatrix, numBins, minCountsPerStratifyValue, concatenatedCategoricalColumns);
        }
    }

    /**
     * Load TruthSet from (possibly gzipped) csv file. The file must represent numeric types only.
     * Lines beginning with "#" will be ignored.
     * Values in CLASS_LABEL_COLUMN of file will be interpreted as class labels.
     */
    public static TruthSet loadCsvFile(final String filename) {
        return loadCsvFile(filename, ",", "#", CLASS_LABEL_COLUMN);
    }

    /**
     * Load TruthSet from (possibly gzipped) file which uses arbitrary delimiter to sparate values.
     * @param delimiter specify separator for values in file
     * @param commentCharacter lines beginning with this value are ignored
     * @param classLabelsColumn values in this column of the file will be interpreted as class labels.
     *                          If classLabelsColumn < 0, then classLabels will be empty
     */
    public static TruthSet loadCsvFile(final String filename, final String delimiter,
                                       final String commentCharacter, final int classLabelsColumn) {
        final int numColumns;
        final List<double[]> rowsList = new ArrayList<>();
        final List<Integer> classLabelsList = new ArrayList<>();
        try (final BufferedReader reader = IOUtil.openFileForBufferedReading(new File(filename))) {
            if(!reader.ready()) {
                throw new GATKException("Unable to read matrix from " + filename);
            }
            getNextCsvFeaturesRow(reader, delimiter, commentCharacter, classLabelsColumn, rowsList, classLabelsList, -1);
            numColumns = rowsList.get(0).length;
            while(reader.ready()) {
                getNextCsvFeaturesRow(reader, delimiter, commentCharacter, classLabelsColumn, rowsList, classLabelsList, numColumns);
            }
        } catch(IOException err) {
            throw new GATKException(err.getMessage());
        }
        return new TruthSet(
                new Array2DRowRealMatrix(rowsList.toArray(new double[0][]), false),
                classLabelsList.stream().mapToInt(Integer::intValue).toArray()
        );
    }

    private static void getNextCsvFeaturesRow(final BufferedReader reader, final String delimiter,
                                              final String commentCharacter, final int classLabelsColumn,
                                              final List<double[]> rowsList, final List<Integer> classLabelsList,
                                              final int numColumns) throws IOException {
        String line = reader.readLine();
        while(line.startsWith(commentCharacter) || line.isEmpty()) {
            line = reader.readLine();
        }
        final String[] words = line.split(delimiter, -1);
        final int numColumnsInRow = classLabelsColumn < 0 ? words.length : words.length - 1;
        if(numColumns >= 0 && numColumnsInRow != numColumns) {
            throw new GATKException("filename does not encode a matrix, rows do not all have the same length");
        }
        final double[] features = new double[numColumnsInRow];
        int featureIndex = 0;
        for(int columnNumber = 0; columnNumber < words.length; ++columnNumber) {
            final String word = words[columnNumber];
            if(columnNumber == classLabelsColumn) {
                classLabelsList.add(Integer.valueOf(word));
            } else {
                features[featureIndex] = Double.valueOf(word);
                ++featureIndex;
            }
        }
        rowsList.add(features);
    }

    /**
     * Form submatrix specified by selecting specified rows
     */
    public static RealMatrix sliceRows(final RealMatrix matrix, final int[] sliceRows) {
        final int[] allColumns = getRange(matrix.getColumnDimension());
        return matrix.getSubMatrix(sliceRows, allColumns);
    }

    /**
     * Abstract base classifier class. Actual classifiers must extend GATKClassifier and wrap whatever 3rd-party package
     * implements the actual numerics. Note that implementing class must be Serializable and KryoSerializble (that is
     * how saving is implemented)
     */
    public static abstract class GATKClassifier implements Serializable, KryoSerializable {
        private static final long serialVersionUID = 1L;

        private final double[][] singlePredictWrapper = new double [1][];

        /* public abstract routines. Implementing classes must override these */

        /**
         * Train classifier.
         * This routine must be overridden by implementing class. It should update "this" to a trained state
         * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
         *                             initialization and/or training routine. Since map values are Objects, the
         *                             classifier must keep track of their type.
         * @param truthSet             Data with matrix of features paired to correct class labels.
         * @return upon success, the function should return "this"
         */
        public abstract GATKClassifier train(final Map<String, Object> classifierParameters,
                                             final TruthSet truthSet);

        /**
         * Train classifier and return trace of quality vs training round.
         * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
         *                             initialization and/or training routine. Since map values are Objects, the
         *                             classifier must keep track of their type.
         * @param trainingSet          Truth data used for updating classifier during training.
         * @param evaluationSet        Truth data used for evaluating classifier performance during training.
         * @param maxTrainingRounds    Train for no more than this many rounds.
         * @param earlyStoppingRounds  Training should stop if the best quality is more than earlyStoppingRounds ago.
         *                             If early stopping is triggered, remaining trace values should be set to last
         *                             obtained training value (i.e. not the best).
         * @return double[] with length maxTrainingRounds, each value storing classifier "quality" as a function of
         *         training round.
         */
        public abstract double[] trainAndReturnQualityTrace(
                final Map<String, Object> classifierParameters, final TruthSet trainingSet,
                final TruthSet evaluationSet, final int maxTrainingRounds, final int earlyStoppingRounds);

        /**
         * return numDataRows x numClasses double[][] with probability each data point is a member of each class
         */
        public abstract double[][] predictProbability(final RealMatrix matrix);

        /**
         * return true if the "quality" returned by trainAndReturnQualityTrace should be maximized, false if it should
         * be minimized
         */
        public abstract boolean getMaximizeEvalMetric(final Map<String, Object> classifierParameters);

        /* public defined methods. Can be overridden if more efficient routines are available for specific implementation */

        /**
         * Predict the probability that an individual data point is a member of each class
         */
        public double[] predictProbability(final double[] featureVector) {
            singlePredictWrapper[0] = featureVector;
            final RealMatrix matrix = new Array2DRowRealMatrix(singlePredictWrapper, false);
            return (predictProbability(matrix)[0]);
        }

        /**
         * Predict the class label for each data point (represented by matrix rows)
         */
        public int[] predictClassLabels(final RealMatrix matrix) {
            final double [][] predictedProbabilities = predictProbability(matrix);
            final int [] predictedLabels = new int[predictedProbabilities.length];
            if(predictedProbabilities.length == 0) {
                return predictedLabels;
            }
            final int numColumns = predictedProbabilities[0].length;
            if(numColumns == 1) {
                // binary classifier, reporting only probability of class == 1 (or "true")
                for(int row = 0; row < predictedProbabilities.length; ++row) {
                    predictedLabels[row] = predictedProbabilities[row][0] >= 0.5 ? 1 : 0;
                }

            } else {
                // multiclass classifier (or at binary independently reporting probability of class 0 or 1)
                for (int row = 0; row < predictedProbabilities.length; ++row) {
                    predictedLabels[row] = argmax(predictedProbabilities[row]);
                }
            }
            return predictedLabels;
        }

        /**
         * Save classifier to file at specified path (via Kryo)
         */
        public void save(final String saveFilePath) throws IOException {
            try(FileOutputStream fileOutputStream = new FileOutputStream(saveFilePath)) {
                save(fileOutputStream);
            }
        }

        /**
         * Save classifier to specified FileOutputStream (via Kryo)
         */
        public void save(final FileOutputStream fileOutputStream) {
            final Kryo kryo = new Kryo();
            final Output output = new Output(fileOutputStream);
            kryo.writeClassAndObject(output, this);
            output.close();
        }

        /**
         * Load classifier from file at specified path (via Kryo)
         */
        public static GATKClassifier load(final String saveFilePath) throws IOException {
            try(FileInputStream fileInputStream = new FileInputStream(saveFilePath)) {
                return load(fileInputStream);
            }
        }

        /**
         * Load classifier from specified FileInputStream (via Kryo)
         */
        public static GATKClassifier load(final FileInputStream fileInputStream) {
            final Kryo kryo = new Kryo();
            final Input input = new Input(fileInputStream);
            final GATKClassifier classifier = (GATKClassifier)kryo.readClassAndObject(input);
            input.close();
            return classifier;
        }

        /**
         * Obtain optimal classifierParameters. Get quality estimates by calling trainAndReturnQualityTrace() with
         * cross-validation. Select optimal number of training rounds (to be used when early stopping is not possible).
         * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
         *                             initialization and/or training routine. Since map values are Objects, the
         *                             classifier must keep track of their type. This must specify any parameters that
         *                             will not be tuned. If there is overlap in keys between classifierParameters and
         *                             tuneClassifierParameters, the tuned value will take precedence.
         * @param tuneClassifierParameters Map (from parameter names to ClassifierParamRange) specifying values that
         *                                 will be tuned, their range of allowed values, and how they are distributed
         *                                 (float vs integer, linear vs log).
         * @param classifierTuningStrategy enum specifying strategy for selecting candidate parameter values
         * @param truthSet                 Data with matrix of features paired to correct class labels.
         * @param random                   Random number generator.
         * @param stratify                 Vector of integers listing stratification class. Cross-validation will split
         *                                 data so that each stratify value has balanced number of instances in each
         *                                 fold. If null, split randomly. It is highly desirable to stratify at least
         *                                 based on class label.
         * @param numCrossvalidationFolds  Must be >= 2
         * @param maxTrainingRounds        Train classifiers with candidate hyperparameters for at most this many rounds.
         * @param earlyStoppingRounds      Stop training if the best training quality was this many rounds ago.
         * @param numTuningRounds          Try this many candidate hyperparameters before stopping and selecting the best.
         * @return bestHyperparameters     classifierParameters updated with optimal values from tuneClassifierParameters.
         */
        public Map<String, Object> tuneClassifierParameters(final Map<String, Object> classifierParameters,
                                                            final Map<String, ClassifierParamRange<?>> tuneClassifierParameters,
                                                            final ClassifierTuningStrategy classifierTuningStrategy,
                                                            final TruthSet truthSet, final Random random,
                                                            final int[] stratify, final int numCrossvalidationFolds,
                                                            final int maxTrainingRounds, final int earlyStoppingRounds,
                                                            final int numTuningRounds) {
            final List<TrainTestSplit> splits = new ArrayList<>(numCrossvalidationFolds);
            TrainTestSplit.getCrossvalidationSplits(
                    numCrossvalidationFolds, truthSet.getNumRows(), random, stratify
            ).forEachRemaining(splits::add);

            final ClassifierTuner classifierTuner;
            switch(classifierTuningStrategy) {
                case RANDOM:
                    classifierTuner = new RandomClassifierTuner(
                            this, classifierParameters, tuneClassifierParameters, truthSet,
                            splits, maxTrainingRounds, earlyStoppingRounds, numTuningRounds, random
                    );
                    break;
                default:
                    throw new IllegalStateException("Invalid ClassifierTuningStrategy: " + classifierTuningStrategy);
            }

            return classifierTuner.getBestParameters();

        }

        /**
         * Predict class labels using cross-validation, so classifier does not predict on data that it was trained on.
         * @param truthSet             Data with matrix of features paired to correct class labels.
         * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
         *                             initialization and/or training routine. Since map values are Objects, the
         *                             classifier must keep track of their type. This must specify any parameters that
         *                             will not be tuned. If there is overlap in keys between classifierParameters and
         *                             tuneClassifierParameters, the tuned value will take precedence.
         * @param random               Random number generator
         * @param stratify             Vector of integers listing stratification class. Cross-validation will split data
         *                             so that each stratify value has balanced number of instances in each fold. If
         *                             null, split randomly. It is highly desirable to stratify at least based on class
         *                             label.
         * @param numCrossvalidationFolds  >= 2
         */
        public int[] crossvalidatePredict(final TruthSet truthSet, final Map<String, Object> classifierParameters,
                                          final Random random, final int[] stratify, final int numCrossvalidationFolds) {
            final int[] predictedLabels = new int[truthSet.getNumRows()];
            final Iterator<TrainTestSplit> splitIterator = TrainTestSplit.getCrossvalidationSplits(
                    numCrossvalidationFolds, truthSet.getNumRows(), random, stratify
            );
            while(splitIterator.hasNext()) {
                final TrainTestSplit split = splitIterator.next();
                // train on training data from this crossvalidation split
                train(classifierParameters, truthSet.sliceRows(split.trainRows));
                // predict on testing data from this split
                final int[] predictedTestLabels = predictClassLabels(sliceRows(truthSet.features, split.testRows));
                // and assign those values into the final predictions
                sliceAssign(predictedLabels, split.testRows, predictedTestLabels);
            }
            return predictedLabels;
        }

        /**
         * Choose the number of threads to obtain fastest training times. Update classifierParameters with the
         * appropriate value for the number of threads.
         * @param classifierParameters Map (from parameter names to values) that will be passed to classifier
         *                             initialization and/or training routine. Since map values are Objects, the
         *                             classifier must keep track of their type. This must specify any parameters that
         *                             will not be tuned. If there is overlap in keys between classifierParameters and
         *                             tuneClassifierParameters, the tuned value will take precedence.
         * @param numThreadsKey  Key name in classifierParameters for number of threads.
         * @param truthSet       Data to obtain training times from. Ideally this data should have similar structure to
         *                       actual data (there is no reason you can't subsequently train on it). Note that a small
         *                       subset of the data will be selected, so that training times are short and this routine
         *                       does not waste a lot of time.
         */
        public void chooseNumThreads(final Map<String, Object> classifierParameters, final String numThreadsKey,
                                     final TruthSet truthSet) {
            final int numCalibrationRows = NUM_CALIBRATION_CLASS_ROWS * 2;
            localLogger.info("numCalibrationRows = " + numCalibrationRows);
            localLogger.info("numTrainingRows = " + truthSet.getNumRows());
            final TruthSet calibrationMatrix;
            if(truthSet.getNumRows() <= numCalibrationRows) {
                calibrationMatrix = truthSet;
            } else {
                final int[] stratify = truthSet.classLabels;
                localLogger.info("trainingFraction = " + numCalibrationRows / (double)truthSet.getNumRows());
                final TrainTestSplit trainTestSplit = TrainTestSplit.getTrainTestSplit(
                        numCalibrationRows / (double)truthSet.getNumRows(),
                        truthSet.getNumRows(), new Random(), stratify
                );
                localLogger.info("numTrain = " + trainTestSplit.trainRows.length);
                localLogger.info("numTest = " + trainTestSplit.testRows.length);
                calibrationMatrix = truthSet.sliceRows(trainTestSplit.trainRows);
            }
            final Map<String, Object> calibrationParams = new HashMap<>(classifierParameters);
            calibrationParams.put(NUM_TRAINING_ROUNDS_KEY, NUM_CALIBRATION_TRAINING_ROUNDS);
            final int maxNumThreads =
                    classifierParameters.containsKey(numThreadsKey) && (int)classifierParameters.get(numThreadsKey) > 0 ?
                            (int)classifierParameters.get(numThreadsKey) : Runtime.getRuntime().availableProcessors();
            if(maxNumThreads == 1) {
                classifierParameters.put(numThreadsKey, 1);
                return;
            }
            long bestElapsedTime = Long.MAX_VALUE;
            for(int numThreads = 1; numThreads < maxNumThreads; ++numThreads) {
                calibrationParams.put(numThreadsKey, numThreads);
                final long elapsedTime = getTrainingTime(calibrationParams, calibrationMatrix);
                if(elapsedTime < bestElapsedTime) {
                    bestElapsedTime = elapsedTime;
                    classifierParameters.put(numThreadsKey, numThreads);
                }
            }

        }

        /* private methods */

        private long getTrainingTime(final Map<String, Object> classifierParameters, final TruthSet truthSet) {
            final long startTime = System.nanoTime();
            train(classifierParameters, truthSet);
            return System.nanoTime() - startTime;
        }

        private double[][] getCrossvalidatedTrainingTraces(final Map<String, Object> classifierParameters,
                                                          final TruthSet truthSet, final List<TrainTestSplit> splits,
                                                          final int maxTrainingRounds, final int earlyStoppingRounds) {
            final int numCrossvalidationFolds = splits.size();
            double[][] trainingTraces = new double[numCrossvalidationFolds][];
            for(int fold = 0; fold < numCrossvalidationFolds; ++fold) {
                final TrainTestSplit split = splits.get(fold);
                trainingTraces[fold] = trainAndReturnQualityTrace(
                        classifierParameters, truthSet.sliceRows(split.trainRows), truthSet.sliceRows(split.testRows),
                        maxTrainingRounds, earlyStoppingRounds
                );
            }
            return trainingTraces;
        }

        private double getTrainingScore(final Map<String, Object> classifierParameters,
                                        final double[][] trainingTraces) {
            // find the index with the best total (i.e. mean) score across rounds. This yields the optimal number of
            // rounds of training
            final boolean maximizeEvalMetric = getMaximizeEvalMetric(classifierParameters);
            double bestTotalScore = maximizeEvalMetric ? Double.MIN_VALUE : Double.MAX_VALUE;
            int bestRoundIndex = -1;
            final int maxTrainingRounds = trainingTraces[0].length;
            for(int roundIndex = 0; roundIndex < maxTrainingRounds; ++roundIndex) {
                double roundScore = trainingTraces[0][roundIndex];
                for(int fold = 1; fold < trainingTraces.length; ++fold) {
                    roundScore += trainingTraces[fold][roundIndex];
                }
                if(maximizeEvalMetric ? roundScore > bestTotalScore : roundScore < bestTotalScore) {
                    // the score at this round of training is the best so far
                    bestTotalScore = roundScore;
                    bestRoundIndex = roundIndex;
                }
            }
            final int numTrainingRounds = bestRoundIndex + 1;

            // report the overall score for this set of parameters as the score of the *worst* trace at the optimal training
            // index. Selecting the worst trace demands high reliability from the classifier across similar data sets.
            double trainingScore = trainingTraces[0][bestRoundIndex];
            if(maximizeEvalMetric) {
                for (int fold = 1; fold < trainingTraces.length; ++fold) {
                    trainingScore = Math.min(trainingScore, trainingTraces[fold][bestRoundIndex]);
                }
            } else {
                for (int fold = 1; fold < trainingTraces.length; ++fold) {
                    trainingScore = Math.max(trainingScore, trainingTraces[fold][bestRoundIndex]);
                }
            }

            classifierParameters.put(NUM_TRAINING_ROUNDS_KEY, numTrainingRounds);
            return trainingScore;
        }
    }

    /**
     * Class that holds partitions of data sets (by row). trainRows and testRows are row indices into original data set,
     * appropriate for passing to sliceRows()
     */
    public static class TrainTestSplit {
        public final int[] trainRows;
        public final int[] testRows;

        TrainTestSplit(final int[] trainRows, final int[] testRows) {
            this.trainRows = trainRows;
            this.testRows = testRows;
        }


        /**
         * Static method to split data into two sets: training set and testing set.
         * @param trainingFraction proportion in [0, 1] of original data that should be placed in training set.
         * @param numRows          number of rows of original data
         * @param random           random number generator
         * @param stratify         Vector of integers listing stratification class. Data will be split so that each
         *                         stratify value has balanced number of instances in each set (e.g. proportion of
         *                         rows with stratify == 3 in training set will be approximately trainingFraction). If
         *                         null, split randomly. It is highly desirable to stratify at least based on class label.
         * @return TrainTestSplit with appropriate assigned indices.
         */
        public static TrainTestSplit getTrainTestSplit(final double trainingFraction, final int numRows,
                                                       final Random random, final int[] stratify) {
            if(stratify != null && stratify.length != numRows) {
                throw new IllegalArgumentException(
                        "stratify.length (" + stratify.length + ") != numRows (" + numRows + ")"
                );
            }

            final long numTrain = Math.round(numRows * trainingFraction);
            final long numTest = numRows - numTrain;
            if(numTrain <= 0) {
                if(trainingFraction < 0) {
                    throw new IllegalArgumentException(
                            "trainingFraction (" + trainingFraction + ") must be in range [0, 1]"
                    );
                }
                return new TrainTestSplit(new int[0], getRange(numRows));
            } else if(numTest <= 0) {
                if(trainingFraction > 1) {
                    throw new IllegalArgumentException(
                            "trainingFraction (" + trainingFraction + ") must be in range[0, 1]"
                    );
                }
                return new TrainTestSplit(getRange(numRows), new int[0]);
            }
            final int[] split_index_ordering = stratify == null ?
                    getRandomPermutation(random, numRows)
                    : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);
            final int[] trainRows = new int[(int)numTrain];
            final int[] testRows = new int[(int)numTest];
            int nextTrainInd = 1;
            int nextTestInd = 1;
            for(final int split_index : split_index_ordering) {
                if(numTrain * nextTestInd >= numTest * nextTrainInd) {
                    // training set gets next index
                    trainRows[nextTrainInd - 1] = split_index;
                    ++nextTrainInd;
                } else {
                    testRows[nextTestInd - 1] = split_index;
                    ++nextTestInd;
                }
            }
            Arrays.sort(trainRows);
            Arrays.sort(testRows);
            return new TrainTestSplit(trainRows, testRows);
        }

        /**
         * Static method to get iterator over partitions of data, for implementing k-fold cross-validation.
         * @param numCrossvalidationFolds >= 2
         * @param numRows          number of rows in original data set
         * @param random           random number generator
         * @param stratify         Vector of integers listing stratification class. Cross-validation will split data so
         *                         that each stratify value has balanced number of instances in each fold. If null,
         *                         split randomly. It is highly desirable to stratify at least based on class label.
         * @return
         */
        public static Iterator<TrainTestSplit> getCrossvalidationSplits(final int numCrossvalidationFolds, final int numRows,
                                                                        final Random random, final int[] stratify) {
            if(numCrossvalidationFolds < 2) {
                throw new IllegalArgumentException("numCrossvalidationFolds (" + numCrossvalidationFolds + ") must be >= 2");
            }
            if(stratify != null && stratify.length != numRows) {
                throw new IllegalArgumentException(
                        "stratify.length (" + stratify.length + ") != numRows (" + numRows + ")"
                );
            }
            final int[] split_index_ordering = stratify == null ?
                    getRandomPermutation(random, numRows)
                    : TrainTestSplit.getStratfiedIndexOrdering(random, stratify);

            return new FoldSplitIterator(split_index_ordering, numCrossvalidationFolds);
        }

        /**
         * Get ordering of data so that evenly distributing rows according to this order causes sampling that is
         * a) balanced in each stratify value
         * b) randomly-ordered within each stratify value
         */
        private static int[] getStratfiedIndexOrdering(final Random random, final int[] stratify) {
            /*
            logical (but memory inefficient) process
            1. make a random permutation, and use it to permute stratify
            final int[] permutation = getRange(stratify.length);
            final int[] permuted_stratify = slice(stratify, permutation);
            2. find the indices that would sort the permuted stratify array. The permutation ensures that entries with
               the same value in stratify will be in random order. However, they point to the wrong stratify values.
            final int[] permuted_sort_indices = argsort(permuted_stratify);
            3. unpermute the sort_indices to permute to the correct stratify values. Because the sort visited the values
               in random order, the ordering of stratify indices will be permuted (between indices pointing to equal
               stratify values).
            final int[] stratify_inds = slice(permutation, permuted_sort_indices);
            */
            final int[] permutation = getRandomPermutation(random, stratify.length);
            return slice(permutation, argsort(slice(stratify, permutation)));
        }

        private static class FoldSplitIterator implements Iterator<TrainTestSplit> {
            private final int[] split_index_ordering;
            private final int numFolds;
            private int fold;

            FoldSplitIterator(final int[] split_index_ordering, final int numFolds) {
                this.split_index_ordering = split_index_ordering;
                this.numFolds = numFolds;
                this.fold = 0;
            }

            @Override
            public boolean hasNext() {
                return fold < numFolds;
            }

            @Override
            public TrainTestSplit next() {
                final int numTest = 1 + (split_index_ordering.length - 1 - fold) / numFolds;
                final int numTrain = split_index_ordering.length - numTest;
                int[] testRows = new int[numTest];
                int[] trainRows = new int[numTrain];
                int trainIndex;
                int orderingIndex;
                if(fold > 0) {
                    for(trainIndex = 0; trainIndex < fold; ++trainIndex) {
                        trainRows[trainIndex] = split_index_ordering[trainIndex];
                    }
                    orderingIndex = trainIndex;
                } else {
                    orderingIndex = 0;
                    trainIndex = 0;
                }
                for(int testIndex = 0; testIndex < testRows.length; ++testIndex) {
                    testRows[testIndex] = split_index_ordering[orderingIndex];
                    final int orderingStop = Math.min(orderingIndex + numFolds, split_index_ordering.length);
                    for(++orderingIndex; orderingIndex < orderingStop; ++orderingIndex, ++trainIndex) {
                        trainRows[trainIndex] = split_index_ordering[orderingIndex];
                    }
                }

                ++fold;
                Arrays.sort(trainRows);
                Arrays.sort(testRows);
                return new TrainTestSplit(trainRows, testRows);
            }
        }
    }

    // To-do: write ClassifierTuner with strategy = BAYES
    public enum ClassifierTuningStrategy { RANDOM }

    /**
     * Abstract base class for tuning classifier hyperparameters. The base class manages keeping track of previously
     * seen / best hyperparameters and corresponding scores. Actual classifier tuners must extend ClassifierTuner
     * and override chooseNextHyperparameters(). Note, they will probably always be called by tuneClassifierParameters()
     * rather than instantiated directly.
     */
    static abstract class ClassifierTuner {
        private final GATKClassifier classifier;
        protected final Map<String, Object> classifierParameters;
        protected final Map<String, ClassifierParamRange<?>> tuneParameters;
        protected final List<Map<String, Object>> hyperparameterSets;
        protected final List<Double> hyperparameterScores;
        protected final boolean maximizeEvalMetric;

        // keep these as class members instead of locals so that chooseNextHyperparameters() can see them if it needs to
        protected Map<String, Object> bestParameters;
        protected double bestScore;
        protected int numTuningRounds;

        private final TruthSet truthSet;
        private final List<TrainTestSplit> splits;
        private final int maxTrainingRounds;
        private final int earlyStoppingRounds;

        ClassifierTuner(final GATKClassifier classifier, final Map<String, Object> classifierParameters,
                        final Map<String, ClassifierParamRange<?>> tuneParameters, final TruthSet truthSet,
                        final List<TrainTestSplit> splits, final int maxTrainingRounds, final int earlyStoppingRounds,
                        final int numTuningRounds) {
            if(numTuningRounds < 1) {
                throw new IllegalArgumentException("numTuningRounds (" + numTuningRounds + ") must be >= 1");
            }

            this.classifier = classifier;
            this.tuneParameters = tuneParameters;
            this.hyperparameterSets = new ArrayList<>();
            this.hyperparameterScores = new ArrayList<>();
            this.classifierParameters = classifierParameters;

            this.truthSet = truthSet;
            this.splits = splits;
            this.maxTrainingRounds = maxTrainingRounds;
            this.earlyStoppingRounds = earlyStoppingRounds;
            this.numTuningRounds = numTuningRounds;

            maximizeEvalMetric = classifier.getMaximizeEvalMetric(classifierParameters);
            bestParameters = null;
            bestScore = maximizeEvalMetric ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        }

        /**
         * Return next candidate set of hyperparameters
         */
        abstract protected Map<String, Object> chooseNextHyperparameters();

        /**
         * Try multiple candidate hyperparameters (chosen according to strategy in non-abstract ClassifierTuner derived
         * class). Keep track of best, and after numTuningRounds candidates have been tried, return the best hyperparameters
         */
        Map<String, Object> getBestParameters() {
            MachineLearningUtils.localLogger.info("Getting best parameters");
            ConsoleProgressBar progress = new ConsoleProgressBar(numTuningRounds);
            for(int i = 0; i < numTuningRounds; ++i) {
                final Map<String, Object> hyperparameters = chooseNextHyperparameters();
                hyperparameterSets.add(hyperparameters);
                final Map<String, Object> testParameters = new HashMap<>(classifierParameters);
                testParameters.putAll(hyperparameters);
                final double[][] trainingTraces = classifier.getCrossvalidatedTrainingTraces(
                        testParameters, truthSet, splits, maxTrainingRounds, earlyStoppingRounds
                );
                final double score = classifier.getTrainingScore(testParameters, trainingTraces);
                hyperparameterScores.add(score);
                // This is the new best score if a) it is better than the previous best score so far -OR-
                //                               b) it exactly ties the best score, but uses fewer rounds of training
                if((maximizeEvalMetric ? score > bestScore : score < bestScore)
                        || (score == bestScore &&
                            (int)testParameters.get(NUM_TRAINING_ROUNDS_KEY) < (int)bestParameters.get(NUM_TRAINING_ROUNDS_KEY))) {
                    bestScore = score;
                    bestParameters = testParameters;
                }
                progress.update(1);
            }
            return bestParameters;
        }
    }

    /**
     * Concrete ClassifierTuner that selects hyperparameters randomly. Specifically:
     * 1) each parameter value is sampled evenly (log or linear) along its range, with numTuningRounds samples.
     * 2) these samples are permuted into random order
     * 3) different parameter values are permuted independently
     * The result is an even sampling in the hyper-rectangle, but with less clumping. Unimportant parameters (e.g. those
     * that don't affect quality) will not decrease sample density in important parameters. Thus this method should
     * outperform grid search with an equivalent numTuningRounds.
     */
    static class RandomClassifierTuner extends ClassifierTuner {
        private final Map<String, Object[]> randomParameters;

        RandomClassifierTuner(final GATKClassifier classifier, final Map<String, Object> classifierParameters,
                              final Map<String, ClassifierParamRange<?>> tuneParameters, final TruthSet truthSet,
                              final List<TrainTestSplit> splits, final int maxTrainingRounds,
                              final int earlyStoppingRounds, final int numTuningRounds, final Random random) {
            super(classifier, classifierParameters, tuneParameters, truthSet, splits, maxTrainingRounds,
                    earlyStoppingRounds, numTuningRounds);
            randomParameters = tuneParameters.entrySet().stream().collect(
                    Collectors.toMap(
                            Map.Entry::getKey,
                            x -> x.getValue().getRandomSamples(random, numTuningRounds)
                    )
            );
        }

        @Override
        protected Map<String, Object> chooseNextHyperparameters() {
            final int index = hyperparameterScores.size();
            return randomParameters.entrySet().stream().collect(
                    Collectors.toMap(
                            Map.Entry::getKey,
                            x -> x.getValue()[index]
                    )
            );
        }

    }


    /**
     * Interface for classifier parameter range.
     */
    public interface ClassifierParamRange<T> {
        T[] getRandomSamples(final Random random, final int numSamples);
    }

    /**
     * Class for double-valued ClassifierParamRange with linear (uniform) sampling across range
     */
    public static class ClassifierLinearParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        public ClassifierLinearParamRange(double low, double high) {
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Double[] getRandomSamples(final Random random, final int numSamples) {
            final Double[] samples = new Double[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (high + low) / 2.0;
                }
                return samples;
            }
            final double delta = (high - low) / (numSamples - 1);
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val += delta;
                samples[permutation[i]] = val;
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Class for double-valued ClassifierParamRange with logarithmic sampling across range
     */
    public static class ClassifierLogParamRange implements ClassifierParamRange<Double> {
        private final double low;
        private final double high;

        public ClassifierLogParamRange(double low, double high) {
            if(low * high <= 0) {
                throw new IllegalArgumentException("low (" + low + ") and high (" + high + ") must be the same sign, and non-zero");
            }
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Double[] getRandomSamples(final Random random, final int numSamples) {
            final Double[] samples = new Double[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = Math.sqrt(high * low);
                }
                return samples;
            }
            final double delta = Math.pow(high / low, 1.0 / (numSamples - 1));
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val *= delta;
                samples[permutation[i]] = val;
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Class for integer-valued ClassifierParamRange with linear (uniform) sampling across range
     */
    public static class ClassifierIntegerLinearParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        public ClassifierIntegerLinearParamRange(int low, int high) {
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Integer[] getRandomSamples(final Random random, final int numSamples) {
            final Integer[] samples = new Integer[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (int)Math.round((high + low) / 2.0);
                }
                return samples;
            }
            final double delta = (high - low) / (double)(numSamples - 1);
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val += delta;
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Class for integer-valued ClassifierParamRange with logarithmic sampling across range
     */
    public static class ClassifierIntegerLogParamRange implements ClassifierParamRange<Integer> {
        private final int low;
        private final int high;

        public ClassifierIntegerLogParamRange(int low, int high) {
            if(low * high <= 0) {
                throw new IllegalArgumentException("low (" + low + ") and high (" + high + ") must be the same sign, and non-zero");
            }
            if(low > high) {
                throw new IllegalArgumentException("low must be <= high");
            }
            this.low = low;
            this.high = high;
        }

        public Integer[] getRandomSamples(final Random random, final int numSamples) {
            final Integer[] samples = new Integer[numSamples];
            if(numSamples < 2) {
                if(numSamples == 1) {
                    samples[0] = (int)Math.round(Math.sqrt((double)high * low));
                }
                return samples;
            }
            final double delta = Math.pow(high / (double)low, 1.0 / (numSamples - 1));
            double val = low;
            final int[] permutation = getRandomPermutation(random, numSamples);
            samples[permutation[0]] = low;
            // abort loop one step early and assign high to last element, to avoid round-off error potentially assigning
            // samples outside allowed range
            for(int i = 1; i < numSamples - 1; ++i) {
                val *= delta;
                samples[permutation[i]] = (int)Math.round(val);
            }
            samples[permutation[numSamples - 1]] = high;
            return samples;
        }
    }

    /**
     * Progress bar for console programs; to show work done, work to-do, elapsed time, and estimated remaining time for
     * long-running tasks. It avoids redrawing too frequently, so as to prevent slowing progress of actual work.
     */
    public static class ConsoleProgressBar {
        private static final long MIN_UPDATE_INTERVAL_NS = 500000000; // = 0.5 sec
        private static final int NUM_BAR_CHARACTERS = 10;

        private static final double SECONDS_IN_MINUTE = 60.0;
        private static final long MINUTES_IN_HOUR = 60;
        private static final long HOURS_IN_DAY = 24;
        private static final String CARRIAGE_RETURN = "\r";

        private final long workToDo;
        private final long bornTime;
        private long nextUpdateTime;
        private long workDone;
        private long workRemaining;
        private int maxOutputLength;
        private final String updateInfoFormat;

        /**
         * Create progress bar
         * @param workToDo sum of amount of work for all tasks. Usually this will just be number of tasks.
         */
        ConsoleProgressBar(final long workToDo) {
            if(workToDo <= 0) {
                throw new IllegalArgumentException("workToDo must be > 0");
            }
            this.workToDo = workToDo;
            workDone = 0;
            workRemaining = workToDo;
            bornTime = System.nanoTime();
            nextUpdateTime = bornTime + MIN_UPDATE_INTERVAL_NS;
            maxOutputLength = 0;
            final int workToDoLength = String.format("%d", workToDo).length();
            final String workDoneFormat = String.format("%%%dd/%%%dd", workToDoLength, workToDoLength);
            updateInfoFormat = workDoneFormat + " %.1f%% elapsed %s, remaining %s";
            drawBar(0.0, Double.NaN);
        }

        /**
         * Update the progress bar. Re-draw if enough time has elapsed.
         * @param workJustCompleted This value can be altered if different tasks have quantifiable differences in the
         *                          amount of work (or time) they will take. Probably this will normally just equal 1.
         */
        public void update(final long workJustCompleted) {
            if(workJustCompleted <= 0) {
                throw new IllegalArgumentException("workJustCompleted must be > 0");
            }
            workRemaining -= workJustCompleted;
            workDone += workJustCompleted;
            final long now = System.nanoTime();
            if(now < nextUpdateTime && workRemaining > 0) {
                return;  // avoid thrashing to the screen
            } else {
                nextUpdateTime = now + MIN_UPDATE_INTERVAL_NS;
            }
            final double elapsedTimeSec = 1.0e-9 * (now - bornTime);
            final double workPerSec = workDone / elapsedTimeSec;
            final double remainingTimeSec = workRemaining / workPerSec;
            System.out.flush();
            drawBar(elapsedTimeSec, remainingTimeSec);
        }

        private void drawBar(final double elapsedTimeSec, final double remainingTimeSec) {
            // NOTE: carriage return means each output will start from the beginning of the line
            final StringBuilder stringBuilder = new StringBuilder(CARRIAGE_RETURN);
            // draw actual progress bar
            final int numBarFilled = (int)(NUM_BAR_CHARACTERS * workDone / workToDo);
            final int numBarUnfilled = NUM_BAR_CHARACTERS - numBarFilled;
            stringBuilder.append("|");
            stringBuilder.append(StringUtils.repeat('#', numBarFilled));
            stringBuilder.append(StringUtils.repeat(' ', numBarUnfilled));
            stringBuilder.append("| ");
            // write summary statistics on completion amount, times
            stringBuilder.append(
                    workDone > 0 ?
                    String.format(
                        updateInfoFormat, workDone, workToDo, workDone * 100.0 / workToDo,
                        secondsToTimeString(elapsedTimeSec), secondsToTimeString(remainingTimeSec)
                    )
                    : String.format(updateInfoFormat, 0, workToDo, 0.0, secondsToTimeString(0.0), "???")
            );
            // do any necessary padding
            if(stringBuilder.length() > maxOutputLength) {
                maxOutputLength = stringBuilder.length();
            } else if(stringBuilder.length() < maxOutputLength) {
                // pad with spaces to obliterate previous message
                stringBuilder.append(StringUtils.repeat(' ', maxOutputLength - stringBuilder.length()));
            }
            // if we're done, add newline
            if(workRemaining <= 0) {
                stringBuilder.append("\n");
            }
            // write out and flush
            System.out.print(stringBuilder.toString());
            System.out.flush();
        }

        private static String secondsToTimeString(double seconds) {
            if(seconds < SECONDS_IN_MINUTE) {
                return String.format("%.1fs", seconds);
            }
            long minutes = (int)Math.floor(seconds / SECONDS_IN_MINUTE);
            seconds = seconds % SECONDS_IN_MINUTE;
            long hours = minutes / MINUTES_IN_HOUR;
            if(hours <= 0) {
                return String.format("%dm %.1fs", minutes, seconds);
            }
            minutes = minutes % MINUTES_IN_HOUR;
            long days = hours / HOURS_IN_DAY;
            if(days <= 0) {
                return String.format("%dh %dm %.1fs", hours, minutes, seconds);
            } else {
                hours = hours % HOURS_IN_DAY;
                return String.format("%dd %dh %dm %.1fs", days, hours, minutes, seconds);
            }
        }
    }

    /* here are multiple static functions that are generally useful for machine-learning tasks. Kept public to encourage re-use */

    /**
     * Return Integer array with elements {0, 1, 2, ..., numElements - 1}
     */
    public static Integer[] getRange(final Integer numElements) {
        if(numElements < 0) {
            throw new IllegalArgumentException("numElements must be >= 0");
        }
        final Integer[] range = new Integer[numElements];
        for(Integer i = 0; i < numElements; ++i) {
            range[i] = i;
        }
        return range;
    }

    /**
     * Return int array with elements {0, 1, 2, ..., numElements - 1}
     */
    public static int[] getRange(final int numElements) {
        if(numElements < 0) {
            throw new IllegalArgumentException("numElements must be >= 0");
        }
        return IntStream.range(0, numElements).toArray();
    }

    /**
     * Dereference array and return new array sampled from the array of supplied indices
     * @return arr[indices]
     */
    public static int[] slice(final int[] arr, final int[] indices) {
        final int[] sliced_arr = new int[indices.length];
        for(int i = 0; i < indices.length; ++i) {
            sliced_arr[i] = arr[indices[i]];
        }
        return sliced_arr;
    }

    /**
     * Dereference array and assign new values into the supplied indices: arr[indices] = newValues
     */
    public static void sliceAssign(final int[] arr, final int[] indices, final int[] newValues) {
        if(indices.length != newValues.length) {
            throw new IllegalArgumentException("length of indices does not match length of newValues");
        }
        for(int i = 0; i < indices.length; ++i) {
            arr[indices[i]] = newValues[i];
        }
    }

    /**
     * Find index to maximum value in array
     */
    public static int argmax(final double[] arr) {
        int bestIndex = -1;
        double bestValue = Double.MIN_VALUE;
        for(int index = 0; index < arr.length; ++index) {
            if(arr[index] > bestValue) {
                bestIndex = index;
                bestValue = arr[index];
            }
        }
        return bestIndex;
    }

    /**
     * Join to RealMatrices together so that columns of A come first, then columns of B.
     */
    public static RealMatrix concatenateColumns(final RealMatrix matrixA, final RealMatrix matrixB) {
        if(matrixA.getRowDimension() != matrixB.getRowDimension()) {
            throw new IllegalArgumentException("matrixA and matrixB do not have same number of rows.");
        }
        final RealMatrix matrixC = matrixA.createMatrix(
                matrixA.getRowDimension(),
                matrixA.getColumnDimension() + matrixB.getColumnDimension()
        );
        matrixC.setSubMatrix(matrixA.getData(), 0, 0);
        matrixC.setSubMatrix(matrixB.getData(), 0, matrixA.getColumnDimension());
        return matrixC;
    }

    /**
     * Return int array of indices that would sort supplied array, e.g. in pseudocode: arr[argsort(arr)] == sort(arr)
     */
    public static int[] argsort(final int[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparingInt(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    /**
     * Return int array of indices that would sort supplied array, e.g. in pseudocode: arr[argsort(arr)] == sort(arr)
     */
    public static <T extends Comparable<? super T>> int[] argsort(final T[] arr) {
        final Integer[] sortIndices = getRange((Integer)arr.length);
        Arrays.sort(sortIndices, Comparator.comparing(ind -> arr[ind]));
        return ArrayUtils.toPrimitive(sortIndices);
    }

    /**
     * Return random permutation of elements {0, 1, ..., numElements - 1}, using Knuth shuffle algorithm.
     */
    public static int[] getRandomPermutation(final Random random, final int numElements) {
        // Knuth shuffle
        final int[] permutation = getRange(numElements);
        for(int i = numElements - 1; i > 0; --i) {
            final int swap_ind = random.nextInt(i);
            final int swap_val = permutation[swap_ind];
            permutation[swap_ind] = permutation[i];
            permutation[i] = swap_val;
        }
        return permutation;
    }

    /**
     * Geven predicted and correct labels, return proportion of predictions that are correct
     */
    public static double getPredictionAccuracy(final int[] predictedLabels, final int[] correctLabels) {
        int numCorrect = 0;
        for(int row = 0; row < correctLabels.length; ++row) {
            if(predictedLabels[row] == correctLabels[row]) {
                ++numCorrect;
            }
        }
        return numCorrect / (double)correctLabels.length;
    }

    /**
     * Convert stratify matrix to stratify array. Columns are binned and stratify array values are obtained by finding
     * unique rows.
     * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
     * @param stratifyMatrix matrix to convert to stratify array
     * @param numBins Number of bins into which to group values in each column
     * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
     *                                  the number of columns in consideration is reduced and the under-count rows are
     *                                  searched again for unique values.
     * @return stratifyArray, int array with length = number of rows in stratifyMatrix
     */
    public static int[] stratifyMatrixToStratifyArray(final RealMatrix stratifyMatrix, final int numBins,
                                                      final int minCountsPerStratifyValue) {
        return stratifyMatrixToStratifyArray(stratifyMatrix, numBins, minCountsPerStratifyValue, new HashSet<>());
    }

    /**
     * Convert stratify matrix to stratify array. Columns are binned and stratify array values are obtained by finding
     * unique rows. This version can specify some columns as categorical (unbinnable).
     * The stratify array is useful for balancing data when splitting into test/train sets or cross-validating.
     * @param stratifyMatrix matrix to convert to stratify array
     * @param numBins Number of bins into which to group values in each column
     * @param minCountsPerStratifyValue Minimum number of instances of each unique row. If any rows occur too few times,
     *                                  the number of columns in consideration is reduced and the under-count rows are
     *                                  searched again for unique values.
     * @param categoricalColumns Columns in categoricalColumns are not binned.
     * @return stratifyArray, int array with length = number of rows in stratifyMatrix
     */
    public static int[] stratifyMatrixToStratifyArray(final RealMatrix stratifyMatrix, final int numBins,
                                                      final int minCountsPerStratifyValue,
                                                      final Set<Integer> categoricalColumns) {
        // first form binned version of each column
        final int[][] binnedStratifyMatrix = new int[stratifyMatrix.getRowDimension()][stratifyMatrix.getColumnDimension()];
        for(int columnIndex = 0; columnIndex < stratifyMatrix.getColumnDimension(); ++columnIndex) {
            final int[] columnBins = categoricalColumns.contains(columnIndex) ?
                    getCategoryCodes(stratifyMatrix.getColumn(columnIndex))
                    : getBinnedColumn(stratifyMatrix.getColumn(columnIndex), numBins);
            for(int row = 0; row < stratifyMatrix.getRowDimension(); ++row) {
                binnedStratifyMatrix[row][columnIndex] = columnBins[row];
            }
        }

        // Get a set of unique rows
        Collection<List<Integer>> resultsToCheck = getUniqueRows(
                binnedStratifyMatrix,
                IntStream.range(0, stratifyMatrix.getRowDimension()).boxed().collect(Collectors.toList()),
                stratifyMatrix.getColumnDimension()
        );
        // While there are rows that have too few instances, decrease the number of columns under consideration (only
        // for those undersized rows)
        final Set<List<Integer>> uniqueResults = new HashSet<>();
        for(int useNumColumns = stratifyMatrix.getColumnDimension() - 1; useNumColumns > 0; --useNumColumns) {
            // Add undersized rows to a set to be reprocessed. Add properly sized rows to final uniqueResults set.
            final Set<List<Integer>> undersizedStratify = new HashSet<>();
            for(final List<Integer> uniqueResult : resultsToCheck) {
                if(uniqueResult.size() < minCountsPerStratifyValue) {
                    undersizedStratify.add(uniqueResult);
                } else {
                    uniqueResults.add(uniqueResult);
                }
            }
            if(undersizedStratify.isEmpty()) {
                break; // no problematic values, we're done
            } else if(useNumColumns == 1) {
                // can't consider fewer columns, lump remaining rows into the stratify value with the fewest members, to
                // make an "odd-ball" value.
                if(uniqueResults.isEmpty()) {
                    // handle edge case where every unique row has too few counts
                    uniqueResults.add(new ArrayList<>());
                }
                final List<Integer> smallest = uniqueResults.stream().min(Comparator.comparing(List::size))
                        .orElseThrow(NoSuchElementException::new);
                for(final List<Integer> uniqueResult : undersizedStratify) {
                    smallest.addAll(uniqueResult);
                }
                Collections.sort(smallest);
            } else {
                // find unique rows from subset that are too small, looking at one fewer column
                final List<Integer> tooSmall = undersizedStratify.stream().flatMap(Collection::stream).sorted()
                        .collect(Collectors.toList());
                resultsToCheck = getUniqueRows(binnedStratifyMatrix, tooSmall, useNumColumns);
            }
        }


        final int[] stratifyArray = new int[binnedStratifyMatrix.length];
        int stratifyValue = 0;
        for(final List<Integer> uniqueResult : uniqueResults) {
            for(final int index : uniqueResult) {
                stratifyArray[index] = stratifyValue;
            }
            ++stratifyValue;
        }
        return stratifyArray;
    }

    /**
     * Find unique values in column array, and return integer codes corresponding to those unique values
     * @return int array of length column.length, with values in range {0, 1, ..., numUniqueValues - 1}
     */
    private static int[] getCategoryCodes(final double[] column) {
        final int[] categoryCodes = new int[column.length];
        final Map<Double, Integer> codesMap = new HashMap<>();
        for(int i = 0; i < column.length; ++i) {
            final double value = column[i];
            final Integer code = codesMap.getOrDefault(value, null);
            if(code == null) {
                categoryCodes[i] = codesMap.size();
                codesMap.put(value, codesMap.size());
            } else {
                categoryCodes[i] = code;
            }
        }
        return categoryCodes;
    }

    /**
     * Bin column by
     * 1) choosing bins by sampling evenly in percentile space
     * 2) using bisection to map each value into its appropriate bin.
     * If NaN values are present, they are placed in highest bin.
     * @param column array of values to bin
     * @param numBins number of bins to distribute values into (e.g. number of unique integer values of output array)
     * @return int array with length column.length and values in {0, 1, ..., numBins-1}
     */
    private static int[] getBinnedColumn(final double[] column, final int numBins) {
        final int numUniqueColumnValues = (int)DoubleStream.of(column).distinct().count();
        if(numUniqueColumnValues <= numBins) {
            // Too much repetition to bin the data into the requested number of bins, just return unique values.
            return getCategoryCodes(column);
        }

        // Attempt to get requested number of percentiles, insisting that all percentiles are unique. If some values are
        // repeated, "percentiles" will not be evenly spaced, and it may not return exactly the requested number.
        // NOTE: the last percentile will not be used for binning (percentiles are "posts", bins are "fence") so request
        // one more percentile than bins
        final double[] percentiles = getUniquePercentiles(column, numBins + 1);

        // bin column to specified percentiles, dumping NaN into the last bin if it is present. Note that the percentiles
        // will be properly sized to account for the presence of NaN values.
        // NOTE: binarySearch is limited to ignore last percentile, because we don't want the max value in the array
        // being mapped past the maximum requested number of bins.
        return DoubleStream.of(column).mapToInt(
                v -> Double.isNaN(v) ? numBins - 1 : Arrays.binarySearch(percentiles, 0, percentiles.length, v)
        ).toArray();
    }

    /**
     * Return samples from Percentile object that are evenly-distributed in percentile space.
     * @param column array with values to sample from
     * @param numPercentiles Number of values to sample. If fewer than requested unique percentiles are found (if some
     *                       values repeat), sampling density will be increased to attempt to find the requested number.
     * @return array of length <= numPercentiles with values drawn from column that are approximately even in percentile
     *         space.
     */
    private static double[] getUniquePercentiles(final double[] column, final int numPercentiles) {
        final int numNaN = (int)Arrays.stream(column).filter(Double::isNaN).count();
        final int numSortablePercentiles = numNaN == 0 ? numPercentiles : numPercentiles - 1;
        final Percentile percentileEvaluator = new Percentile().withEstimationType(Percentile.EstimationType.R_1);
        percentileEvaluator.setData(column);
        double[] percentiles = percentileSpace(percentileEvaluator, numSortablePercentiles);
        if(percentiles.length == numSortablePercentiles) {
            return percentiles; // should be the case for typical non-repeating data
        }

        final int numSortable = column.length - numNaN;
        int numRequestHigh = numSortablePercentiles;
        while(percentiles.length < numSortablePercentiles && numRequestHigh < numSortable) {
            numRequestHigh = Math.min(2 * numRequestHigh, numSortable);
            percentiles = percentileSpace(percentileEvaluator, numRequestHigh);
        }
        if(percentiles.length == numSortablePercentiles) {
            return percentiles;
        }
        int numRequestLow = numSortablePercentiles;
        int range = numRequestHigh - numRequestLow;
        while(range > 1) {
            final int numRequest = numRequestLow + range / 2;
            percentiles = percentileSpace(percentileEvaluator, numRequest);
            if(percentiles.length < numSortablePercentiles) {
                numRequestLow = numRequest;
            } else if(percentiles.length > numSortablePercentiles){
                numRequestHigh = numRequest;
            } else {
                return percentiles;
            }
            range = numRequestHigh - numRequestLow;
        }

        // I'm not convinced that it's possible to get down to this point. It implies that you checked for one more
        // percentile but got 2 or more new unique values back. Just in case, insist on having *fewer* than requested,
        // since other approximate cases always yield that outcome.
        percentiles = percentileSpace(percentileEvaluator, numRequestLow);

        return percentiles;
    }

    /**
     * Return samples from Percentile object that are evenly-distributed in percentile space.
     * @param percentileEvaluator Percentile object that previously had setData() called to assign data.
     * @param numPercentiles >= 2
     */
    private static final double[] percentileSpace(final Percentile percentileEvaluator, final int numPercentiles) {
        if(numPercentiles < 2) {
            throw new IllegalArgumentException("numPercentiles must be >= 2");
        }
        final double low = 50.0 / percentileEvaluator.getData().length;
        final double high = 100.0 - low;
        final double coef = (high - low) / (numPercentiles - 1);
        return IntStream.range(0, numPercentiles).mapToDouble(i -> percentileEvaluator.evaluate(low + i * coef))
                .distinct().toArray();
    }


    /**
     * Find unique rows, and return set with indices to rows in original data
     * @param binnedStratifyMatrix  int[][] representation of data matrix
     * @param checkRows       only look for unique rows indexed by checkRows
     * @param useNumColumns   when checking for uniqueness, only consider columns in {0, 1, , ..., numUseColumns - 1}
     * @return set where each element is a list of row indices (from checkRows) so that for each row in the list,
     *         binnedStratifyMatrix[row] is identical.
     */
    private static Collection<List<Integer>> getUniqueRows(final int[][] binnedStratifyMatrix,
                                                    final List<Integer> checkRows,
                                                    final int useNumColumns) {
        // Use columns 0 to numUseColumns - 1, convert to string, and use as a key to collect indices into lists of unique rows
        return checkRows.stream().collect(Collectors.groupingBy(
                rowIndex -> Arrays.toString(Arrays.copyOfRange(binnedStratifyMatrix[rowIndex], 0, useNumColumns))
        )).values();
    }
}
