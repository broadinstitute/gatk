package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.MachineLearningUtils.DEFAULT_NUM_STRATIFY_BINS;


/**
 * <p>(Internal) Demo of XGBoostUtils + MachineLearningUtils: trains classifier on data in .csv file and saves classifier to binary file</p>
 *
 * <p><h3>Inputs</h3>
 * <ul>
 *     <li>A csv file of numeric data, with first column being class labels.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A binary file that stores a trained classifier.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk ExampleTrainXGBoostClassifier \
 *     -I input_data.csv \
 *     -O classifier.model
 * </pre></p>
 *
 * <p><h3>Caveats</h3>
 * <li>This tool uses xgboost's multi-threading. The user can pass --nthread to set the number of threads.</li>
 * <li>If --auto-select-nthread is passed, then nthread sets the upper bound on threads to use, and the actual value
 * will be obtained by finding the best performance when repeatedly solving a small test problem using varying numbers
 * of threads.</li>
 * <li>Because this tool is multi-threaded, it must be used on a single computer (not a multi-computer cluster).
 * However, the classifiers are KryoSerializable, so trained classifiers can be used to predict probabilities or class
 * labels in a spark environment.</li>
 * Coverage much lower than that probably won't work well.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Demo of XGBoostUtils + MachineLearningUtils: trains classifier on data in .csv file and saves classifier to binary file",
        summary = "This tool tunes hyperparameters, assesses classifier quality, and trains a classifier to predict class" +
                        " membership from data provided in a .csv file. The trained classifier is saved to a binary model file.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class ExampleTrainXGBoostClassifier extends CommandLineProgram {
    private static final long serialVersionUID = 1L;
    private static final Logger localLogger = LogManager.getLogger(ExampleTrainXGBoostClassifier.class);

    @Argument(doc = "path to .csv file storing data to train classifier. It is expected that the file will have all numeric" +
            "values, and the first column will contain class label.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    private String demoDataFile;

    @Argument(doc = "full path to save output (binary) classifier model file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private String classifierModelFile;

    @Argument(doc="Stop classifier training if score does not improve for this many consecutive rounds.",
            fullName = "early-stopping-rounds", optional = true)
    private final int earlyStoppingRounds = XGBoostUtils.DEFAULT_EARLY_STOPPING_ROUNDS;

    @Argument(doc="Train classifier for at most this many rounds.",
            fullName = "max-training-rounds", optional = true)
    private final int maxTrainingRounds = XGBoostUtils.DEFAULT_NUM_TRAINING_ROUNDS;

    @Argument(doc="When performing cross-validation, use this many folds.",
            fullName = "num-crossvalidation-folds", optional = true)
    private final int numCrossvalidationFolds = MachineLearningUtils.DEFAULT_NUM_CROSSVALIDATION_FOLDS;

    @Argument(doc="When optimizing hyperparameters, search for this many rounds.",
            fullName = "num-hyperparameter-optimization-rounds", optional = true)
    private final int numTuningRounds = MachineLearningUtils.DEFAULT_NUM_TUNING_ROUNDS;

    @Argument(doc="When optimizing hyperparameters, reserve this proportion of data for tuning hyperparameters.",
            fullName = "hyperparameter-tuning-proportion", optional = true)
    private Double hyperparameterTuningProportion = null;

    @Argument(doc="Use this metric to evaluate performance of the classifier.",
            fullName = "eval-metric", optional = true)
    private final String evalMetric = XGBoostUtils.DEFAULT_EVAL_METRIC;

    @Argument(doc="Seed for random numbers. If null, initialize randomly",
            fullName = "random-seed", optional = true)
    private final Long seed = XGBoostUtils.DEFAULT_SEED;

    @Advanced
    @Argument(doc="Auto-select optimal number of threads to use for training classifier? If false, use nthread as"
                  + " specified by user. If false, interpret that number as an upper bound, and select number of threads"
                  + " with maximal throughput on small test problem",
              fullName = "auto-select-nthread", optional = true)
    private final Boolean autoSelectNumThreads = false;

    @Argument(doc="Number of threads to use for training classifier. If ",
            fullName = "nthread", optional = true)
    private final int numThreads = autoSelectNumThreads ? Runtime.getRuntime().availableProcessors()
            : (int)Math.round(XGBoostUtils.GUESS_OPTIMAL_NUM_THREADS_PROPORTION * Runtime.getRuntime().availableProcessors());

    @Argument(doc="Tuning strategy for choosing classifier hyperparameters",
            fullName = "classifier-tuning-strategy", optional = true)
    final MachineLearningUtils.ClassifierTuningStrategy classifierTuningStrategy
            = MachineLearningUtils.ClassifierTuningStrategy.RANDOM;

    private final Random random = (seed == null ? new Random() : new Random(seed));
    /**
     * Demo stuff ENDS here....
     */

    @Override
    protected Object doWork() {
        localLogger.info("Loading demo data");
        final MachineLearningUtils.TruthSet truthSet = MachineLearningUtils.loadCsvFile(demoDataFile);

        // tune hyperparameters, assess cross-validated classifier accuracy, train final classifier
        final MachineLearningUtils.GATKClassifier classifier = trainClassifierAndAssessAccuracy(truthSet);

        localLogger.info("Saving final classifier to " + classifierModelFile);
        try {
            classifier.save(classifierModelFile);
        } catch(IOException err) {
            throw new GATKException(err.getClass() + ": " + err.getMessage());
        }

        // load trained classifier from disk, and use it to predict class membership probabilities, or class labels
        loadClassifierAndPredictClass(truthSet.features);

        return classifier;
    }

    /**
     * Demo tuning hyperparameters, assessing classifier accuracy, and training final classifier
     */
    private MachineLearningUtils.GATKClassifier trainClassifierAndAssessAccuracy(final MachineLearningUtils.TruthSet truthSet) {
        if(hyperparameterTuningProportion == null) {
            // This proportion produces the fastest training times overall, not necessarily the most useful results.
            hyperparameterTuningProportion = 1.0 / (1.0 + numTuningRounds);
        }

        // Choose min counts per stratify value so that after train-test split, each fold of cross-validation has at least
        // one element of each stratify value.
        final int minCountsPerStratifyValue = (hyperparameterTuningProportion == 0 || hyperparameterTuningProportion == 1) ?
                numCrossvalidationFolds
                : (int)Math.ceil(numCrossvalidationFolds / Math.min(hyperparameterTuningProportion, 1.0 - hyperparameterTuningProportion));
        localLogger.info("Stratifying data matrix to balance data splitting / cross-validation folds.");
        final int[] stratify = truthSet.getStratifyArray(
                 DEFAULT_NUM_STRATIFY_BINS, minCountsPerStratifyValue
        );


        localLogger.info("Splitting data");
        final MachineLearningUtils.TrainTestSplit hyperSplit = MachineLearningUtils.TrainTestSplit.getTrainTestSplit(
                hyperparameterTuningProportion, truthSet.getNumRows(), random, stratify
        );
        final MachineLearningUtils.TruthSet tuneSet = truthSet.sliceRows(hyperSplit.trainRows);
        final int[] tuneStratify = MachineLearningUtils.slice(stratify, hyperSplit.trainRows);
        final MachineLearningUtils.TruthSet validateSet = truthSet.sliceRows(hyperSplit.testRows);
        final int[] validateStratify = MachineLearningUtils.slice(stratify, hyperSplit.testRows);

        localLogger.info("Creating classifier object");
        final XGBoostUtils.GATKXGBooster classifier = new XGBoostUtils.GATKXGBooster();

        // set the number of threads
        final Map<String, Object> classifierParameters = new HashMap<> (XGBoostUtils.DEFAULT_CLASSIFIER_PARAMETERS);
        classifierParameters.put(XGBoostUtils.NUM_THREADS_KEY, numThreads);
        classifierParameters.put(XGBoostUtils.EVAL_METRIC_KEY, evalMetric);
        if(autoSelectNumThreads) {
            classifier.chooseNumThreads(classifierParameters, XGBoostUtils.NUM_THREADS_KEY, tuneSet);
            localLogger.info("Chose " + classifierParameters.get(XGBoostUtils.NUM_THREADS_KEY) + " threads");
        } else {
            localLogger.info("User selected " + classifierParameters.get(XGBoostUtils.NUM_THREADS_KEY) + " threads");
        }

        // Note: tuning hyperparameters is basically always necessary: if data has changed enough to merit training a
        // new classifier, then it merits finding the new best hyperparameters
        final Map<String, Object> bestClassifierParameters;
        if(tuneSet.getNumRows() > 0) {
            localLogger.info("Tuning hyperparameters");
            bestClassifierParameters = classifier.tuneClassifierParameters(
                    classifierParameters, XGBoostUtils.DEFAULT_TUNING_PARAMETERS, classifierTuningStrategy,
                    tuneSet, random, tuneStratify, numCrossvalidationFolds, maxTrainingRounds, earlyStoppingRounds,
                    numTuningRounds
            );
            localLogger.info("bestClassifierParameters: " + bestClassifierParameters.toString());
        } else {
            localLogger.info("skipping tuning hyperparameters, using default parameters");
            bestClassifierParameters = classifierParameters;
        }

        // note cross-validation is necessary if you want to estimate the accuracy of the new classifier. If you don't
        // care (e.g. are committed to using the classifier regardless) then you can set hyperparameterTuningProportion
        // to 1.0, or skip splitting the data and just tune on the whole set.
        if(validateSet.getNumRows() > 0) {
            localLogger.info("Cross-val predicting");
            final int[] predictedTestLabels = classifier.crossvalidatePredict(
                    validateSet, bestClassifierParameters, random, validateStratify, numCrossvalidationFolds
            );

            final double accuracy = MachineLearningUtils.getPredictionAccuracy(
                    predictedTestLabels, validateSet.classLabels
            );
            localLogger.info("Crossvalidated accuracy = " + String.format("%.1f%%", 100.0 * accuracy));
        } else {
            localLogger.info("skipping evaluating crossvalidated accuracy");
        }

        localLogger.info("Training final classifier");
        classifier.train(bestClassifierParameters, truthSet);
        return classifier;
    }


    /**
     * Demo loading saved classifier and using it to make predictions on data.
     */
    private void loadClassifierAndPredictClass(final RealMatrix dataMatrix) {
        localLogger.info("Re-loading saved classifier");
        final MachineLearningUtils.GATKClassifier loadedClassifier;
        try {
            loadedClassifier = MachineLearningUtils.GATKClassifier.load(classifierModelFile);
        } catch(IOException err) {
            throw new GATKException(err.getClass() +": " + err.getMessage());
        }

        // These aren't used here, but presumably any real application would do something based on probability or label
        localLogger.info("predicting probability that data is in each class");
        final double[][] probabilities = loadedClassifier.predictProbability(dataMatrix);

        localLogger.info("predicting class labels");
        final int[] classLabels = loadedClassifier.predictClassLabels(dataMatrix);
    }
}
