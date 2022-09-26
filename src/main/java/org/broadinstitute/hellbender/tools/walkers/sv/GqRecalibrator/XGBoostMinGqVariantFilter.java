package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import ml.dmlc.xgboost4j.LabeledPoint;
import ml.dmlc.xgboost4j.java.*;
import ml.dmlc.xgboost4j.java.util.BigDenseMatrix;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.Filter;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

@CommandLineProgramProperties(
        summary = "Extract matrix of properties for each variant. Also extract, num_variants x num_trios x 3 tensors of" +
                "allele count and genotype quality. These data will be used to train a variant filter based on min GQ" +
                "(and stratified by other variant properties) that maximizes the admission of variants with Mendelian" +
                "inheritance pattern while omitting non-Mendelian variants." +
                "Derived class must implement abstract method train_filter()",
        oneLineSummary = "Extract data for training min GQ variant filter from .vcf and .ped file with trios.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class XGBoostMinGqVariantFilter extends MinGqVariantFilterBase {
    @Argument(fullName="max-training-rounds", shortName="tr", doc="Maximum number of rounds of training", optional=true, minValue=1)
    public int maxTrainingRounds = 100;
    @Argument(fullName="early-stopping-rounds", shortName="e", doc="Stop training if no improvement is made in validation set for this many rounds. Set <= 0 to disable.", optional=true)
    public int earlyStoppingRounds = 10;
    @Argument(fullName="learning-rate", shortName="lr", doc="Learning rate for xgboost", optional=true)
    public double eta = 0.1;
    @Argument(fullName="max-depth", shortName="d", doc="Max depth of boosted decision tree", optional=true, minValue=1)
    public int maxDepth = 10;
    @Argument(fullName="gamma", doc="Regularization factor for xgboost", optional=true)
    public double gamma = 1.0e-3;
    @Argument(fullName="subsample", doc="Proportion of data selected for each tree", optional=true)
    public double subsample = 0.9;
    @Argument(fullName="colsample-by-tree", doc="Proportion of columns selected for each tree", optional=true)
    public double colsampleByTree = 0.9;
    @Argument(fullName="colsample-by-level", doc="Proportion of columns selected for each level of each tree", optional=true)
    public double colsampleByLevel = 1.0;
    @Argument(fullName="min-child-weight", doc="Minimum sum of hessian weight for each node of tree", optional=true, minValue=0.0)
    public double minChildWeight = 10000;
    @Argument(fullName="n-threads-xgboost", doc="Number of threads to use for xgboost training", optional=true)
    public int numThreadsXgBoost = -1;
    @Argument(fullName="temp-dir", doc="path to directory for temporary files", optional=true)
    File tempDir = null;
    enum EarlyStoppingLoss {Training, Diagnostic}
    @Argument(fullName="early-stopping-loss", doc="Stop early based on non-decreasing loss in either: \"Training\" or \"Diagnostic\"", optional=true)
    public EarlyStoppingLoss earlyStoppingLoss = EarlyStoppingLoss.Diagnostic;

    private Booster booster = null;

    private static final String TRAIN_MAT_KEY = "train";
    private static final String VALIDATION_MAT_KEY = "validation";

    private enum DMatrixFillMethod {RowMajorArray, BigDenseMatrix, LabeledPointIterator}
    @Argument(fullName="dmatrix-fill-method", doc="Method of filling dMatrix: \"RowMajorArray\", \"BigDenseMatrix\" or \"LabeledPointIterator\"")
    private final DMatrixFillMethod dMatrixFillMethod = DMatrixFillMethod.LabeledPointIterator;

    private DMatrix getDMatrixRowMajorArray(int[] variantIndices) throws XGBoostError {
        // Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
        final int numRows = getNumTrainableSampleVariants(variantIndices);

        final int numProperties = getNumProperties();
        final float[] rowMajorVariantProperties = new float[numRows * numProperties];
        // Loop over variants and filterable samples. Store properties for each sample in a single flat array
        final int numSamples = getNumSamples();
        int flatIndex = 0;
        for(final int variantIndex : variantIndices) {
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                flatIndex = propertiesTable.copyPropertiesRow(
                    rowMajorVariantProperties, flatIndex, variantIndex, sampleIndex, needsNormalizedProperties()
                );
            }
        }

        if(progressVerbosity > 1) {
            System.out.format("\t\tgot %d property entries\n", rowMajorVariantProperties.length);
        }

        return new DMatrix(
            rowMajorVariantProperties, numRows, numProperties, Float.NaN
        );
    }

    private DMatrix getDMatrixBigDenseMatrix(int[] variantIndices) throws XGBoostError {
        // Get number of rows, account for the fact that unfilterable (e.g. already HOMREF) samples will not be used
        final int numRows = getNumTrainableSampleVariants(variantIndices);

        final int numProperties = getNumProperties();
        final BigDenseMatrix variantSampleProperties = new BigDenseMatrix(numRows, numProperties);
        // Loop over variants and filterable samples. Store properties for each sample in a single flat array
        final int numSamples = getNumSamples();
        int row = 0;
        for(final int variantIndex : variantIndices) {
            for(int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
                if(!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    continue;
                }
                final float[] sampleProperties = propertiesTable.getPropertiesRow(variantIndex, sampleIndex,
                                                                                  needsNormalizedProperties());
                for(int column = 0; column < sampleProperties.length; ++column) {
                    variantSampleProperties.set(row, column, sampleProperties[column]);
                }
                ++row;
            }
        }
        if(progressVerbosity > 1) {
            System.out.format("\t\tgot %d property entries\n",
                    variantSampleProperties.nrow * (long)variantSampleProperties.ncol);
        }

        return new DMatrix(variantSampleProperties, Float.NaN);
    }

    private int fillVariantWeights(final int variantIndex, final int numSamples,
                                   final boolean[] sampleVariantTruth, int row, final float[] weights) {
        final int propertyBin = propertyBins[variantIndex];
        final Set<Integer> inheritanceTrainableSampleIndices = getInheritanceTrainableSampleIndices(variantIndex);
        final Set<Integer> truthTrainableSampleIndices = new HashSet<>(
                getGoodSampleIndices(variantIndex)
        );
        truthTrainableSampleIndices.addAll(
                getBadSampleIndices(variantIndex)
        );

        final float minGqWeight = propertyBinMinGqWeights[propertyBin];
        final float goodTruthWeight = minGqWeight + propertyBinGoodTruthWeights[propertyBin];
        final float badTruthWeight = minGqWeight + propertyBinBadTruthWeights[propertyBin];
        final float truthWeight = minGqWeight + (float)propertyBinTruthWeights[propertyBin];
        final float goodInheritanceWeight = minGqWeight + propertyBinGoodInheritanceWeights[propertyBin];
        final float badInheritanceWeight = minGqWeight + propertyBinBadInheritanceWeights[propertyBin];
        final float inheritanceWeight = minGqWeight + (float)propertyBinInheritanceWeights[propertyBin];
        final float goodBothWeights = minGqWeight + (
                propertyBinGoodTruthWeights[propertyBin] +
                        propertyBinGoodInheritanceWeights[propertyBin]
        );
        final float badBothWeights = minGqWeight + (
                propertyBinGoodTruthWeights[propertyBin] +
                        propertyBinGoodInheritanceWeights[propertyBin]
        );
        final float bothWeights = minGqWeight + (float)(propertyBinTruthWeights[propertyBin] +
                                                        propertyBinInheritanceWeights[propertyBin]);

        for(int sampleNum = 0; sampleNum < numSamples; ++sampleNum) {
            if(getSampleVariantIsTrainable(variantIndex, sampleNum)) {
                if(inheritanceTrainableSampleIndices.contains(sampleNum)) {
                    if(truthTrainableSampleIndices.contains(sampleNum)) {
                        weights[row] = sampleVariantTruth[row] ? goodBothWeights : badBothWeights;
                        //weights[row] = bothWeights;
                    } else {
                        weights[row] = sampleVariantTruth[row] ? goodInheritanceWeight : badInheritanceWeight;
                        //weights[row] = inheritanceWeight;
                    }
                } else if(truthTrainableSampleIndices.contains(sampleNum)) {
                    weights[row] = sampleVariantTruth[row] ? goodTruthWeight : badTruthWeight;
                    //weights[row] = truthWeight;
                } else {  // this sample is only trainable because we've set minGQ. Take the mean of the two weights
                    weights[row] = minGqWeight;
                }
                ++row;
            }
        }
        return row;
    }

    private float[] getDMatrixWeights(final int[] variantIndices, final boolean[] sampleVariantTruth) {
        final int numRows = sampleVariantTruth.length;
        final float[] sampleVariantWeights = new float[numRows];
        final int numSamples = getNumSamples();
        int row = 0;
        for(int variantIndex : variantIndices) {
            row = fillVariantWeights(variantIndex, numSamples, sampleVariantTruth, row, sampleVariantWeights);
        }
        if(progressVerbosity > 1) {
            System.out.format("\t\tgot %d row weights\n", row);
        }
        return sampleVariantWeights;
    }

    private float[] getDMatrixLabels(final boolean[] sampleVariantTruth) {
        final float[] labels = new float[sampleVariantTruth.length];
        for(int i = 0; i < sampleVariantTruth.length; ++i) {
            labels[i] = sampleVariantTruth[i] ? 1F : -1F;
        }
        return labels;
    }

    class DMatrixDataIterator implements Iterator<LabeledPoint> {
        private final int[] variantIndices;
        private final int numSamples;
        private final int numProperties;
        private int variantIterIndex;
        private int variantIndex;
        private int sampleIndex;
        private final float[] weights;
        private int trainableSamplesIndex;
        private final boolean[] samplePasses;
        private final float[][] sparseValuesBuffer;
        private final int[][] sparseIndicesBuffer;
        private final int[] copySparseIndicesBuffer;
        private final float[] copySparseValuesBuffer;
        private int valuesBufferIndex;

        private static final int defaultBufferLen = 1 + 32 << 10;

        DMatrixDataIterator(final int[] variantIndices) {
            this.variantIndices = variantIndices;
            this.variantIterIndex = 0;
            this.variantIndex = variantIndices.length > 0 ? variantIndices[this.variantIterIndex] : -1;
            this.sampleIndex = 0;
            this.numSamples = getNumSamples();
            this.numProperties = getNumProperties();
            this.samplePasses = new boolean[numSamples];
            this.weights = new float[numSamples];
            trainableSamplesIndex = 0;
            if(runMode == RunMode.Train && variantIndices.length > 0) {
                fillSamplePasses(variantIndex, 0, samplePasses);
                fillVariantWeights(variantIndex, numSamples, samplePasses, 0, weights);
            }

            this.copySparseIndicesBuffer = new int[numProperties];
            this.copySparseValuesBuffer = new float[numProperties];
            final int bufferLen = (int) FastMath.min(defaultBufferLen, variantIndices.length * (long) numSamples);
            this.sparseValuesBuffer = new float[bufferLen][];
            this.sparseIndicesBuffer = new int[bufferLen][];
            valuesBufferIndex = 0;
        }

        private boolean increment() {
            if(variantIterIndex >= variantIndices.length) {
                return false;
            }
            ++sampleIndex;
            if(sampleIndex == numSamples) {
                ++variantIterIndex;
                if(variantIterIndex >= variantIndices.length) {
                    return false;
                }
                variantIndex = variantIndices[variantIterIndex];
                sampleIndex = 0;
                trainableSamplesIndex = 0;
                if(runMode == RunMode.Train) {
                    fillSamplePasses(variantIndex, 0, samplePasses);
                    fillVariantWeights(variantIndex, numSamples, samplePasses, 0, weights);
                }
            }
            return true;
        }

        @Override public boolean hasNext() {
            if(variantIterIndex >= variantIndices.length) {
                return false;
            }
            if(runMode == RunMode.Train) {
                while (!getSampleVariantIsTrainable(variantIndex, sampleIndex)) {
                    if (!increment()) {
                        return false;
                    }
                }
            } else {
                while (!getSampleVariantIsFilterable(variantIndex, sampleIndex)) {
                    if (!increment()) {
                        return false;
                    }
                }
            }
            return true;
        }

        @Override
        public LabeledPoint next() {
            if(!hasNext()) {
                throw new RuntimeException("Iterated past end");
            }
            getSparseValues();
            final float[] values = sparseValuesBuffer[valuesBufferIndex];
            final int[] indices = sparseIndicesBuffer[valuesBufferIndex];
            ++valuesBufferIndex;
            if(valuesBufferIndex >= sparseValuesBuffer.length) {
                valuesBufferIndex = 0;
            }
            final float label;
            final float weight;
            final int group;
            if(runMode == RunMode.Train) {
                label = samplePasses[trainableSamplesIndex] ? 1F : -1F;
                weight = weights[trainableSamplesIndex];
                group = propertyBins[variantIndex];
            } else {
                label = 0F;
                weight = 1F;
                group = 0;
            }
            ++trainableSamplesIndex;
            final LabeledPoint labeledPoint = new LabeledPoint(label, numProperties, indices, values,
                                                               weight, group, 0);
            increment();
            return labeledPoint;
        }

        private void getSparseValues() {
            int flatIndex = 0;
            int propertyIndex = 0;
            for (final String propertyName : propertiesTable.getPropertyNames()) {
                final PropertiesTable.Property property = propertiesTable.get(propertyName);
                if (property instanceof PropertiesTable.BooleanArrProperty) {
                    final boolean value = ((PropertiesTable.BooleanArrProperty) property).values[variantIndex];
                    if (value) {
                        copySparseIndicesBuffer[flatIndex] = propertyIndex;
                        copySparseValuesBuffer[flatIndex] = 1F;
                        ++flatIndex;
                    }
                } else if (property instanceof PropertiesTable.BooleanMatProperty) {
                    final boolean value = ((PropertiesTable.BooleanMatProperty) property).values[variantIndex][sampleIndex];
                    if (value) {
                        copySparseIndicesBuffer[flatIndex] = propertyIndex;
                        copySparseValuesBuffer[flatIndex] = 1F;
                        ++flatIndex;
                    }
                } else {
                    copySparseIndicesBuffer[flatIndex] = propertyIndex;
                    copySparseValuesBuffer[flatIndex] = property.getAsFloat(variantIndex, sampleIndex);
                    ++flatIndex;
                }
                ++propertyIndex;
            }
            sparseIndicesBuffer[valuesBufferIndex] = Arrays.copyOf(copySparseIndicesBuffer, flatIndex);
            sparseValuesBuffer[valuesBufferIndex] = Arrays.copyOf(copySparseValuesBuffer, flatIndex);
        }
    }

    private DMatrix getDMatrixLabeledPoint(final int[] variantIndices, boolean cacheDMatrix) throws XGBoostError {
        if(cacheDMatrix) {
            final File cacheFile;
            try {
                cacheFile = File.createTempFile("temp", ".dmatrix.cache", tempDir);
            } catch (IOException e) {
                throw new GATKException("Error creating temp file for DMatrix cache", e);
            }
            cacheFile.deleteOnExit();

            return new DMatrix(new DMatrixDataIterator(variantIndices), cacheFile.getAbsolutePath());
        } else {
            return new DMatrix(new DMatrixDataIterator(variantIndices), null);
        }
    }

    private DMatrix getDMatrix(final int[] variantIndices) throws XGBoostError {
        if(progressVerbosity > 1) {
            System.out.println("\tgetDMatrix");
        }
        final DMatrix dMatrix;
        switch(dMatrixFillMethod) {
            case BigDenseMatrix:
                dMatrix = getDMatrixBigDenseMatrix(variantIndices);
                break;
            case RowMajorArray:
                dMatrix = getDMatrixRowMajorArray(variantIndices);
                break;
            case LabeledPointIterator:
                dMatrix = getDMatrixLabeledPoint(variantIndices, true);
                if(progressVerbosity > 2) {
                    System.out.println("\tgetDMatrix complete");
                }
                return dMatrix;  // don't try to set other properties, they're already handled
            default:
                throw new IllegalArgumentException("Unknown DMatrixFillMethod: " + dMatrixFillMethod);
        }
        final boolean[] sampleVariantTruth = getSampleVariantTruth(variantIndices);
        dMatrix.setWeight(getDMatrixWeights(variantIndices, sampleVariantTruth));
        dMatrix.setLabel(getDMatrixLabels(sampleVariantTruth));
        dMatrix.setBaseMargin(new float[sampleVariantTruth.length]);
        if(progressVerbosity > 2) {
            System.out.println("\tgetDMatrix complete");
        }
        return dMatrix;
    }

    final int[] variantIndicesForFiltering = new int[1];  // always the first variant index
    DMatrix getDMatrixForFiltering() throws XGBoostError {
        return getDMatrixLabeledPoint(variantIndicesForFiltering, false);
    }

    private Map<String, Object> getXgboostParams() {
        return new HashMap<String, Object>() {
            private static final long serialVersionUID = 0L;
            {
                put("eta", eta);
                put("max_depth", maxDepth);
                put("gamma", gamma);
                put("subsample", subsample);
                put("colsample_bytree", colsampleByTree);
                put("colsample_bylevel", colsampleByLevel);
                put("min_child_weight", minChildWeight);
                put("validate_parameters", true);
                put("objective", "binary:logistic");
                put("eval_metric", "logloss");
                if(numThreadsXgBoost > 0) {
                    put("nthread", numThreadsXgBoost);
                }
            }
        };
    }

    static class FilterLosses {
        final FilterLoss trainingLoss;
        final FilterLoss diagnosticLoss;

        FilterLosses(final FilterLoss trainingLoss, final FilterLoss diagnosticLoss) {
            this.trainingLoss = trainingLoss;
            this.diagnosticLoss = diagnosticLoss;
        }
    }

    private class DataSubset {
        final String name;
        final int[] variantIndices;
        final boolean isTrainingSet;
        final int numTrainableSampleVariants;

        final List<Float> scores;
        private int bestScoreInd;
        private float bestScore;

        final DMatrix dMatrix;
        final float[] pSampleVariantGood;
        final float[] d1Loss;
        final float[] d2Loss;

        DataSubset(final String name, final int[] variantIndices, final boolean isTrainingSet) {
            System.out.println("Building DataSubset " + name);
            this.name = name;
            this.variantIndices = variantIndices;
            this.isTrainingSet = isTrainingSet;
            numTrainableSampleVariants = getNumTrainableSampleVariants(variantIndices);
            if(progressVerbosity > 1) {
                System.out.format("\t%d variants, %d sample x variants, %d properties, %d entries\n",
                                   variantIndices.length, numTrainableSampleVariants, getNumProperties(),
                                   numTrainableSampleVariants * (long)getNumProperties());
            }

            scores = new ArrayList<>();
            bestScore = Float.POSITIVE_INFINITY;
            bestScoreInd = 0;

            try {
                this.dMatrix = getDMatrix(variantIndices);
            } catch(XGBoostError xgBoostError) {
                throw new GATKException("Error constructing DMatrix", xgBoostError);
            }
            pSampleVariantGood = new float[numTrainableSampleVariants];

            // pre-allocate derivative arrays if they are needed (i.e. this is the training set)
            d1Loss = isTrainingSet ? new float[numTrainableSampleVariants] : null;
            d2Loss = isTrainingSet ? new float[numTrainableSampleVariants] : null;

            final boolean[] sampleVariantTruth = getSampleVariantTruth(variantIndices);
            final float[] tempTruthProbs = new float[numTrainableSampleVariants];
            for(int idx = 0; idx < numTrainableSampleVariants; ++idx) {
                tempTruthProbs[idx] = sampleVariantTruth[idx] ? 1F : 0F;
            }
            System.out.println("Best possible " + name + " loss:\n\t" +
                               getLoss(tempTruthProbs, variantIndices).toString().replaceAll("\n", "\n\t"));
        }

        public int size() { return numTrainableSampleVariants; }

        private float[][] getRawPredictions(final Booster booster, final DMatrix dMatrix) {
            try {
                return booster.predict(dMatrix, true, 0);
            } catch(XGBoostError xgBoostError) {
                throw new GATKException("In " + name + " DataSubset: Predict error", xgBoostError);
            }
        }

        private void displayQuantilesFloat(final float[] arr, final String description) {
            DoubleStream.Builder builder = DoubleStream.builder();
            for (float v : arr) {
                builder.add(v);
            }
            displayPercentiles(description, builder.build());
        }

        private FilterLosses getRoundLoss(final Booster booster, final boolean training) {
            final float[][] rawPredictions = getRawPredictions(booster, dMatrix);
            for(int idx = 0; idx < rawPredictions.length; ++idx) {
                pSampleVariantGood[idx] = scaledLogitsToP(rawPredictions[idx][0]);
            }

            final FilterLoss lossForTraining;
            if(isTrainingSet && training) {
                lossForTraining = getTrainingLoss(pSampleVariantGood, d1Loss, d2Loss, variantIndices);
                System.out.format("\t%s\n\t\t%s\n", "optimizer objective",
                                  lossForTraining.toString().replaceAll("\n", "\n\t\t"));
                if(progressVerbosity > 0) {
                    System.out.format("\tBoosting round %d on %s\n", 1 + getRound(), name);
                }
                if(progressVerbosity > 1) {
                    displayQuantilesFloat(pSampleVariantGood, "pSampleVariantIsGood");
                    displayQuantilesFloat(d1Loss, "d1Loss");
                    displayQuantilesFloat(d2Loss, "d2Loss");
                }
                try {
                    booster.boost(dMatrix, d1Loss, d2Loss);
                } catch(XGBoostError xgBoostError) {
                    throw new GATKException("In " + name + " DataSubset: Boost error", xgBoostError);
                }
            } else {
                lossForTraining = getTrainingLoss(pSampleVariantGood, variantIndices);
            }

            final FilterLoss lossForDiagnostics = getLoss(pSampleVariantGood, variantIndices);


            if(progressVerbosity > 0) {
                System.out.format("\t%s\n\t\t%s\n", name, lossForDiagnostics.toString().replaceAll("\n", "\n\t\t"));
            }
            return new FilterLosses(lossForTraining, lossForDiagnostics);
        }

        public DoubleStream streamPredictions(final Booster booster) {
            return Arrays.stream(getRawPredictions(booster, dMatrix))
                .mapToDouble(rawPrediction -> scaledLogitsToP(rawPrediction[0]));
        }


        public float getLastScore() { return scores.get(scores.size() - 1); }

        public void appendEarlyStoppingScore(final FilterLosses filterLosses) {
            if(earlyStoppingLoss == EarlyStoppingLoss.Diagnostic) {
                appendScore(filterLosses.diagnosticLoss.toFloat());
            } else {
                appendScore(filterLosses.trainingLoss.toFloat());
            }
        }

        public void appendScore(final float score) {
            scores.add(score);
            if(score < bestScore) {
                bestScore = score;
                bestScoreInd = getRound();
            }
        }

        public int getRound() {
            return scores.size() - 1;
        }

        public boolean isBestScore() {
            return bestScoreInd == getRound();
        }

        public boolean stop() {
            return getRound() == maxTrainingRounds || stopEarly();
        }

        public boolean stopEarly() {
            return earlyStoppingRounds > 0 && getRound() - bestScoreInd > earlyStoppingRounds;
        }
    }

    @Override
    protected boolean needsNormalizedProperties() { return false; }

    @Override
    protected int predictBatch(final float[] outputProbabilities) {
        final float[][] rawPredictions;
        try {
            final DMatrix dMatrix = getDMatrixForFiltering();
            if(dMatrix.rowNum() == 0) {
                return 0;  // don't try to predict on empty DMatrix, XGBoost doesn't like it.
            }
            rawPredictions = booster.predict(
                    getDMatrixForFiltering(), true, 0
            );
            for(int row = 0; row < rawPredictions.length; ++ row) {
                outputProbabilities[row] = scaledLogitsToP(rawPredictions[row][0]);
            }
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error predicting", xgBoostError);
        }
        return rawPredictions.length;
    }

    @Override
    protected void loadModel(final InputStream inputStream) {
        try {
            ByteArrayOutputStream os = new ByteArrayOutputStream();
            int readValue = inputStream.read();
            while(readValue >= 0) {
                os.write(readValue);
                readValue = inputStream.read();
            }

            booster = XGBoost.loadModel(Base64.getDecoder().decode(os.toByteArray()));
        } catch(XGBoostError | IOException xgBoostError) {
            throw new GATKException("Error loading XGBoost model", xgBoostError);
        }
    }

    @Override
    protected void saveModel(final OutputStream outputStream) {
        try {
            outputStream.write(Base64.getEncoder().encode(booster.toByteArray()));
        } catch(XGBoostError | IOException xgBoostError) {
            throw new GATKException("Error saving XGBoost model", xgBoostError);
        }
    }

    private Booster initializeBooster(final List<DataSubset> dataSubsets) {
        final DataSubset trainingSubset = dataSubsets.stream()
            .filter(dataSubset -> dataSubset.isTrainingSet)
            .findFirst()
            .orElseThrow(() -> new GATKException("dataSubsets does not contain a training set"));
        final Map<String, DMatrix> watches = new HashMap<>();
        // add watches for any DataSubset that a) is not the main training set, and
        //                                     b) does not use mini batches (it can comfortably hang out in memory)
        for(final DataSubset dataSubset : dataSubsets) {
            if(!dataSubset.equals(trainingSubset)) {
                watches.put(dataSubset.name, dataSubset.dMatrix);
            }
        }

        try {
            // Do 0 rounds of training on the first mini-batch (if using, otherwise whole dMatrix) of training DataSubset
            return XGBoost.train(trainingSubset.dMatrix, getXgboostParams(), 0, watches, null, null);
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error creating Booster", xgBoostError);
        }
    }

    private boolean trainOneRound(final Booster booster, final List<DataSubset> dataSubsets) {
        // evaluate booster on all data sets, calculate derivatives on training set
        if(progressVerbosity > 0) {
            System.out.println("Training round " + (dataSubsets.get(0).getRound() + 1));
        }

        boolean needSaveCheckpoint = !dataSubsets.isEmpty();
        boolean needStop = !dataSubsets.isEmpty();
        for(final DataSubset dataSubset : dataSubsets) {
            final FilterLosses roundLosses = dataSubset.getRoundLoss(booster, true);
            dataSubset.appendEarlyStoppingScore(roundLosses);
            if (!dataSubset.isTrainingSet) {
                // check if it's time to stop iterating or save model checkpoint
                needSaveCheckpoint = needSaveCheckpoint && dataSubset.isBestScore();
                needStop = needStop && dataSubset.stop();
            }
        }

        // check if booster needs to be saved, or if early stopping is necessary
        if(needSaveCheckpoint) {
            saveModelCheckpoint();
        }
        return !needStop;
    }

    void displayFeatureImportance(final Booster booster) {
        final List<String> propertyNames = propertiesTable.getPropertyNames();
        final Map<String, Integer> featureScore;
        try {
            featureScore = booster.getFeatureScore((String) null).entrySet().stream()
                .collect(
                    Collectors.toMap(
                        entry -> propertyNames.get(Integer.parseInt(entry.getKey().substring(1))),
                        Map.Entry::getValue
                    )
                );
        } catch(XGBoostError xgBoostError) {
            throw new GATKException("Error getting feature score", xgBoostError);
        }

        System.out.println("Feature importance:");
        featureScore.entrySet().stream()
            .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
            .forEachOrdered(entry -> System.out.format("\t%s: %d\n", entry.getKey(), entry.getValue()));
    }

    @Override
    protected void trainFilter() {
        final List<DataSubset> dataSubsets = new ArrayList<DataSubset>() {
            private static final long serialVersionUID = 0L;
            {
                add(new DataSubset(TRAIN_MAT_KEY, getTrainingIndices(), true));
                add(new DataSubset(VALIDATION_MAT_KEY, getValidationIndices(), false));
            }
        };

        clearUnneededProperties();

        booster = initializeBooster(dataSubsets);
        //noinspection StatementWithEmptyBody
        while (trainOneRound(booster, dataSubsets))
            ;

        loadModelCheckpoint();
        if(progressVerbosity > 0) {
            System.out.println("Losses of final trained classifier:");
            for (final DataSubset dataSubset : dataSubsets) {
                dataSubset.getRoundLoss(booster, false);
            }
        }
        displayHistogram("Final training adjusted GQ histogram",
                          dataSubsets.get(0).streamPredictions(booster).mapToInt(this::probToScaledLogits),true);
        displayHistogram("Final training probability histogram",
                          dataSubsets.get(0).streamPredictions(booster),20, 0.0, 1.0);
        displayFeatureImportance(booster);
    }
}
