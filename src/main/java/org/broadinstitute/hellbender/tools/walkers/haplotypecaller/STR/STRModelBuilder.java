package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * STRModel builder utility.
 */
public final class STRModelBuilder extends CommandLineProgram {

    private static final double MIN_INITIAL_PI = 0.001;
    private static final double MIN_INITIAL_LAMBDA = 1.0;
    private static final double INITIAL_LEARNING_RATE = 1.0;
    private static final int MAX_ITERATIONS = 1000;
    private static final double MIN_LAMBDA = 0.001;
    private static final double MAX_LAMBDA = 2.0; // beyond this lambda and rcd of 2 would be more likelily than 1 which is non-sensical.
    private static final double MIN_PI = 0.00001;
    private static final double MAX_PI = 0.99999;
    private static final double LOG10_EPSILON = 0.00001;

    /**
     * 95% max average repeat count exepcted.
     * <p>
     * sapply(1:100, function(n) { quantile(sapply(1:1000, function(x) {mean(rtpois(n,2)) }), 0.95) })
     */
    private static final double[] MAX_AVG_RC = new double[]{0,
            5.000000, 4.000000, 3.666667, 3.262500, 3.210000, 3.166667, 3.142857, 3.125000,
            3.111111, 3.000000, 2.909091, 3.000000, 2.923077, 2.928571, 2.866667, 2.875000,
            2.823529, 2.833333, 2.789474, 2.750000, 2.809524, 2.772727, 2.739130, 2.750000,
            2.720000, 2.730769, 2.740741, 2.750000, 2.689655, 2.700000, 2.709677, 2.687500,
            2.696970, 2.676471, 2.658571, 2.666667, 2.702703, 2.657895, 2.666667, 2.650000,
            2.658537, 2.642857, 2.674419, 2.613636, 2.622222, 2.608696, 2.618085, 2.604167,
            2.632653, 2.620000, 2.607843, 2.596154, 2.603774, 2.611111, 2.600000, 2.589286,
            2.614035, 2.586207, 2.576271, 2.583333, 2.573770, 2.580645, 2.587302, 2.578125,
            2.553846, 2.560606, 2.567164, 2.558824, 2.551449, 2.557857, 2.549296, 2.569444,
            2.575342, 2.567568, 2.560000, 2.552632, 2.559091, 2.564103, 2.556962, 2.538125,
            2.530864, 2.548780, 2.542169, 2.559524, 2.529412, 2.534884, 2.540230, 2.534091,
            2.528090, 2.522778, 2.527473, 2.532609, 2.537634, 2.532447, 2.526316, 2.521354,
            2.536082, 2.520918, 2.525253, 2.530500};

    private static final org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();


    @Argument(doc = "the reference fasta file", fullName = "reference", shortName = "R")
    protected File referenceFile;

    @Argument(doc = "the input file", fullName = "input", shortName = "I")
    protected File inputFile;

    @Argument(doc = "the output file", fullName = "output", shortName = "O")
    protected File outputFile;

    @Argument(doc = "count output directory", fullName = "countOutputDirectory", shortName = "countDir", optional = true)
    protected File countDir;

    private List<ImmutablePair<STRContext, Map<String, List<STRAllele>>>> trainingData;

    public STRModelBuilder() {
        this(10000);
    }

    public STRModelBuilder(final int expectedTrainingDataSize) {
        trainingData = new ArrayList<>(expectedTrainingDataSize);

    }

    public void addTrainingData(final STRContext context, final Map<String, List<STRAllele>> truth) {
        if (validLikelihoods(context)) {
            trainingData.add(new ImmutablePair<>(context, truth));
        } else {
            logger.warn("STR Context contain some invalid likelihoods will be ignored: " + context.toLongString());
        }
    }

    private boolean validLikelihoods(final STRContext context) {
        final ReadLikelihoods<STRAllele> likelihoods = context.getLikelihoods();
        final int alleleCount = likelihoods.numberOfAlleles();
        for (final String sample : likelihoods.samples()) {
            final int sampleIndex = likelihoods.indexOfSample(sample);
            final LikelihoodMatrix<STRAllele> sampleMatrix = likelihoods.sampleMatrix(sampleIndex);
            final int readCount = sampleMatrix.numberOfReads();
            for (int i = 0; i < readCount; i++) {
                for (int j = 0; j < alleleCount; j++) {
                    if (!Double.isFinite(sampleMatrix.get(j, i))) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    public STRModelParameterCollection fitTrainingData(final STRModelCovariate.List covariates) {
        final int setCount = covariates.stream().mapToInt(STRModelCovariate::size).reduce(1, (a, b) -> a * b);
        final Map<Integer, List<ImmutablePair<STRContext, Map<String, List<STRAllele>>>>> trainingDataBySetIndex =
                trainingData.stream().collect(Collectors.groupingBy(pair -> covariates.indexFor(pair.getLeft())));
        final List<STRModelParameters> sets = new ArrayList<>(setCount);
        logger.info("Training elements: " + setCount);
        for (int i = 0; i < setCount; i++) {
            final List<ImmutablePair<STRContext, Map<String, List<STRAllele>>>> setData = trainingDataBySetIndex.getOrDefault(i, Collections.emptyList());
            if (setData.isEmpty()) {
                sets.add(new STRModelParameters(0, Double.NaN, Double.NaN));
            } else {
                sets.add(fitTrainingData(covariates, setData));
            }
            logger.info("Trained: " + i + " with results " + sets.get(sets.size() - 1));
            if (countDir != null) {
                writeCounts(covariates, i, setData);
            }
        }
        return new STRModelParameterCollection(covariates, sets);
    }

    private void writeCounts(final STRModelCovariate.List covariates, final int index, final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData) {
        countDir.mkdir();
        java.util.List<Object> covariateValues = covariates.valuesAt(index);
        final File fileName = new File(countDir, "counts_" + covariateValues.stream().map(Object::toString).collect(Collectors.joining("_")) + ".list");
        try (final PrintWriter writer = new PrintWriter(new FileWriter(fileName))) {
            for (final STRModelCovariate covariate : covariates) {
                writer.println("#Covariate " + covariate.getClass().getSimpleName() + " " + covariates.valueAt(index, covariate));
            }
            writer.println("# First column is the count for the called allele, ");
            writer.println("# second is the combined count for alleles with one repeat unit difference, ");
            writer.println("# third for 2 repeat units, ");
            writer.println("# and so forth");
            for (final Pair<STRContext, Map<String, List<STRAllele>>> pair : setData) {
                int maxNonTruthCount = 0;
                final STRContext context = pair.getLeft();
                final int truth = context.getAlleles().indexOf(pair.getRight().values().stream().findFirst().get().get(0));
                final int truthRepeatCount = context.getAlleles().get(truth).repeatCount;
                final int maximumDifference = Math.max(context.getAlleles().getMaximumRepeatCount() - truthRepeatCount,
                        truthRepeatCount - context.getAlleles().getMinimumRepeatCount());
                final int[] counts = new int[maximumDifference + 1];
                final int[] alleleDepths = context.getAlleleDepths();
                long sum = 0;
                int total = 0;
                for (int i = 0; i < context.getAlleles().size(); i++) {
                    final int diff = Math.abs(truthRepeatCount - context.getAlleles().get(i).repeatCount);
                    counts[diff] += alleleDepths[i];
                    if (i != truth && alleleDepths[i] > maxNonTruthCount) {
                        maxNonTruthCount = alleleDepths[i];
                        sum += alleleDepths[i] * diff;
                        total += alleleDepths[i];
                    }
                }
                if (maxNonTruthCount > alleleDepths[truth]) {
                    logger.debug("Skipping counts due to infrequent truth: " + Utils.join("\t", context.getRepeatLocus(), Utils.join("\t", counts)));
                    continue;
                }
                if (sum / (double) total > MAX_AVG_RC[total >= MAX_AVG_RC.length ? MAX_AVG_RC.length - 1 : total]) {
                    logger.debug("Skipping counts due to extreme average non-truth read-count diff : " + Utils.join("\t", context.getRepeatLocus(), Utils.join("\t", counts)));
                    continue;
                }
                writer.println(Utils.join("\t", context.getRepeatLocus(), Utils.join("\t", counts)));
            }
            writer.close();

        } catch (final IOException ex) {
            throw new GATKException(ex.getMessage(), ex);
        }
    }

    private STRModelParameters fitTrainingData(final STRModelCovariate.List covariates, final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData) {
        STRModelParameters result = new STRModelParameters(setData.size(), estimateErrorRate(setData), initialLambda(setData));

        double likelihood = calculateLog10Likelihood(setData, result);
        double learningRate = INITIAL_LEARNING_RATE;
        //boolean hasEverFailedToImprove = false;
        int iteration;
        outter:
        for (iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
            final double errorRateDerivate = calculateErrorRateDerivative(setData, result);
            final double lambdaDerivate = calculateLambdaDerivative(setData, result);
            while (iteration < MAX_ITERATIONS) {
                final double updatedErrorRate = Math.min(Math.max(result.pi + learningRate * errorRateDerivate, MIN_PI), MAX_PI);
                final double updatedLambda = Math.min(Math.max(result.lambda + learningRate * lambdaDerivate, MIN_LAMBDA), MAX_LAMBDA);
                final STRModelParameters updatedResult = new STRModelParameters(result.size, updatedErrorRate, updatedLambda);
                final double updatedLikelihood = calculateLog10Likelihood(setData, updatedResult);
                if (updatedLikelihood >= likelihood && Double.isFinite(updatedLikelihood)) {
                    learningRate *= 1.5;
                    result = updatedResult;
                    if (updatedLikelihood - likelihood < LOG10_EPSILON) {
                        likelihood = updatedLikelihood;
                        break outter;
                    }
                    likelihood = updatedLikelihood;
                    break;
                } else {
                    //           hasEverFailedToImprove = true;
                    learningRate *= .5;
                    iteration++;
                }
            }
        }
        result.iterations = iteration;
        result.likelihood = likelihood;
        result.modelFreeLikelihood = calculateLog10Likelihood(setData, new STRModelParameters(0, 0, 1));
        if (logger.isDebugEnabled()) {
            logger.debug("Trained " + result.toString() + " with the following data: ");
            for (final Pair<STRContext, Map<String, List<STRAllele>>> datum : setData) {
                logger.debug("    " + datum.getLeft());
            }
            logger.debug("Trained " + result.toString() + " end");
        }
        return result;
    }

    private double calculateErrorRateDerivative(final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData, final STRModelParameters result) {
        return calculateDerivative(setData, result::piDerivativeMatrix, result::log10TransformationMatrix);
    }

    private double calculateLambdaDerivative(final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData, final STRModelParameters result) {
        return calculateDerivative(setData, result::lambdaDerivativeMatrix, result::log10TransformationMatrix);
    }

    private double calculateLog10Likelihood(final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData,
                                            final STRModelParameters paramSet) {
        final double[] dataLikelihoods = new double[setData.size()];
        int nextDataIndex = 0;
        for (final Pair<STRContext, Map<String, List<STRAllele>>> setDatum : setData) {
            final ReadLikelihoods<STRAllele> likelihoods = setDatum.getLeft().getLikelihoods();
            final RealMatrix log10ProductMatrix = paramSet.log10TransformationMatrix(likelihoods.alleles());
            //final ReadLikelihoods<STRAllele> likelihoods = originalLikelihoods.transform(log10TransformationMatrix);
            final int likelhoodAlleleCount = likelihoods.numberOfAlleles();
            final double[] likelihoodAlleleContributions = new double[likelhoodAlleleCount];
            final double[] readLikelihoods = new double[likelihoods.readCount()];
            int nextReadIndex = 0;
            for (final String sample : likelihoods.samples()) {
                final List<STRAllele> truth = setDatum.getRight().get(sample);
                final int truthAlleleCount = truth.size();
                final double[] truthAlleleLikelihoods = new double[truthAlleleCount];
                final LikelihoodMatrix<STRAllele> sampleLikelihoods = likelihoods.sampleMatrix(likelihoods.indexOfSample(sample));
                final int readCount = sampleLikelihoods.numberOfReads();
                for (int readIndex = 0; readIndex < readCount; readIndex++) {
                    for (int truthAlleleIndex = 0; truthAlleleIndex < truthAlleleCount; truthAlleleIndex++) {
                        final int truthAlleleLikelihoodIndex = likelihoods.indexOfAllele(truth.get(truthAlleleIndex));
                        for (int likelihoodAlleleIndex = 0; likelihoodAlleleIndex < likelhoodAlleleCount; likelihoodAlleleIndex++) {
                            likelihoodAlleleContributions[likelihoodAlleleIndex] =
                                    log10ProductMatrix.getEntry(truthAlleleLikelihoodIndex, likelihoodAlleleIndex) + sampleLikelihoods.get(likelihoodAlleleIndex, readIndex);
                        }
                        final double log10TruthAlleleReadLikelihood = MathUtils.log10sumLog10(likelihoodAlleleContributions) - MathUtils.log10(truthAlleleCount);
                        truthAlleleLikelihoods[truthAlleleIndex] = log10TruthAlleleReadLikelihood;
                    }
                    readLikelihoods[nextReadIndex++] = MathUtils.log10sumLog10(truthAlleleLikelihoods);
                }
            }
            dataLikelihoods[nextDataIndex++] = MathUtils.sum(readLikelihoods);
        }
        return MathUtils.sum(dataLikelihoods);
    }

    private double calculateDerivative(final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData,
                                       final Function<List<STRAllele>, RealMatrix> derivativeFactorsProducer,
                                       final Function<List<STRAllele>, RealMatrix> factorsProducer) {
        final double[] dataDerivative = new double[setData.size()];
        int nextDataIndex = 0;
        for (final Pair<STRContext, Map<String, List<STRAllele>>> setDatum : setData) {
            final ReadLikelihoods<STRAllele> likelihoods = setDatum.getLeft().getLikelihoods();
            final RealMatrix derivativeMatrix = derivativeFactorsProducer.apply(likelihoods.alleles());
            final RealMatrix log10FactorsMatrix = factorsProducer.apply(likelihoods.alleles());
            final int likelhoodAlleleCount = likelihoods.numberOfAlleles();
            final double[] likelihoodAlleleContributions = new double[likelhoodAlleleCount];
            final double[] derivativeAlleleContributions = new double[likelhoodAlleleCount];
            final double[] readDerivative = new double[likelihoods.readCount()];
            final double[] readLikelihood = new double[likelihoods.readCount()];
            int nextReadIndex = 0;
            for (final String sample : likelihoods.samples()) {
                final LikelihoodMatrix<STRAllele> sampleLikelihoods = likelihoods.sampleMatrix(likelihoods.indexOfSample(sample));
                final List<STRAllele> truth = setDatum.getRight().get(sample);
                final int truthAlleleCount = truth.size();
                final double[] truthAlleleDerivative = new double[truthAlleleCount];
                final double[] truthAlleleLikelihood = new double[truthAlleleCount];
                final int readCount = sampleLikelihoods.numberOfReads();
                for (int readIndex = 0; readIndex < readCount; readIndex++) {
                    for (int truthAlleleIndex = 0; truthAlleleIndex < truthAlleleCount; truthAlleleIndex++) {
                        final int truthAlleleLikelihoodIndex = likelihoods.indexOfAllele(truth.get(truthAlleleIndex));
                        for (int likelihoodAlleleIndex = 0; likelihoodAlleleIndex < likelhoodAlleleCount; likelihoodAlleleIndex++) {
                            likelihoodAlleleContributions[likelihoodAlleleIndex] = log10FactorsMatrix.getEntry(truthAlleleLikelihoodIndex, likelihoodAlleleIndex)
                                    + sampleLikelihoods.get(likelihoodAlleleIndex, readIndex);
                            derivativeAlleleContributions[likelihoodAlleleIndex] = Math.pow(10, sampleLikelihoods.get(likelihoodAlleleIndex, readIndex))
                                    * derivativeMatrix.getEntry(truthAlleleLikelihoodIndex, likelihoodAlleleIndex);
                        }
                        truthAlleleLikelihood[truthAlleleIndex] = MathUtils.log10sumLog10(likelihoodAlleleContributions) - Math.log10(truthAlleleCount);
                        truthAlleleDerivative[truthAlleleIndex] = MathUtils.sum(derivativeAlleleContributions) / truthAlleleCount;
                    }
                    readLikelihood[nextReadIndex] = MathUtils.log10sumLog10(truthAlleleLikelihood);
                    readDerivative[nextReadIndex] = MathUtils.sum(truthAlleleDerivative) / Math.pow(10, readLikelihood[nextReadIndex]);
                    nextReadIndex++;
                }
            }
            dataDerivative[nextDataIndex] = MathUtils.sum(readDerivative);
            nextDataIndex++;
        }
        return MathUtils.sum(dataDerivative);
    }

    private double initialLambda(final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData) {
        double totalSum = 1.0;
        double totalWeight = 1.0;
        for (final Pair<STRContext, Map<String, List<STRAllele>>> setDatum : setData) {
            final ReadLikelihoods<STRAllele> likelihoods = setDatum.getLeft().getLikelihoods();
            for (final String sample : setDatum.getRight().keySet()) {
                final List<STRAllele> truth = setDatum.getRight().get(sample);
                final OptionalDouble averageDistance = likelihoods.bestAllelesBreakingTies().stream()
                        .filter(ba -> ba.sample.equals(sample))
                        .filter(ba -> !truth.contains(ba.allele))
                        .flatMap(ba -> truth.stream().map(l -> Math.abs(l.repeatCount - ba.allele.repeatCount)))
                        .mapToInt(i -> i)
                        .average();
                if (averageDistance.isPresent()) {
                    final double weight = Math.log(likelihoods.sampleMatrix(likelihoods.indexOfSample(sample)).numberOfReads()); // avoid giving to much weight to sites with many reads.
                    totalWeight += weight;
                    totalSum += averageDistance.getAsDouble() * weight;
                }
            }
        }
        return Math.max(MIN_INITIAL_LAMBDA, totalSum / totalWeight);
    }

    private double estimateErrorRate(final List<? extends Pair<STRContext, Map<String, List<STRAllele>>>> setData) {


        return MIN_INITIAL_PI;
  //      double totalSum = 0.0;
  //      double totalWeight = 0.0;
  //      for (final Pair<STRContext, Map<String, List<STRAllele>>> setDatum : setData) {
  //          final ReadLikelihoods<STRAllele> likelihoods = setDatum.getRight().getLikelihoods();
  //          for (final String sample : setDatum.getRight().keySet()) {
    //            final List<STRAllele> truth = setDatum.getRight().get(sample);
 //               final double weight = Math.log(likelihoods.sampleMatrix(likelihoods.indexOfSample(sample)).numberOfReads()); // avoid giving to much weight to sites with many reads.
 //               totalWeight += weight;
  //              final long mismatches = likelihoods.bestAllelesBreakingTies().stream()
  //                      .filter(ba -> ba.sample.equals(sample))
 //                       .filter(ba -> !truth.contains(ba.allele)).count();
 //               totalSum += (mismatches / likelihoods.readCount()) * weight;
 //           }
 //       }
///        return Math.max(MIN_INITIAL_PI, totalSum / totalWeight);
    }

    private GATKRead readFactory(final String name) {
        final String[] parts = name.split("\\s+");
        if (parts.length == 0) {
            throw new GATKException("cannot handle empty read-name");
        } else {
            final SAMRecord read = new SAMRecord(ArtificialReadUtils.createArtificialSamHeader());
            read.setReadName(parts[0]);
            if (parts.length > 1 && parts[1].matches("\\d+\\/\\d+")) {
                final int numerator = Integer.parseInt(parts[1].substring(0, parts[1].indexOf('/')));
                final int denominator = Integer.parseInt(parts[1].substring(parts[1].indexOf('/') + 1));
                if (denominator > 1) {
                    read.setReadPairedFlag(true);
                    read.setFirstOfPairFlag(numerator <= 1);
                }
            }
            return new SAMRecordToGATKReadAdapter(read);
        }
    }

    private static final Pattern SAMPLE_NAME_PATTERN = Pattern.compile(".*GT=\\[\\[(\\S+).*");

    @Override
    protected String doWork() {

        int maxMeanRepeatCount = 0;
        int maxRepeatUnitLength = 0;
        int minMeanRepeatCount = 100;


        try (final BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
            final ReferenceSequenceFile reference = new CachingIndexedFastaSequenceFile(referenceFile.toPath());
            String line;
            logger.info("Reading training data...");
            while ((line = reader.readLine()) != null) {
                final Matcher sampleNameMatcher = SAMPLE_NAME_PATTERN.matcher(line);
                if (!sampleNameMatcher.matches()) {
                    throw new GATKException("cannot find sample name in input line: " + line);
                }
                final String sampleName = sampleNameMatcher.group(1);
                final String[] parts = line.split("\t");

                final STRContext strContext = STRContext.fromLongString(parts[5], new GenomeLocParser(reference.getSequenceDictionary()), this::readFactory);
//                final int[] alleleDepths = strContext.getAlleleDepths();
//                if (MathUtils.maxElementIndex(alleleDepths) != Integer.parseInt(parts[3])) {
//                    continue;
//                }
                final STRAlleleSet alleles = strContext.getAlleles();
                final int callAlleleRefLengthDiff = Integer.parseInt(parts[4]);
                final STRAllele callAllele = alleles.getByRepeatCount(
                        alleles.getReferenceRepeatCount() + callAlleleRefLengthDiff / alleles.getRepeatUnitLength());
                addTrainingData(strContext, Collections.singletonMap(sampleName, Collections.nCopies(2, callAllele)));
                if (alleles.getRepeatUnitLength() > maxRepeatUnitLength) {
                    maxRepeatUnitLength = alleles.getRepeatUnitLength();
                }
                final double meanRepeatUnitCount = alleles.stream().mapToInt(STRAllele::getRepeatCount).average().orElseThrow(() -> new GATKException("invalid str-allele set " + alleles));
                if (meanRepeatUnitCount > maxMeanRepeatCount) {
                    maxMeanRepeatCount = (int) Math.ceil(meanRepeatUnitCount);
                }
                if (meanRepeatUnitCount < minMeanRepeatCount) {
                    minMeanRepeatCount = (int) Math.floor(meanRepeatUnitCount);
                }
                if (trainingData.size() % 100 == 0) {
                    logger.info("Loaded " + trainingData.size() + " sites");
                }
            }
            logger.info("Training data size: " + trainingData.size());
            final STRModelParameterCollection parameters = fitTrainingData(new STRModelCovariate.List(Arrays.asList(new STRModelCovariate.STRUnitLength(10), new STRModelCovariate.STRMeanRepeatCount(minMeanRepeatCount, maxMeanRepeatCount))));
            parameters.write(outputFile);
        } catch (final IOException ex) {
            throw new GATKException(ex.getMessage(), ex);
        }
        return "SUCCESS";
    }
}
