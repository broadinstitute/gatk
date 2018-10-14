package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.denoising.HDF5SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.denoising.SVDDenoisedCopyRatioResult;
import org.broadinstitute.hellbender.tools.copynumber.denoising.SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Integration test for {@link CreateReadCountPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CreateReadCountPanelOfNormalsIntegrationTest extends CommandLineProgramTest {
    private static final int RANDOM_SEED = 1;

    private static final int NUM_GOOD_SAMPLES = 95;
    private static final int NUM_BAD_SAMPLES_WITH_TOO_MANY_ZEROS = 5;
    private static final int NUM_SAMPLES = NUM_GOOD_SAMPLES + NUM_BAD_SAMPLES_WITH_TOO_MANY_ZEROS;

    private static final int NUM_GOOD_INTERVALS = 95;
    private static final int NUM_BAD_INTERVALS_WITH_TOO_MANY_ZEROS = 5;
    private static final int NUM_INTERVALS = NUM_GOOD_INTERVALS + NUM_BAD_INTERVALS_WITH_TOO_MANY_ZEROS;

    private static final double DEPTH_SCALE_FACTOR = 10000.;
    private static final double ALPHA_MIN = 0.01;
    private static final double ALPHA_MAX = 100.;
    private static final double MEAN_BIAS_SHAPE = 10.;
    private static final double MEAN_BIAS_SCALE = 0.1;

    //we test only for filtering of samples and intervals with too many zeros
    private static final double MINIMUM_INTERVAL_MEDIAN_PERCENTILE = 0.;
    private static final double MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE = 5.;
    private static final double MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE = 5.;
    private static final double EXTREME_SAMPLE_MEDIAN_PERCENTILE = 0.;

    //test that number of eigenvalues is recovered for a few different values using fraction of variance as a heuristic
    private static final int NUMBER_OF_EIGENVALUES_REQUESTED = NUM_SAMPLES;
    private static final List<Integer> TRUE_NUMBER_OF_EIGENVALUES_LIST = Arrays.asList(1, 4);
    private static final double FRACTION_OF_VARIANCE_EXPLAINED_THRESHOLD = 0.95;        //only a rough threshold---generating different test data may cause failures

    //test that denoised log2 copy ratios are sufficiently denoised
    private static final double DENOISED_LOG2CR_STANDARD_DEVIATION_THRESHOLD = 0.15;    //generating different test data may cause failures

    //a reasonable default GC bias curve (borrowed from GCBiasCorrectorUnitTest)
    private static final Function<Double, Double> QUADRATIC_GC_BIAS_CURVE = gc -> 0.5 + 2 * gc * (1 - gc);

    private static final SAMSequenceDictionary SEQUENCE_DICTIONARY = new SAMSequenceDictionary(Collections.singletonList(
            new SAMSequenceRecord("1", NUM_INTERVALS)));

    @DataProvider(name = "dataPanelOfNormals")
    public Object[][] dataPanelOfNormals() {
        final RandomDataGenerator rng = new RandomDataGenerator();
        rng.reSeed(RANDOM_SEED);

        //make fake intervals
        final List<SimpleInterval> intervals = IntStream.range(1, NUM_INTERVALS + 1)
                .mapToObj(i -> new SimpleInterval("1", i, i))
                .collect(Collectors.toList());

        //randomly generate GC content in each interval (borrowed from GCBiasCorrectorUnitTest)
        // gc_i ~ Bound(Normal(0.5, 0.2), 0.05, 0.95)
        final double[] intervalGCContent = IntStream.range(0, NUM_INTERVALS)
                .mapToDouble(n -> rng.nextGaussian(0.5, 0.2))
                .map(x -> Math.min(x, 0.95)).map(x -> Math.max(x, 0.05))
                .toArray();
        //assume all samples have the same GC bias (this does not add any principal components)
        final double[] intervalGCBias = Arrays.stream(intervalGCContent)
                .map(QUADRATIC_GC_BIAS_CURVE::apply)
                .toArray();
        //write GC-content annotations to file
        final AnnotatedIntervalCollection annotatedIntervals = new AnnotatedIntervalCollection(
                new SimpleLocatableMetadata(SEQUENCE_DICTIONARY),
                IntStream.range(0, NUM_INTERVALS)
                        .mapToObj(i -> new AnnotatedInterval(intervals.get(i),
                                new AnnotationMap(Collections.singletonList(Pair.of(CopyNumberAnnotations.GC_CONTENT, intervalGCContent[i])))))
                        .collect(Collectors.toList()));
        final File annotatedIntervalsFile = createTempFile("annotated-intervals", ".tsv");
        annotatedIntervals.write(annotatedIntervalsFile);

        final List<List<Object>> data = new ArrayList<>();
        for (final Integer trueNumberOfEigenvalues : TRUE_NUMBER_OF_EIGENVALUES_LIST) {

            //simulate data from simple version of gCNV coverage model
            // alpha_j ~ Uniform(alpha_min, alpha_max)
            // sigma_j = 1 / sqrt(alpha_j)
            // W_ji ~ Normal(0, sigma_j)
            // z_sj ~ Normal(0, 1)
            // m_i ~ Gamma(m_shape, m_scale)
            // n_si ~ Poisson(d * gc_bias_i * exp(z_sj W_ji + m_i))

            final List<Double> sigmas = IntStream.range(0, trueNumberOfEigenvalues)
                    .mapToObj(i -> rng.nextUniform(ALPHA_MIN, ALPHA_MAX))
                    .map(alpha -> 1. / Math.sqrt(alpha))
                    .collect(Collectors.toList());

            final RealMatrix W = new Array2DRowRealMatrix(trueNumberOfEigenvalues, NUM_INTERVALS);
            W.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int row, int column, double value) {
                    return rng.nextGaussian(0., sigmas.get(row));
                }
            });

            final RealMatrix z = new Array2DRowRealMatrix(NUM_SAMPLES, trueNumberOfEigenvalues);
            z.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int row, int column, double value) {
                    return rng.nextGaussian(0., 1.);
                }
            });

            final RealVector m = new ArrayRealVector(IntStream.range(0, NUM_INTERVALS)
                    .mapToDouble(i -> rng.nextGamma(MEAN_BIAS_SHAPE, MEAN_BIAS_SCALE))
                    .toArray());

            final RealMatrix zdotW = z.multiply(W);
            final RealMatrix bias = new Array2DRowRealMatrix(NUM_SAMPLES, NUM_INTERVALS);
            bias.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int row, int column, double value) {
                    return Math.exp(zdotW.getEntry(row, column) + m.getEntry(column));
                }
            });

            final RealMatrix counts = new Array2DRowRealMatrix(NUM_SAMPLES, NUM_INTERVALS);
            counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int row, int column, double value) {
                    return rng.nextPoisson(DEPTH_SCALE_FACTOR * intervalGCBias[column] * bias.getEntry(row, column));
                }
            });

            //corrupt first NUM_BAD_SAMPLES_WITH_TOO_MANY_ZEROS samples by randomly adding zeros
            //to 5 * MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE / 100. of intervals
            for (int sampleIndex = 0; sampleIndex < NUM_BAD_SAMPLES_WITH_TOO_MANY_ZEROS; sampleIndex++) {
                for (int intervalIndex = 0; intervalIndex < NUM_INTERVALS; intervalIndex++) {
                    if (rng.nextUniform(0., 1.) < 5 * MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE / 100.) {
                        counts.setEntry(sampleIndex, intervalIndex, 0.);
                    }
                }
            }

            //corrupt first NUM_BAD_INTERVALS_WITH_TOO_MANY_ZEROS intervals by randomly adding zeros
            //to 5 * MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE / 100. of samples
            for (int intervalIndex = 0; intervalIndex < NUM_BAD_INTERVALS_WITH_TOO_MANY_ZEROS; intervalIndex++) {
                for (int sampleIndex = 0; sampleIndex < NUM_SAMPLES; sampleIndex++) {
                    if (rng.nextUniform(0., 1.) < 5 * MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE / 100.) {
                        counts.setEntry(sampleIndex, intervalIndex, 0.);
                    }
                }
            }

            //make input files from counts matrix
            final List<File> inputTSVFiles = new ArrayList<>(NUM_SAMPLES);
            final List<File> inputHDF5Files = new ArrayList<>(NUM_SAMPLES);
            for (int sampleIndex = 0; sampleIndex < NUM_SAMPLES; sampleIndex++) {
                final File inputTSVFile = createTempFile("sample-" + sampleIndex, ".tsv");
                final File inputHDF5File = createTempFile("sample-" + sampleIndex, ".hdf5");
                final double[] sampleCounts = counts.getRow(sampleIndex);
                final SimpleCountCollection scc = new SimpleCountCollection(
                        new SimpleSampleLocatableMetadata(
                                "sample_" + sampleIndex,
                                SEQUENCE_DICTIONARY),
                        IntStream.range(0, NUM_INTERVALS)
                                .mapToObj(i -> new SimpleCount(intervals.get(i), (int) sampleCounts[i]))
                                .collect(Collectors.toList()));
                scc.write(inputTSVFile);
                inputTSVFiles.add(inputTSVFile);
                scc.writeHDF5(inputHDF5File);
                inputHDF5Files.add(inputHDF5File);
            }

            for (final File inputAnnotatedIntervalsFile : Arrays.asList(annotatedIntervalsFile, null)) {
                //counts for all samples as TSV files
                data.add(Arrays.asList(
                        inputTSVFiles,
                        inputAnnotatedIntervalsFile,
                        trueNumberOfEigenvalues));

                //counts for all samples as HDF5 files
                data.add(Arrays.asList(
                        inputHDF5Files,
                        inputAnnotatedIntervalsFile,
                        trueNumberOfEigenvalues));

                //mix of TSV and HDF5 files
                data.add(Arrays.asList(
                        ListUtils.union(
                                inputTSVFiles.subList(0, NUM_SAMPLES / 2),
                                inputHDF5Files.subList(NUM_SAMPLES / 2, NUM_SAMPLES)),
                        inputAnnotatedIntervalsFile,
                        trueNumberOfEigenvalues));
            }
        }
        return data.stream().map(List::toArray).toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataPanelOfNormals")
    public void test(final List<File> inputFiles,
                     final File annotatedIntervalsFile,
                     final int expectedNumberOfEigenvalues) {
        final File resultOutputFile = createTempFile("create-read-count-panel-of-normals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME, Double.toString(MINIMUM_INTERVAL_MEDIAN_PERCENTILE))
                .addArgument(CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME, Double.toString(MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE))
                .addArgument(CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME, Double.toString(MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE))
                .addArgument(CreateReadCountPanelOfNormals.EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME, Double.toString(EXTREME_SAMPLE_MEDIAN_PERCENTILE))
                .addArgument(CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME, Integer.toString(NUMBER_OF_EIGENVALUES_REQUESTED))
                .addOutput(resultOutputFile);
        if (annotatedIntervalsFile != null) {
            argsBuilder.addFileArgument(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME, annotatedIntervalsFile);
        }
        inputFiles.forEach(argsBuilder::addInput);
        runCommandLine(argsBuilder);
        testPanelOfNormals(annotatedIntervalsFile, expectedNumberOfEigenvalues, resultOutputFile);
    }

    @Test(dataProvider = "dataPanelOfNormals")
    public void testSingleSample(final List<File> inputFiles,
                                 final File annotatedIntervalsFile,
                                 final int expectedNumberOfEigenvalues) {   //ignored in this test
        final File resultOutputFile = createTempFile("create-read-count-panel-of-normals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME, Double.toString(MINIMUM_INTERVAL_MEDIAN_PERCENTILE))
                .addArgument(CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME, Double.toString(MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE))
                .addArgument(CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME, Double.toString(MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE))
                .addArgument(CreateReadCountPanelOfNormals.EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME, Double.toString(EXTREME_SAMPLE_MEDIAN_PERCENTILE))
                .addArgument(CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME, Integer.toString(NUMBER_OF_EIGENVALUES_REQUESTED))
                .addOutput(resultOutputFile);
        if (annotatedIntervalsFile != null) {
            argsBuilder.addFileArgument(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME, annotatedIntervalsFile);
        }
        argsBuilder.addInput(inputFiles.get(NUM_BAD_SAMPLES_WITH_TOO_MANY_ZEROS));  //use only first good sample
        runCommandLine(argsBuilder);                                                //just check that we can build the panel; no other assertions checked
    }

    @Test(dataProvider = "dataPanelOfNormals")
    public void testZeroEigensamples(final List<File> inputFiles,
                                     final File annotatedIntervalsFile,
                                     final int expectedNumberOfEigenvalues) {   //ignored in this test
        final File resultOutputFile = createTempFile("create-read-count-panel-of-normals-test", ".tsv");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME, Double.toString(MINIMUM_INTERVAL_MEDIAN_PERCENTILE))
                .addArgument(CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME, Double.toString(MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE))
                .addArgument(CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME, Double.toString(MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE))
                .addArgument(CreateReadCountPanelOfNormals.EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME, Double.toString(EXTREME_SAMPLE_MEDIAN_PERCENTILE))
                .addArgument(CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME, "0")
                .addOutput(resultOutputFile);
        if (annotatedIntervalsFile != null) {
            argsBuilder.addFileArgument(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME, annotatedIntervalsFile);
        }
        inputFiles.forEach(argsBuilder::addInput);
        runCommandLine(argsBuilder);
        try (final HDF5File hdf5PanelOfNormalsFile = new HDF5File(resultOutputFile)) {
            final SVDReadCountPanelOfNormals panelOfNormals = HDF5SVDReadCountPanelOfNormals.read(hdf5PanelOfNormalsFile);

            Assert.assertEquals(panelOfNormals.getNumEigensamples(), 0);

            //denoise last sample (which is not a bad sample) in original counts using zero eigensamples
            final RealMatrix counts = new Array2DRowRealMatrix(panelOfNormals.getOriginalReadCounts());
            final List<SimpleInterval> originalIntervals = panelOfNormals.getOriginalIntervals();
            final SimpleCountCollection sampleCounts = new SimpleCountCollection(
                    new SimpleSampleLocatableMetadata("test-sample", SEQUENCE_DICTIONARY),
                    IntStream.range(0, NUM_INTERVALS)
                            .mapToObj(i -> new SimpleCount(originalIntervals.get(i), (int) counts.getEntry(counts.getRowDimension() - 1, i)))
                            .collect(Collectors.toList()));
            final SVDDenoisedCopyRatioResult denoisedResult = panelOfNormals.denoise(sampleCounts, 0);
            //check that the denoised copy ratios are identical to the standardized copy ratios
            Assert.assertEquals(denoisedResult.getDenoisedCopyRatios(), denoisedResult.getStandardizedCopyRatios());

            //check that exceptions are thrown when attempting to get singular values and eigenvectors
            Assert.assertThrows(UnsupportedOperationException.class, panelOfNormals::getSingularValues);
            Assert.assertThrows(UnsupportedOperationException.class, panelOfNormals::getEigensampleVectors);
        }
    }

    private void testPanelOfNormals(final File annotatedIntervalsFile,
                                    final int expectedNumberOfEigenvalues,
                                    final File resultOutputFile) {
        try (final HDF5File hdf5PanelOfNormalsFile = new HDF5File(resultOutputFile)) {
            final SVDReadCountPanelOfNormals panelOfNormals = HDF5SVDReadCountPanelOfNormals.read(hdf5PanelOfNormalsFile);

            //check dimensions of original counts and intervals
            final RealMatrix counts = new Array2DRowRealMatrix(panelOfNormals.getOriginalReadCounts());
            Assert.assertEquals(counts.getRowDimension(), NUM_SAMPLES);
            Assert.assertEquals(counts.getColumnDimension(), NUM_INTERVALS);
            final List<SimpleInterval> originalIntervals = panelOfNormals.getOriginalIntervals();
            Assert.assertEquals(originalIntervals.size(), NUM_INTERVALS);

            //check that GC annotations are present (missing) if explicit GC correction was (was not) performed
            if (annotatedIntervalsFile != null) {
                Assert.assertEquals(panelOfNormals.getOriginalIntervalGCContent().length, NUM_INTERVALS);
            } else {
                Assert.assertEquals(panelOfNormals.getOriginalIntervalGCContent(), null);
            }

            //check filtering of samples and intervals with too many zeros
            Assert.assertEquals(panelOfNormals.getNumEigensamples(), NUM_GOOD_SAMPLES);
            Assert.assertEquals(panelOfNormals.getPanelIntervals().size(), NUM_GOOD_INTERVALS);
            Assert.assertEquals(panelOfNormals.getPanelIntervalFractionalMedians().length, NUM_GOOD_INTERVALS);

            //check that correct number of significant eigenvalues is found (this is a bit heuristic and may fail if test data is changed
            final double totalVariance = DoubleStream.of(panelOfNormals.getSingularValues()).map(x -> x * x).sum();
            final double fractionOfVarianceExplainedMissingLastEigenvalue = IntStream.range(0, expectedNumberOfEigenvalues - 1)
                    .mapToDouble(i -> panelOfNormals.getSingularValues()[i]).map(x -> x * x).sum() / totalVariance;
            Assert.assertTrue(fractionOfVarianceExplainedMissingLastEigenvalue < FRACTION_OF_VARIANCE_EXPLAINED_THRESHOLD);
            final double fractionOfVarianceExplained = IntStream.range(0, expectedNumberOfEigenvalues)
                    .mapToDouble(i -> panelOfNormals.getSingularValues()[i]).map(x -> x * x).sum() / totalVariance;
            Assert.assertTrue(fractionOfVarianceExplained > FRACTION_OF_VARIANCE_EXPLAINED_THRESHOLD);

            //check dimensions of eigenvectors
            final RealMatrix eigensampleVectors = new Array2DRowRealMatrix(panelOfNormals.getEigensampleVectors());
            Assert.assertEquals(eigensampleVectors.getRowDimension(), NUM_GOOD_INTERVALS);
            Assert.assertEquals(eigensampleVectors.getColumnDimension(), Math.min(NUMBER_OF_EIGENVALUES_REQUESTED, NUM_GOOD_SAMPLES));

            //denoise last sample (which is not a bad sample) in original counts using true number of eigenvalues
            final SimpleCountCollection sampleCounts = new SimpleCountCollection(
                    new SimpleSampleLocatableMetadata("test-sample", SEQUENCE_DICTIONARY),
                    IntStream.range(0, NUM_INTERVALS)
                            .mapToObj(i -> new SimpleCount(originalIntervals.get(i), (int) counts.getEntry(counts.getRowDimension() - 1, i)))
                            .collect(Collectors.toList()));
            final SVDDenoisedCopyRatioResult denoisedResult = panelOfNormals.denoise(sampleCounts, expectedNumberOfEigenvalues);
            //check that the denoised log2 copy ratios are sufficiently denoised
            final CopyRatioCollection denoisedCopyRatios = denoisedResult.getDenoisedCopyRatios();
            final double denoisedLog2CRStandardDeviation = new StandardDeviation().evaluate(
                    denoisedCopyRatios.getLog2CopyRatioValues().stream()
                            .mapToDouble(x -> x)
                            .toArray());
            Assert.assertTrue(denoisedLog2CRStandardDeviation < DENOISED_LOG2CR_STANDARD_DEVIATION_THRESHOLD);
            //check that the denoised log2 copy ratios are less noisy than the standardized log2 copy ratios
            final CopyRatioCollection standardizedCopyRatios = denoisedResult.getStandardizedCopyRatios();
            final double standardizedLog2CRStandardDeviation = new StandardDeviation().evaluate(
                    standardizedCopyRatios.getLog2CopyRatioValues().stream()
                            .mapToDouble(x -> x)
                            .toArray());
            Assert.assertTrue(denoisedLog2CRStandardDeviation < standardizedLog2CRStandardDeviation);

            //denoise first sample (which is a bad sample) in original counts using true number of eigenvalues
            final SimpleCountCollection badSampleCounts = new SimpleCountCollection(
                    new SimpleSampleLocatableMetadata("bad-test-sample", SEQUENCE_DICTIONARY),
                    IntStream.range(0, NUM_INTERVALS)
                            .mapToObj(i -> new SimpleCount(originalIntervals.get(i), (int) counts.getEntry(0, i)))
                            .collect(Collectors.toList()));
            final SVDDenoisedCopyRatioResult badDenoisedResult = panelOfNormals.denoise(badSampleCounts, expectedNumberOfEigenvalues);
            //check that the denoised log2 copy ratios are not sufficiently denoised
            final CopyRatioCollection badDenoisedCopyRatios = badDenoisedResult.getDenoisedCopyRatios();
            final double badDenoisedLog2CRStandardDeviation = new StandardDeviation().evaluate(
                    badDenoisedCopyRatios.getLog2CopyRatioValues().stream()
                            .mapToDouble(x -> x)
                            .toArray());
            Assert.assertFalse(badDenoisedLog2CRStandardDeviation < DENOISED_LOG2CR_STANDARD_DEVIATION_THRESHOLD);
        }
    }
}