package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.hipparchus.stat.fitting.EmpiricalDistribution;

import java.io.File;
import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * File interface for passing annotations to a scoring backend and returning scores.
 */
public interface VariantAnnotationsScorer {

    String SCORES_PATH = "/data/scores"; // our HDF5 library does not allow writing to a bare/root path (e.g., /scores)

    /**
     * @param inputAnnotationsFile  Annotations to be scored in HDF5 format, containing at least the directory structure
     *
     *                              <p>
     *                                |--- annotations<br>
     *                                |    |--- chunk_0<br>
     *                                |    |--- ...<br>
     *                                |    |--- chunk_{num_chunks - 1}<br>
     *                                |    |--- names<br>
     *                                |    |--- num_chunks<br>
     *                                |    |--- num_columns<br>
     *                                |    |--- num_rows<br>
     *                              </p>
     *
     *                              Here, each chunk is a double matrix, with dimensions given by
     *                              (number of sites in the chunk) x (number of annotations).
     *                              See {@link LabeledVariantAnnotationsData#writeHDF5}.
     *
     *                              Scoring backends are responsible for consuming annotations in this format and
     *                              outputting a double array of scores to file. This responsibility includes
     *                              validation of annotation names.
     *
     * @param outputScoresFile      Output file in HDF5 format, containing scores at {@link VariantAnnotationsScorer#SCORES_PATH}.
     */
    void score(final File inputAnnotationsFile,
               final File outputScoresFile);

    /**
     * Given scores for a calibration set, returns a function for converting a subsequent score to a
     * sensitivity to that calibration set. This function is simply given by 1 - ECDF,
     * where ECDF is the empirical cumulative distribution function of the calibration scores;
     * see <a href='https://en.wikipedia.org/wiki/Empirical_distribution_function'>here</a>.
     * For example, a score that is very low relative to the calibration scores would yield a
     * high calibration sensitivity; that is, using this score as the minimum allowable threshold for filtering
     * would result in a high sensitivity to the calibration set.
     *
     * @param calibrationScores must all be finite
     */
    static Function<Double, Double> createScoreToCalibrationSensitivityConverter(final double[] calibrationScores) {
        Utils.validateArg(Arrays.stream(calibrationScores).allMatch(Double::isFinite),
                "Calibration scores must all be finite.");
        final EmpiricalDistribution empiricalDistribution = new EmpiricalDistribution();
        empiricalDistribution.load(calibrationScores);
        return score -> 1. - empiricalDistribution.cumulativeProbability(score);
    }

    /**
     * Reads a double array of scores from {@value SCORES_PATH} in an HDF5 file.
     */
    static double[] readScores(final File inputFile) {
        try (final HDF5File inputHDF5File = new HDF5File(inputFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(inputHDF5File.getFile());
            return inputHDF5File.readDoubleArray(SCORES_PATH);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                    inputFile.getAbsolutePath(), exception));
        }
    }

    /**
     * Writes a double array of scores to {@value SCORES_PATH} in an HDF5 file.
     */
    static void writeScores(final File outputFile,
                            final double[] scores) {
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            outputHDF5File.makeDoubleArray(SCORES_PATH, scores);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of scores (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
    }

    /**
     * Yields a VQSR-style positive-negative scorer that returns the difference of the positive score and the negative score.
     */
    static VariantAnnotationsScorer combinePositiveAndNegativeScorer(final VariantAnnotationsScorer positiveScorer,
                                                                     final VariantAnnotationsScorer negativeScorer) {
        return (inputAnnotationsFile, outputScoresFile) -> {
            final File tempPositiveScoresFile = IOUtils.createTempFile("positive", "scores.hdf5");
            final File tempNegativeScoresFile = IOUtils.createTempFile("negative", "scores.hdf5");
            positiveScorer.score(inputAnnotationsFile, tempPositiveScoresFile);
            final double[] positiveScores = VariantAnnotationsScorer.readScores(tempPositiveScoresFile);
            negativeScorer.score(inputAnnotationsFile, tempNegativeScoresFile);
            final double[] negativeScores = VariantAnnotationsScorer.readScores(tempNegativeScoresFile);
            final double[] scores = IntStream.range(0, positiveScores.length).mapToDouble(i -> positiveScores[i] - negativeScores[i]).toArray();
            VariantAnnotationsScorer.writeScores(outputScoresFile, scores);
        };
    }
}
