package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.hipparchus.stat.fitting.EmpiricalDistribution;

import java.io.File;
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
     *                                ├── annotations
     *                                     ├── chunk_0
     *                                     ├── ...
     *                                     ├── chunk_{num_chunks - 1}
     *                                     ├── names
     *                                     ├── num_chunks
     *                                     ├── num_columns
     *                                     └── num_rows
     *
     *                              Scoring backends are responsible for consuming annotations in this format.
     * @param outputScoresFile      Output file in HDF5 format, containing scores at {@link VariantAnnotationsScorer#SCORES_PATH}.
     */
    void score(final File inputAnnotationsFile,
               final File outputScoresFile);

    // TODO document that this is 1 - ECDF
    static Function<Double, Double> createScoreToCalibrationSensitivityConverter(final double[] calibrationScores) {
        // TODO validate
        final EmpiricalDistribution empiricalDistribution = new EmpiricalDistribution();
        empiricalDistribution.load(calibrationScores);
        return score -> 1. - empiricalDistribution.cumulativeProbability(score);
    }

    static double[] readScores(final File inputFile) {
        try (final HDF5File inputHDF5File = new HDF5File(inputFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(inputHDF5File.getFile());
            return inputHDF5File.readDoubleArray(SCORES_PATH);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                    inputFile.getAbsolutePath(), exception));
        }
    }

    static void writeScores(final File outputFile,
                            final double[] scores) {
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            outputHDF5File.makeDoubleArray(SCORES_PATH, scores);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of scores (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
    }

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
