package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

final class SerializationUtils {

    private static final Logger logger = LogManager.getLogger(SerializationUtils.class);

    static <T> void serialize(final File outputFile,
                              final T object) {
        try (final FileOutputStream fileOutputStream = new FileOutputStream(outputFile);
             final ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream)) {
            objectOutputStream.writeObject(object);
        } catch (final IOException e) {
            throw new GATKException(String.format("Exception encountered during serialization of %s to %s: %s",
                    object.getClass(), outputFile.getAbsolutePath(), e));
        }
    }

    static <T> T deserialize(final File inputFile,
                             final Class<T> clazz) {
        try (final FileInputStream fileInputStream = new FileInputStream(inputFile);
             final ObjectInputStream objectInputStream = new ObjectInputStream(fileInputStream)) {
            return clazz.cast(objectInputStream.readObject());
        } catch (final IOException | ClassNotFoundException e) {
            throw new GATKException(String.format("Exception encountered during deserialization from %s: %s",
                    inputFile.getAbsolutePath(), e));
        }
    }

    static double[] readScores(final File inputFile) {
        try (final HDF5File inputHDF5File = new HDF5File(inputFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(inputHDF5File.getFile());
            return inputHDF5File.readDoubleArray("/data/scores");
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of scores from %s: %s",
                    inputFile.getAbsolutePath(), exception));
        }
    }
    
    static void writeScores(final File outputFile,
                            final double[] scores) {
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            outputHDF5File.makeDoubleArray("/data/scores", scores);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of scores (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
        logger.info(String.format("Scores written to %s.", outputFile.getAbsolutePath()));
    }

//    // TODO we could just get all this stuff from the VariantDataManager, but let's clean that up later
//    static void writeTruthSensitivityTranches(final File outputFile,
//                                              final File scoresFile,
//                                              final File annotationsFile,
//                                              final List<Double> truthSensitivityTranches,
//                                              final VariantType mode) {
//        // TODO validate
//        final List<Double> scores = Doubles.asList(readScores(scoresFile));
//        final List<Boolean> isTruth = VariantLabeledAnnotationsData.readLabel(annotationsFile, VariantLabeledAnnotationsData.TRUTH_LABEL);
//        try {
//            final PrintStream tranchesStream = new PrintStream(outputFile);
//
//            // Find the score cutoff values which correspond to the various tranches of calls requested by the user
//            final int nCallsAtTruth = TruthSensitivityTranche.countCallsAtTruth(scores, isTruth, Double.NEGATIVE_INFINITY);
//            final TruthSensitivityTranche.TruthSensitivityMetric metric = new TruthSensitivityTranche.TruthSensitivityMetric(nCallsAtTruth);
//            final List<TruthSensitivityTranche> tranches = TruthSensitivityTranche.findTranches(scores, isTruth, truthSensitivityTranches, metric, mode);
//            tranchesStream.print(TruthSensitivityTranche.printHeader());
//            tranchesStream.print(TruthSensitivityTranche.tranchesString(tranches));
//        } catch (final FileNotFoundException e) {
//            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
//        }
//        logger.info(String.format("Truth-sensitivity tranches written to %s.", outputFile.getAbsolutePath()));
//    }
}
