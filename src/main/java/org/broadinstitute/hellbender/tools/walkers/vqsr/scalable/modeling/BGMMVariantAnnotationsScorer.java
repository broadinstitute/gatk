package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import com.google.common.collect.ImmutableList;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.List;

public final class BGMMVariantAnnotationsScorer implements VariantAnnotationsScorer, Serializable {
    private static final long serialVersionUID = 1L;

    private final List<String> annotationNames;     // for validation
    private final BGMMVariantAnnotationsModel.Preprocesser preprocesser;
    private final BayesianGaussianMixtureModeller bgmm;

    public BGMMVariantAnnotationsScorer(final List<String> annotationNames,
                                        final BGMMVariantAnnotationsModel.Preprocesser preprocesser,
                                        final BayesianGaussianMixtureModeller bgmm) {
        this.annotationNames = ImmutableList.copyOf(annotationNames);
        this.preprocesser = preprocesser;
        this.bgmm = bgmm;
    }

    @Override
    public void score(final File inputAnnotationsFile,
                      final File outputScoresFile) {
        final List<String> inputAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        Utils.validateArg(inputAnnotationNames.equals(annotationNames), "Annotation names must be identical.");
        final double[][] annotations = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);
        final double[][] preprocessedAnnotations = preprocesser.transform(annotations);
        final double[] scores = bgmm.scoreSamples(preprocessedAnnotations);
        VariantAnnotationsScorer.writeScores(outputScoresFile, scores);
    }

    public double[][] preprocess(final double[][] annotations) {
        return preprocesser.transform(annotations);
    }

    public void serialize(final File scorerFile) {
        serialize(scorerFile, this);
    }

    public static BGMMVariantAnnotationsScorer deserialize(final File scorerFile) {
        return deserialize(scorerFile, BGMMVariantAnnotationsScorer.class);
    }

    private static <T> void serialize(final File outputFile,
                                      final T object) {
        try (final FileOutputStream fileOutputStream = new FileOutputStream(outputFile);
             final ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream)) {
            objectOutputStream.writeObject(object);
        } catch (final IOException e) {
            throw new GATKException(String.format("Exception encountered during serialization of %s to %s: %s",
                    object.getClass(), outputFile.getAbsolutePath(), e));
        }
    }

    private static <T> T deserialize(final File inputFile,
                                     final Class<T> clazz) {
        try (final FileInputStream fileInputStream = new FileInputStream(inputFile);
             final ObjectInputStream objectInputStream = new ObjectInputStream(fileInputStream)) {
            return clazz.cast(objectInputStream.readObject());
        } catch (final IOException | ClassNotFoundException e) {
            throw new GATKException(String.format("Exception encountered during deserialization from %s: %s",
                    inputFile.getAbsolutePath(), e));
        }
    }
}
