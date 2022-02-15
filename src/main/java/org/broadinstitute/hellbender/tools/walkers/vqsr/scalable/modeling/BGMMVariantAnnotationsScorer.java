package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling;

import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clustering.BayesianGaussianMixtureModeller;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
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
    public void scoreSamples(final File inputAnnotationsFile,
                             final File outputScoresFile) {
        final List<String> inputAnnotationNames = LabeledVariantAnnotationsData.readAnnotationNames(inputAnnotationsFile);
        Utils.validateArg(inputAnnotationNames.equals(annotationNames), "Annotation names must be identical.");
        final double[][] data = LabeledVariantAnnotationsData.readAnnotations(inputAnnotationsFile);
        final double[] scores = preprocessAndScoreSamples(data).getRight();
        VariantAnnotationsScorer.writeScores(outputScoresFile, scores);
    }

    public static BGMMVariantAnnotationsScorer deserialize(final File scorerFile) {
        return deserialize(scorerFile, BGMMVariantAnnotationsScorer.class);
    }

    private Pair<double[][], double[]> preprocessAndScoreSamples(final double[][] data) {
        final double[][] preprocessedData = preprocesser.transform(data);
        final double[] scores = bgmm.scoreSamples(preprocessedData);
        return Pair.of(preprocessedData, scores);
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
