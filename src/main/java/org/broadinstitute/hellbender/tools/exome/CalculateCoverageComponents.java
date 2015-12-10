package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.pca.PCA;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Perform PCA on coverage counts.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Perform PCA on coverage counts",
        oneLineSummary = "Perform PCA on coverage counts",
        programGroup = ExomeAnalysisProgramGroup.class)
@SuppressWarnings("serial")
public class CalculateCoverageComponents extends SparkToggleCommandLineProgram {

    public static final String SAMPLES_FILE_FULL_NAME = "samplesFile";

    public static final String SAMPLES_FILE_SHORT_NAME = "samples";

    @ArgumentCollection
    protected TargetArgumentCollection targets =
            new TargetArgumentCollection(() -> this.inputFile);

    @Argument(
            doc = "PCA output file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    @Argument(
            doc = "List of samples to be considered, by default all samples are used",
            fullName = SAMPLES_FILE_FULL_NAME,
            shortName = SAMPLES_FILE_SHORT_NAME,
            optional = true
    )
    protected File sampleListFile;

    @Argument(
            doc = "Input coverage file",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputFile;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        final ReadCountCollection analysisReadCounts = composeAnalysisReadCounts();
        final List<String> targetNames = analysisReadCounts.targets().stream().map(Target::getName).collect(Collectors.toList());
        final List<String> sampleNames = analysisReadCounts.columnNames();
        final double[][] counts = analysisReadCounts.counts().getData();
        //TODO Mysterious bug when trying to use the Spark based SVD
        //TODO Corresponding issue number #242.
        //TODO for now force to use the non-spark implementation:
        //TODO once solved uncomment the following 2 lines and remove the other two that follow.
        //final PCA pca = PCA.createPCA(targetNames, sampleNames, new Array2DRowRealMatrix(counts),
        //        ctx == null ? SVDFactory::createSVD : (dm) -> SVDFactory.createSVD(dm, ctx));
        final PCA pca = PCA.createPCA(targetNames, sampleNames, new Array2DRowRealMatrix(counts),
                SVDFactory::createSVD);
        try (final HDF5File output = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            PCA.writeHDF5(pca, output);
        } catch (final GATKException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    private ReadCountCollection composeAnalysisReadCounts() {
        final ReadCountCollection inputReadCounts;
        try {
            inputReadCounts = ReadCountCollectionUtils.parse(inputFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, "cannot reach or read the input target coverage file");
        }
        final List<String> inputSamples = inputReadCounts.columnNames();
        final Set<String> outputSamples = composeAnalysisSampleSet(inputSamples);

        final List<Target> inputTargets = inputReadCounts.targets();

        final Set<Target> outputTargets = composeAnalysisTargetsSet(inputTargets);
        return inputReadCounts.subsetColumns(outputSamples).subsetTargets(outputTargets);
    }

    private Set<Target> composeAnalysisTargetsSet(final List<Target> inputTargets) {
        final TargetCollection<Target> requestedTargets = targets.readTargetCollection(false);
        if (requestedTargets == null) {
            return new HashSet<>(inputTargets);
        } else if (requestedTargets.targets().stream().anyMatch(t -> !inputTargets.contains(t))) {
            throw new UserException.BadInput(
                    "some of the requested targets are not present in the input read counts; e.g. " +
                requestedTargets.targets().stream()
                        .filter(t -> !inputTargets.contains(t))
                        .limit(5)
                        .map(Target::getName)
                        .collect(Collectors.joining(", ")));
        } else {
            return requestedTargets.targets().stream().collect(Collectors.toSet());
        }
    }

    private Set<String> composeAnalysisSampleSet(final List<String> inputSamples) {
        if (sampleListFile == null) {
            return new HashSet<>(inputSamples);
        } else {
            final Set<String> listedSamples;
            try (final BufferedReader reader = new BufferedReader(new FileReader(sampleListFile))) {
                listedSamples = reader.lines().collect(Collectors.toSet());
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(sampleListFile, "cannot read the sample list file", ex);
            }
            if (listedSamples.stream().anyMatch(s -> !inputSamples.contains(s))) {
                throw new UserException.BadInput("the sample list provided contains sample names not present in the input counts; e.g. " +
                    listedSamples.stream().filter(s -> !inputSamples.contains(s)).limit(5).collect(Collectors.joining(", ")));
            }
            return listedSamples;
        }
    }
}
