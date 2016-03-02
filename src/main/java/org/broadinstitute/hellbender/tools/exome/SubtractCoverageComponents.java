package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pca.PCA;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Tool to remove principal components from coverage data.
 * <p>
 * This tool takes two inputs: the coverage to be normalized (argument {@link #inputFile}) and the {@link CalculateCoverageComponents PCA} run
 * tool output file (argument {@link #pcaFile}).
 * </p>
 * <p>
 * Syntax example:
 * <pre>
 *     java -jar gatk-protected.jar CalculateCoverageComponents -I my-panel-coverage.tab -o pca.hd5
 *     java -jar gatk-protected.jar SubtractCoverageComponents -I my-case-coverage.tab -pca pca.hd5 -o my-case-coverage-normalized.tab
 * </pre>
 * </p>
 * <p>
 * By default this tool will remove all the components present in the PCA input file. Nevertheless often the user will want to
 * indicate what components to subtract and what components to keep.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
    summary = "Normalizes coverage counts by subtracting principal coverage components",
    oneLineSummary = "Normalizes coverage by subtracting principal components",
    programGroup = CopyNumberProgramGroup.class
)
public final class SubtractCoverageComponents extends CommandLineProgram {

    public static final String PCA_INPUT_FULL_NAME = "principalComponentsFile";

    public static final String PCA_INPUT_SHORT_NAME = "pca";

    public static final String NUM_COMPONENTS_FULL_NAME = "numberOfComponents";

    public static final String NUM_COMPONENTS_SHORT_NAME = "top";

    @Argument(
            doc = "Input target coverage to normalize",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputFile;

    @Argument(
            doc = "PCA input file",
            fullName = PCA_INPUT_FULL_NAME,
            shortName = PCA_INPUT_SHORT_NAME,
            optional = false
    )
    protected File pcaFile;

    @Argument(
            doc = "Normalized target coverage output",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    @Argument(
            doc = "Number of principal components to use",
            fullName = NUM_COMPONENTS_FULL_NAME,
            shortName = NUM_COMPONENTS_SHORT_NAME,
            optional = true
    )
    protected int numberOfComponentsToUse = Integer.MAX_VALUE;

    @Override
    protected Object doWork() {
        ParamUtils.isPositiveOrZero(numberOfComponentsToUse, "number of components to use must be positive or zero.");
        final PCA pca = readPcaInputFile(pcaFile);
        final ReadCountCollection coverage = readInputCoverage(inputFile);

        // Match up the targets in the PCA and the input counts file.
        final Map<String, Integer> pcaTargetIndexByName = composeComponentTargetIndexMap(pca);
        final Map<String, Integer> coverageTargetIndexByName = composeCoverageTargetIndexMap(coverage);
        compareCompomentsAndCoverageTargetNames(pcaTargetIndexByName, coverageTargetIndexByName);

        // The actual work:

        // T x C matrix.
        final RealMatrix eigenVectors = pca.getEigenVectors();
        final int componentCount = eigenVectors.getColumnDimension();
        final RealVector centers = pca.getCenters();
        // S x T matrix with the counts ready to be multiplied by the eigen-vectors T x C
        final RealMatrix transposedCoverage = composeCoverageMatrixReadyForNormalization(pcaTargetIndexByName, coverageTargetIndexByName, centers, coverage);
        // S x C matrix where each row indicate the projection value of the sample onto a component.
        final RealMatrix projection = transposedCoverage.multiply(eigenVectors);

        // T x S matrix with the normalized new counts.
        final RealMatrix resultMatrix = subtractComponentProjectionsAndTranspose(projection, transposedCoverage, eigenVectors);

        // Creates the output read count collection.
        final ReadCountCollection normalizedCoverage = new ReadCountCollection(
                SetUniqueList.setUniqueList(composeTargetList(pca, coverage, coverageTargetIndexByName)),
                SetUniqueList.setUniqueList(new ArrayList<>(coverage.columnNames())), resultMatrix);

        // Output the result normalized coverage.
        writeOutputCoverage(normalizedCoverage);

        return "SUCCESS";
    }

    private ReadCountCollection readInputCoverage(final File inputFile) {
        final ReadCountCollection coverage;
        try {
            coverage = ReadCountCollectionUtils.parse(inputFile);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, ex);
        }
        return coverage;
    }

    private PCA readPcaInputFile(final File pcaFile) {
        final PCA pca;
        try (final HDF5File pcaHd5 = new HDF5File(pcaFile, HDF5File.OpenMode.READ_ONLY)) {
            pca = PCA.readHDF5(pcaHd5);
        }
        return pca;
    }

    private void writeOutputCoverage(final ReadCountCollection normalizedCoverage) {
        try {
            ReadCountCollectionUtils.write(outputFile, normalizedCoverage);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }
    }

    private List<Target> composeTargetList(final PCA pca, final ReadCountCollection coverage, final Map<String, Integer> coverageTargetIndexByName) {
        return pca.getVariables().stream()
                .mapToInt(coverageTargetIndexByName::get)
                .mapToObj(i -> coverage.targets().get(i))
                .collect(Collectors.toList());
    }

    /**
     * @param projection S (sample) x C (component) matrix with the length of the projection
     * @param transposedCoverage S (sample) x T (target) matrix with the coverage pre normalization.
     * @return never {@code null}, a S(sample) x T
     */
    private RealMatrix subtractComponentProjectionsAndTranspose(final RealMatrix projection, final RealMatrix transposedCoverage, final RealMatrix eigenVectors) {
        final int sampleCount = projection.getRowDimension();
        final int targetCount = transposedCoverage.getColumnDimension();
        final int componentCount = Math.min(numberOfComponentsToUse, projection.getColumnDimension());
        final RealMatrix result = new Array2DRowRealMatrix(targetCount, sampleCount);
        for (int i = 0; i < sampleCount; i++) {
            final double[] newCoverage = transposedCoverage.getRow(i);
            for (int j = 0; j < componentCount; j++) {
                final double[] eigenVector = eigenVectors.getColumn(j);
                final double factor = projection.getEntry(i, j);
                for (int k = 0; k < targetCount; k++) {
                    newCoverage[k] -= eigenVector[k] * factor;
                }
            }
            result.setColumn(i, newCoverage);
        }
        return result;
    }

    /**
     * Creates a map between target name and row index in the PCA result file.
     * @param pca the PCA result object.
     * @return never {@code null}.
     * @throws UserException.BadInput if:
     * <ul>
     *     <li>the input PCA object does not have assigned variable names,</li>
     *     <li>or the number of variables is 0,</li>
     *     <li>or any of the variable names is null</li>
     *     <li>or any of the variable names is duplicated.</li>
     * </ul>
     */
    private Map<String, Integer> composeComponentTargetIndexMap(final PCA pca) {
        final Map<String, Integer> result = composeTargetIndexByNameMap(pca.getVariables(), "input pca file");
        if (result.size() == 0) {
            throw new UserException.BadInput("there are 0 targets in the PCA file");
        }
        return result;
    }

    /**
     * Creates a map between target name and row index in the coverage input file.
     * @param coverage the input coverage collection.
     * @return never {@code null}.
     * @throws UserException.BadInput if:
     * <ul>
     *     <li>the input coverage file does not have assigned targets,</li>
     *     <li>or the number of targets is 0,</li>
     *     <li>or any of the target names is null</li>
     *     <li>or any of the target names is duplicated.</li>
     * </ul>
     */
    private Map<String, Integer> composeCoverageTargetIndexMap(final ReadCountCollection coverage) {
        final List<Target> targets = coverage.targets();
        return composeTargetIndexByNameMap(targets == null ? null :
                targets.stream().map(Target::getName).collect(Collectors.toList()), "input coverage file");
    }

    /**
     * Returns the coverage values that correspond to a count column as a 1-row matrix
     * ready to be normalized.
     *
     * <p>
     *     The resulting matrix has one row in the order they appear in the PCA file,
     *     removing targets that are not present in such a file.
     * </p>
     * @param pcaTargetIndexByName a map from the target name to its index in the PCA.
     * @param coverageTargetIndexByName a map from the target name to its index in the coverage matrix.
     * @param coverage the coverage matrix.
     * @return never {@code null}.
     */
    private RealMatrix composeCoverageMatrixReadyForNormalization(final Map<String, Integer> pcaTargetIndexByName,
        final Map<String, Integer> coverageTargetIndexByName, final RealVector centers, final ReadCountCollection coverage) {
        final RealMatrix coverageValues = coverage.counts();

        final RealMatrix result = new Array2DRowRealMatrix(coverage.columnNames().size(), pcaTargetIndexByName.size());
        final int[] coverageToComponentTargetIndexes = pcaTargetIndexByName.keySet().stream()
                .sorted(Comparator.comparing(pcaTargetIndexByName::get))
                .mapToInt(coverageTargetIndexByName::get)
                .toArray();
        for (int i = 0; i < result.getColumnDimension(); i++) {
           final int coverageIndex = coverageToComponentTargetIndexes[i];
           final double[] coverageRowValues = coverageValues.getRow(coverageIndex);
           final double center = centers.getEntry(i);
           for (int j = 0; j < coverageRowValues.length; j++) {
               coverageRowValues[j] -= center;
           }
           result.setColumn(i, coverageRowValues);
        }
        return result;
    }

    private void compareCompomentsAndCoverageTargetNames(final Map<String, Integer> componentTargets,
                                                         final Map<String, Integer> coverageTargetIndexByName) {
        final List<String> componentTargetsMissingCoverage = componentTargets.keySet().stream()
                .filter(n -> !coverageTargetIndexByName.containsKey(n)).collect(Collectors.toList());

        if (componentTargetsMissingCoverage.size() > 0) {
            throw new UserException.BadInput("the input coverage file is missing some of the targets used to calculate "
                    + "the components; e.g. " + componentTargetsMissingCoverage.stream().limit(5).collect(Collectors.joining(", ")));
        }

        final List<String> coverageTargetsMissingComponents = coverageTargetIndexByName.keySet().stream()
                .filter(n -> !componentTargets.containsKey(n)).collect(Collectors.toList());
        if (coverageTargetsMissingComponents.size() > 0) {
            logger.warn(
                    String.format("there are some targets in the coverage input file that are not present in the components file (a total of %d). These won't be normalized nor output", coverageTargetsMissingComponents.size()));
        }
    }

    private Map<String, Integer> composeTargetIndexByNameMap(final List<String> targetNames, final String inputName) {
        if (targetNames == null) {
            throw new UserException.BadInput(String.format("the %s does not contain target names", inputName));
        }
        final Map<String, Integer> result = IntStream.range(0, targetNames.size())
                .boxed()
                .collect(Collectors.toMap(targetNames::get, i -> i));

        if (result.containsKey(null)) {
            throw new UserException.BadInput(String.format("the %s contains null target names", inputName));
        } else if (result.size() != targetNames.size()) {
            throw new UserException.BadInput("the %s contains repeated target names");
        }
        return result;
    }
}
