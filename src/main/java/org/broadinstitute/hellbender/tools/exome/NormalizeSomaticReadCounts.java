package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.Level;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.HDF5Reader;
import org.broadinstitute.hellbender.utils.hdf5.PoN;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Normalizes read counts given the PanelOfNormals.
 * <p/>
 * <p>
 * Currently this tool only performs target factor normalization.
 * </p>
 * <p/>
 * <p>
 * Subsequent version will perform tangent normalization as well.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        usage = "Count overlapping reads exon by exon",
        usageShort = "Count overlapping reads exon by exon",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class NormalizeSomaticReadCounts extends CommandLineProgram {

    protected final static String READ_COUNTS_FILE_FULL_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    protected final static String READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    protected final static String TARGET_FILE_FULL_NAME = "targets";
    protected final static String TARGET_FILE_SHORT_NAME = "T";

    protected final static String PON_FILE_FULL_NAME = "panelOfNormals";
    protected final static String PON_FILE_SHORT_NAME = "pon";

    protected final static String FACTOR_NORMALIZED_COUNTS_FULL_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    protected final static String FACTOR_NORMALIZED_COUNTS_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    @Argument(
            doc = "read counts input file.",
            shortName = READ_COUNTS_FILE_SHORT_NAME,
            fullName = READ_COUNTS_FILE_FULL_NAME,
            optional = false
    )
    protected File readCountsFile;

    @Argument(
            doc = "target BED file.",
            shortName = TARGET_FILE_SHORT_NAME,
            fullName = TARGET_FILE_FULL_NAME,
            optional = true
    )
    protected File targetFile;

    @Argument(
            doc = "panel Of normals HDF5 file",
            shortName = PON_FILE_SHORT_NAME,
            fullName = PON_FILE_FULL_NAME,
            optional = false
    )
    protected File ponFile;

    @Argument(
            doc = "Factor normalized counts output",
            shortName = FACTOR_NORMALIZED_COUNTS_SHORT_NAME,
            fullName = FACTOR_NORMALIZED_COUNTS_FULL_NAME,
            optional = false
    )
    protected File fntOutFile;

    @Override
    protected Object doWork() {

        Utils.regularReadableUserFile(ponFile);
        try (final HDF5Reader ponReader = new HDF5Reader(ponFile)) {
            final PoN pon = new HDF5PoN(ponReader);
            final ExonCollection<? extends BEDFeature> exonCollection = readExonCollection(targetFile);
            final ReadCountCollection readCountCollection = readInputReadCounts(readCountsFile, exonCollection);
            final ReadCountCollection targetFactorNormalized = targetFactorNormalization(readCountCollection, pon);
            writeOutput(targetFactorNormalized);
            return "SUCCESS";
        }
    }

    private ExonCollection<? extends BEDFeature> readExonCollection(final File targetFile) {
        if (targetFile == null) {
            return null;
        } else {
            Utils.regularReadableUserFile(targetFile);
            final FeatureCodec<?, ?> codec = FeatureManager.getCodecForFile(targetFile);
            logger.log(Level.INFO, String.format("Reading target intervals from exome file '%s' ...", new Object[]{targetFile.getAbsolutePath()}));
            final Class<?> featureType = codec.getFeatureType();
            if (BEDFeature.class.isAssignableFrom(featureType)) {
                @SuppressWarnings("unchecked")
                final FeatureCodec<? extends BEDFeature, ?> bedCodec = (FeatureCodec<? extends BEDFeature, ?>) codec;
                final ExonCollection<? extends BEDFeature> result = ExonCollections.fromBEDFeatureFile(targetFile, bedCodec);
                logger.log(Level.INFO, String.format("Found %d targets to analyze.", new Object[]{Integer.valueOf(result.exonCount())}));
                return result;
            } else {
                throw new UserException.BadInput(String.format("currently only BED formatted target file are supported. '%s' does not seem to be a BED file",
                        new Object[]{targetFile.getAbsolutePath()}));
            }
        }
    }

    /**
     * Reads the read-counts from the input file using an exon collection (if provided) to resolve target names if this
     * are missing.
     *
     * @param readCountsFile the read-counts file.
     * @param exonCollection the input exon collection. {@code null} indicates that no collection was provided by the user.
     * @param <E>            the element type of {@code exonCollection}, it must be a {@link BEDFeature}.
     * @return never {@code null}.
     * @throws UserException.CouldNotReadInputFile if there was some problem
     *                                             trying to read from {@code readCountsFile}.
     * @throws UserException.BadInput              if there is some format issue with the input file {@code readCountsFile},
     *                                             or it does not contain target names and {@code exonCollection} is {@code null}.
     */
    private <E extends BEDFeature> ReadCountCollection readInputReadCounts(final File readCountsFile,
                                                                           final ExonCollection<E> exonCollection) {
        try {
            return ReadCountCollectionUtils.parse(readCountsFile, exonCollection);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(readCountsFile, ex.getMessage(), ex);
        }
    }

    private void writeOutput(final ReadCountCollection targetFactorNormalized) {
        try {
            ReadCountCollectionUtils.write(fntOutFile, targetFactorNormalized,
                    "#fileFormat = tsv",
                    "#commandLine = " + getCommandLine(),
                    "#title = Factor normalized target counts");
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(fntOutFile, ex.getMessage());
        }
    }

    /**
     * Normalize a read-count collection by the PoN target factors.
     * <p>
     * The result read-count collection only contains entries for target contained in the input {@code PoN},
     * and target are sorted as they are in the PoN.
     * </p>
     * <p>
     * The number of columns and their order correspond to the input read-count collection.
     * </p>
     *
     * @param input list of output target information.
     * @return never {@code null}.
     */
    private ReadCountCollection targetFactorNormalization(final ReadCountCollection input, final PoN pon) {

        final List<String> ponTargetNames = pon.targetNames();
        final List<Target> inputTargets = input.targets();
        final Map<String, Integer> inputTargetIndexes = IntStream.range(0, inputTargets.size()).boxed()
                .collect(Collectors.toMap(i -> inputTargets.get(i).getName(), Function.identity()));
        final RealMatrix inputCounts = input.counts();

        // Reduce the input counts to match the PoN targets with the same order.
        final RealMatrix normalizationInput = new Array2DRowRealMatrix(ponTargetNames.size(), inputCounts.getColumnDimension());
        for (int i = 0; i < ponTargetNames.size(); i++) {
            final String targetName = ponTargetNames.get(i);
            final Integer inputIndex = inputTargetIndexes.get(targetName);
            if (inputIndex == null) {
                throw new UserException.BadInput(String.format("missing PoN target name %s in input read counts %s", targetName, readCountsFile));
            }
            normalizationInput.setRow(i, inputCounts.getRow(inputIndex));
        }
        // Actual normalization step.
        final RealMatrix normalizationOutput = pon.factorNormalization(normalizationInput);

        // Calculate final target output order to match the inputs target order:
        final List<Target> outputTargets = ponTargetNames.stream()
                .map(inputTargetIndexes::get).sorted().map(inputTargets::get).collect(Collectors.toList());

        final RealMatrix outputCounts = rearrangeNormalizedOutputCountRows(normalizationOutput, ponTargetNames, outputTargets);

        // Compose the output normalized counts and return.
        return new ReadCountCollection(SetUniqueList.setUniqueList(outputTargets), SetUniqueList.setUniqueList(new ArrayList<>(input.columnNames())), outputCounts);
    }

    private RealMatrix rearrangeNormalizedOutputCountRows(final RealMatrix normalizationOutput, final List<String> ponTargetNames, final List<Target> outputTargets) {
        final Map<String,Integer> ponTargetsIndexes = IntStream.range(0, ponTargetNames.size()).boxed()
                .collect(Collectors.toMap(i -> ponTargetNames.get(i), Function.identity()));
        final RealMatrix outputCounts = new Array2DRowRealMatrix(normalizationOutput.getRowDimension(),normalizationOutput.getColumnDimension());
        for (int i = 0; i < outputTargets.size(); i++) {
            final Target target = outputTargets.get(i);
            outputCounts.setRow(i,normalizationOutput.getRow(ponTargetsIndexes.get(target.getName())));
        }
        return outputCounts;
    }
}
