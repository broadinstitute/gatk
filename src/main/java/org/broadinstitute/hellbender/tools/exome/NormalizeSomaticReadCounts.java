package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.logging.log4j.Level;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExomeAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.hdf5.PoN;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

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
        summary = "Normalizes PCOV read counts using a panel of normals",
        oneLineSummary = "Normalizes proportional coverage (PCOV) read counts using a panel of normals",
        programGroup = ExomeAnalysisProgramGroup.class
)
public final class NormalizeSomaticReadCounts extends CommandLineProgram {

    /**
     * Name of the column that contains the I.D. of the PoN <i>"eigen sample"</i>.
     */
    public static final String PON_SAMPLE_BETA_HAT_COLUMN_NAME = "PON.SAMPLE";

    /**
     * Minimum target normalized and column centered count possible.
     *
     * <p>
     *     It must be small yet greater than 0 to avoid -Inf problems in the calculations.
     * </p>
     */
    public static double EPSILON = Math.pow(10,-10);

    /**
     * Cached inverse of the natural logarithm of 2.
     */
    private static double INV_LN2 = 1.0 / Math.log(2.0);

    public static final String READ_COUNTS_FILE_FULL_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String TARGET_FILE_FULL_NAME = "targets";
    public static final String TARGET_FILE_SHORT_NAME = "T";

    public static final String PON_FILE_FULL_NAME = "panelOfNormals";
    public static final String PON_FILE_SHORT_NAME = "pon";

    public static final String TANGENT_NORMALIZED_COUNTS_FULL_NAME = StandardArgumentDefinitions.OUTPUT_LONG_NAME;
    public static final String TANGENT_NORMALIZED_COUNTS_SHORT_NAME = StandardArgumentDefinitions.OUTPUT_SHORT_NAME;

    public static final String FACTOR_NORMALIZED_COUNTS_FULL_NAME = "factorNormalizedOutput";
    public static final String FACTOR_NORMALIZED_COUNTS_SHORT_NAME = "FNO";

    public static final String TANGENT_BETA_HATS_FULL_NAME = "betaHatsOutput";
    public static final String TANGENT_BETA_HATS_SHORT_NAME = "BHO";

    public static final String PRE_TANGENT_NORMALIZATION_FULL_NAME = "preTangentNormalizationOutput";
    public static final String PRE_TANGENT_NORMALIZATION_SHORT_NAME = "PTNO";

    @Argument(
            doc = "read counts input file.  This can only contain one sample at a time.",
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
            optional = true
    )
    protected File fntOutFile;

    @Argument(
            doc = "Tangent normalized counts output",
            shortName = TANGENT_NORMALIZED_COUNTS_SHORT_NAME,
            fullName = TANGENT_NORMALIZED_COUNTS_FULL_NAME,
            optional = false
    )
    protected File outFile;

    @Argument(
            doc = "Tangent normalization Beta Hats output file",
            shortName = TANGENT_BETA_HATS_SHORT_NAME,
            fullName = TANGENT_BETA_HATS_FULL_NAME,
            optional = true
    )
    protected File betaHatsOutFile;

    @Argument(
            doc = "Pre-tangent normalization counts",
            shortName = PRE_TANGENT_NORMALIZATION_SHORT_NAME,
            fullName = PRE_TANGENT_NORMALIZATION_FULL_NAME,
            optional = true
    )
    protected File preTangentNormalizationOutFile;

    @Override
    protected Object doWork() {

        Utils.regularReadableUserFile(ponFile);
        try (final HDF5File ponReader = new HDF5File(ponFile)) {
            final PoN pon = new HDF5PoN(ponReader);
            final TargetCollection<? extends BEDFeature> exonCollection = readTargetCollection(targetFile);
            final ReadCountCollection readCountCollection = readInputReadCounts(readCountsFile, exonCollection);
            final ReadCountCollection targetFactorNormalized = applyTargetFactorNormalization(pon, readCountCollection);
            applyTangentNormalization(pon, readCountCollection, targetFactorNormalized);
            return "SUCCESS";
        }
    }

    /**
     * Perform tangent normalization.
     *
     * <p>
     * This is done on target-factor normalized counts provided as the input.
     * </p>
     *
     * @param pon the panel-of-normals to use to tangent-normalize
     * @param originalCounts the original input read-counts.
     * @param targetFactorNormalizedCounts the target-factor normalized counts.
     */
    private void applyTangentNormalization(final PoN pon, final ReadCountCollection originalCounts, final ReadCountCollection targetFactorNormalizedCounts) {

        final Case2PoNTargetMapper targetMapper = new Case2PoNTargetMapper(targetFactorNormalizedCounts.targets(), pon.getPanelTargetNames(), readCountsFile);

        // The input counts with row (targets) sorted so that the match the PoN's matrix order.
        final RealMatrix tangentNormalizationRawInputCounts = targetMapper.fromCaseToPoNCounts(targetFactorNormalizedCounts.counts());

        // We prepare the counts for tangent normalization.
        final RealMatrix tangentNormalizationInputCounts = composeTangentNormalizationInputMatrix(tangentNormalizationRawInputCounts);

        // Calculate the beta-hats for the input read count columns (samples).
        final RealMatrix tangentBetaHats = pon.betaHats(tangentNormalizationInputCounts, true, EPSILON);

        // Actual tangent normalization step.
        final RealMatrix tangentNormalizedCounts = pon.tangentNormalization(tangentNormalizationInputCounts, tangentBetaHats, true);

        // Output the tangent normalized counts.
        final ReadCountCollection tangentNormalized = targetMapper.fromPoNtoCaseCountCollection(tangentNormalizedCounts, originalCounts.columnNames());
        writeTangentNormalizedOutput(tangentNormalized);
        final ReadCountCollection preTangentNormalized = targetMapper.fromPoNtoCaseCountCollection(tangentNormalizationInputCounts, originalCounts.columnNames());
        writePreTangentNormalizationOutput(preTangentNormalized);
        writeTangentBetaHats(tangentBetaHats, originalCounts.columnNames());
    }

    /**
     * Writes the pre-tangent-normalization read counts if a file was provided for it.
     *
     * Please note that this file is written in the target format originally used in recapseg.  In other words, it is
     *  written as if it was target coverage. If the output file is null, we do not write a file.
     *
     * @param preTangentNormalized the read count collection to write.
     */
    private void writePreTangentNormalizationOutput(final ReadCountCollection preTangentNormalized) {
        /**
         * If the output file is null, we do not write a file.
         */
        if (preTangentNormalizationOutFile != null) {
            try {
                ReadCountCollectionUtils.writeAsTargetCoverage(preTangentNormalizationOutFile, preTangentNormalized, "fileFormat = tsv",
                        "commandLine = " + getCommandLine(),
                        "title = Pre tangent normalized coverage profile");
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(preTangentNormalizationOutFile, ex.getMessage());
            }
        }
    }

    /**
     * Perform target-factor normalization of the original read counts.
     * @param pon the panel-of-normals.
     * @param inputReadCounts the input/original read counts.
     * @return never {@code null} a new read counts collections.
     */
    private ReadCountCollection applyTargetFactorNormalization(final PoN pon, final ReadCountCollection inputReadCounts) {

        final Case2PoNTargetMapper targetMapper = new Case2PoNTargetMapper(inputReadCounts.targets(), pon.getTargetNames(), readCountsFile);
        final RealMatrix inputCounts = targetMapper.fromCaseToPoNCounts(inputReadCounts.counts());
        final RealMatrix targetNormalizedCounts = pon.factorNormalization(inputCounts);

        final ReadCountCollection targetFactorNormalized = targetMapper.fromPoNtoCaseCountCollection(targetNormalizedCounts, inputReadCounts.columnNames());

        writeTargetFactorNormalizedOutput(targetFactorNormalized);
        return targetFactorNormalized;
    }

    /**
     * Write the beta-hats if an output file was provided for it.
     *
     * @param tangentBetaHats the beta-hats to write.
     * @param countColumnNames the input read count column names.
     */
    private void writeTangentBetaHats(final RealMatrix tangentBetaHats, final List<String> countColumnNames) {
        if (betaHatsOutFile == null) {
            return;
        }

        final List<String> columnNames = new ArrayList<>(countColumnNames.size() + 1);
        columnNames.add(PON_SAMPLE_BETA_HAT_COLUMN_NAME);
        columnNames.addAll(countColumnNames);
        final TableColumnCollection columns = new TableColumnCollection(columnNames);
        try (final TableWriter<Integer> writer = TableUtils.writer(betaHatsOutFile, columns,
                (i,dataLine) -> dataLine.append(Integer.toString(i)).append(tangentBetaHats.getRow(i)) )) {
            writer.writeComment("fileFormat  = tsv");
            writer.writeComment("commandLine = " + getCommandLine());
            writer.writeComment("title = Tangent normalization Beta Hats");
            for (int i = 0; i < tangentBetaHats.getRowDimension(); i++) {
                writer.writeRecord(i);
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(betaHatsOutFile,ex.getMessage());
        }
    }

    /**
     * Prepares the data to perform tangent factor normalization.
     * <p>
     * This is done count group or column:
     *   <ol>
     *     </li>we change counts to their ratio versus the column mean,</li>
     *     </li>then we transform value to their log_2,</li>
     *     </li>and finally we center then around the median.</li>
     *   </ol>
     * </p>
     *
     * @param matrix input matrix.
     * @return never {@code null}.
     */
    private RealMatrix composeTangentNormalizationInputMatrix(final RealMatrix matrix) {

        final double[] columnInverseMeans = calculateColumnInverseMeans(matrix);

        final double[][] result = log2ColumnMeanCenteredCounts(matrix, columnInverseMeans);

        centerAroundColumnMedian(result);

        return new Array2DRowRealMatrix(result,false);
    }

    /**
     * Center the values around the median per column.
     * @param result the input and output matrix where the first dimension are rows and the second columns.
     */
    private void centerAroundColumnMedian(final double[][] result) {
        for (int i = 0; i < result[0].length; i++) {
            final DescriptiveStatistics stats = new DescriptiveStatistics();
            for (final double[] aResult : result) {
                stats.addValue(aResult[i]);
            }
            final double median = stats.getPercentile(50);
            for (int j = 0; j < result.length; j++) {
                result[j][i] -= median;
            }
        }
    }

    private double[][] log2ColumnMeanCenteredCounts(RealMatrix matrix, double[] columnInverseMeans) {
        final double[][] result = new double[matrix.getRowDimension()][columnInverseMeans.length];
        // First we divide by column mean then becoming a ratio;
        // We also log_2 transform it at this point.
        for (int i = 0; i < result.length; i++) {
            final double[] inputRow = matrix.getRow(i);
            for (int j = 0; j < columnInverseMeans.length; j++) {
                if ((result[i][j] = inputRow[j] * columnInverseMeans[j]) < EPSILON) {
                    result[i][j] = EPSILON;
                }
                result[i][j] = Math.log(result[i][j]) * INV_LN2;
            }
        }
        return result;
    }

    private double[] calculateColumnInverseMeans(RealMatrix matrix) {
        return IntStream.range(0, matrix.getColumnDimension())
                    .mapToDouble(i ->
                                    1.0 / IntStream.range(0, matrix.getRowDimension())
                                            .mapToDouble(j -> matrix.getEntry(j, i))
                                            .average().orElseThrow(
                                                    () -> new IllegalArgumentException("cannot calculate a average for column " + i)
                                            )
                    ).toArray();
    }

    /**
     * Reads the exon collection from a file.
     * @param targetFile the input target file.
     * @return never {@code null}.
     */
    private TargetCollection<? extends BEDFeature> readTargetCollection(final File targetFile) {
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
                final TargetCollection<? extends BEDFeature> result = TargetCollections.fromBEDFeatureFile(targetFile, bedCodec);
                logger.log(Level.INFO, String.format("Found %d targets to analyze.", result.targetCount()));
                return result;
            } else {
                throw new UserException.BadInput(String.format("currently only BED formatted target file are supported. '%s' does not seem to be a BED file",
                        targetFile.getAbsolutePath()));
            }
        }
    }

    /**
     * Reads the read-counts from the input file using an exon collection (if provided) to resolve target names if this
     * are missing.
     *
     * @param readCountsFile the read-counts file.
     * @param targetCollection the input exon collection. {@code null} indicates that no collection was provided by the user.
     * @param <E>            the element type of {@code exonCollection}, it must be a {@link BEDFeature}.
     * @return never {@code null}.
     * @throws UserException.CouldNotReadInputFile if there was some problem
     *                                             trying to read from {@code readCountsFile}.
     * @throws UserException.BadInput              if there is some format issue with the input file {@code readCountsFile},
     *                                             or it does not contain target names and {@code exonCollection} is {@code null}.
     */
    private <E extends BEDFeature> ReadCountCollection readInputReadCounts(final File readCountsFile,
                                                                           final TargetCollection<E> targetCollection) {
        try {
            return ReadCountCollectionUtils.parse(readCountsFile, targetCollection, false);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(readCountsFile, ex.getMessage(), ex);
        }
    }

    /**
     * Writes the target-factor normalized outputs.
     * @param targetFactorNormalized the read count collection to write.
     */
    private void writeTargetFactorNormalizedOutput(final ReadCountCollection targetFactorNormalized) {
        writeOutput(fntOutFile, targetFactorNormalized, "Target factor normalized target counts");
    }

    /**
     * Writes the tangent-factor normalized outputs.
     *
     * Please note that this file is written in the target format originally used in recapseg.  In other words, it is
     *  written as if it was target coverage.
     *
     * @param tangentNormalized the read count collection to write.
     */
    private void writeTangentNormalizedOutput(final ReadCountCollection tangentNormalized) {
        try {
            ReadCountCollectionUtils.writeAsTargetCoverage(outFile, tangentNormalized, "fileFormat = tsv",
                    "commandLine = " + getCommandLine(),
                    "title = Tangent normalized coverage profile");
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outFile, ex.getMessage());
        }
    }


    private void writeOutput(final File file, final ReadCountCollection counts, final String title) {
        if (file == null) {
            return;
        }
        try {
            ReadCountCollectionUtils.write(file, counts,
                    "fileFormat = tsv",
                    "commandLine = " + getCommandLine(),
                    "title = " + title);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, ex.getMessage());
        }
    }


    /**
     * Utility class use to map matrices rows between two different target lists.
     * <p>
     *     The "case" target list is provided by a {@link ReadCountCollection}, whereas the "pon" target list
     *     is provided by a {@link PoN}.
     * </p>
     * <p>
     *     This class contains method to re-organized a matrix rows that is based on one list to match the order
     *     impose by the other list.
     * </p>
     * <p>
     *     It also provided methods to create a PoN re-organized count matrix directly from a {@link ReadCountCollection}.
     *     object and the reverse operation, that is to create a new {@link ReadCountCollection} from a count matrix
     *     whose rows have been reorganized to satisfy the {@link PoN} target order.
     * </p>
     */
    private static final class Case2PoNTargetMapper {

        private final List<String> ponTargetNames;
        private final Map<String, Integer> caseTargetIndexes;
        private final List<Target> outputTargets;
        private final File readCountsFile;

        /**
         * Creates a new mapper.
         * @param caseTargets the case sample targets as they appear in the input read counts.
         * @param ponTargetNames the PoN target names in the order they are present in the PoN.
         */
        private Case2PoNTargetMapper(final List<Target> caseTargets, final List<String> ponTargetNames, final File readCountsFile)  {
            this.readCountsFile = readCountsFile;
            this.ponTargetNames = ponTargetNames;
            this.caseTargetIndexes = IntStream.range(0, caseTargets.size()).boxed()
                    .collect(Collectors.toMap(i -> caseTargets.get(i).getName(), Function.identity()));
            if (this.ponTargetNames.stream().anyMatch(n -> ! caseTargetIndexes.containsKey(n))) {
                throw new UserException.BadInput(
                        String.format("the read-count input file '%s' is missing some target in the PoN: e.g. '%s'",
                                readCountsFile,this.ponTargetNames.stream().filter(n -> ! caseTargetIndexes.containsKey(n)).limit(10).collect(Collectors.joining(", "))));
            }
            this.outputTargets =  ponTargetNames.stream()
                    .map(caseTargetIndexes::get).sorted().map(caseTargets::get).collect(Collectors.toList());
        }

        /**
         * Re-arrange case read counts counts based on PoN target order.
         * @param caseCounts the input counts to rearrange.
         * @return never {@code null}, a new matrix with row sorted according to the PoN target order.
         */
        private RealMatrix fromCaseToPoNCounts(final RealMatrix caseCounts) {
            final RealMatrix result = new Array2DRowRealMatrix(ponTargetNames.size(), caseCounts.getColumnDimension());
            for (int i = 0; i < ponTargetNames.size(); i++) {
                final String targetName = ponTargetNames.get(i);
                final Integer inputIndex = caseTargetIndexes.get(targetName);
                if (inputIndex == null) {
                    throw new UserException.BadInput(String.format("missing PoN target name %s in input read counts %s", targetName, readCountsFile));
                }
                result.setRow(i, caseCounts.getRow(inputIndex));
            }
            return result;
        }

        /**
         * Re-arrange the input rows from the PoN to the case data target order.
         * @param ponCounts count matrix with row organized using the PoN target order.
         * @return never {@code null} a new matrix with the row order changed according to the case read count target order.
         */
        private RealMatrix fromPoNToCaseCounts(final RealMatrix ponCounts) {
            final Map<String,Integer> ponTargetsIndexes = IntStream.range(0, ponTargetNames.size()).boxed()
                    .collect(Collectors.toMap(ponTargetNames::get, Function.identity()));
            final RealMatrix result = new Array2DRowRealMatrix(ponCounts.getRowDimension(),ponCounts.getColumnDimension());
            for (int i = 0; i < outputTargets.size(); i++) {
                final Target target = outputTargets.get(i);
                result.setRow(i, ponCounts.getRow(ponTargetsIndexes.get(target.getName())));
            }
            return result;
        }

        /**
         * Given a read counts using the PoN target order to a read-count collection
         * based on the case sample target order.
         *
         * @param counts PoN target sorted count matrix.
         * @param countColumnNames count column names.
         * @return never {@code null} a read-count collection with the row order according to the case read count target order.
         */
        private ReadCountCollection fromPoNtoCaseCountCollection(final RealMatrix counts, final List<String> countColumnNames) {
            return new ReadCountCollection(SetUniqueList.setUniqueList(new ArrayList<>(outputTargets)),SetUniqueList.setUniqueList(new ArrayList<>(countColumnNames)),
                    fromPoNToCaseCounts(counts));
        }
    }

}
