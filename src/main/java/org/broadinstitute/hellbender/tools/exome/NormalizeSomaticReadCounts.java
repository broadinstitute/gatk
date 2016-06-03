package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.Level;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
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

/**
 * Normalizes read counts given the PanelOfNormals.
 *
 * <p> Dev note:  If this is extended to use spark, please be wary that the parallelization in tangent normalization is
 * by case sample, which may not yield benefits for most use cases (which are one sample)  </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Normalizes PCOV read counts using a panel of normals",
        oneLineSummary = "Normalizes proportional coverage (PCOV) read counts using a panel of normals",
        programGroup = CopyNumberProgramGroup.class
)
public final class NormalizeSomaticReadCounts extends CommandLineProgram {

    /**
     * Name of the column that contains the I.D. of the PoN <i>"eigen sample"</i>.
     */
    public static final String PON_SAMPLE_BETA_HAT_COLUMN_NAME = "PON.SAMPLE";

    public static final String READ_COUNTS_FILE_FULL_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String FACTOR_NORMALIZED_COUNTS_LONG_NAME = "factorNormalizedOutput";
    public static final String FACTOR_NORMALIZED_COUNTS_SHORT_NAME = "FNO";

    public static final String TANGENT_BETA_HATS_LONG_NAME = "betaHatsOutput";
    public static final String TANGENT_BETA_HATS_SHORT_NAME = "BHO";

    @Argument(
            doc = "read counts input file.  This can only contain one sample.",
            shortName = READ_COUNTS_FILE_SHORT_NAME,
            fullName = READ_COUNTS_FILE_FULL_NAME,
            optional = false
    )
    protected File readCountsFile;

    @Argument(
            doc = "target BED file.",
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            optional = true
    )
    protected File targetFile;

    @Argument(
            doc = "panel of normals HDF5 file",
            shortName = ExomeStandardArgumentDefinitions.PON_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.PON_FILE_LONG_NAME,
            optional = false
    )
    protected File ponFile;

    @Argument(
            doc = "Factor normalized counts output",
            shortName = FACTOR_NORMALIZED_COUNTS_SHORT_NAME,
            fullName = FACTOR_NORMALIZED_COUNTS_LONG_NAME,
            optional = true
    )
    protected File fntOutFile;

    @Argument(
            doc = "Tangent normalized counts output",
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected File outFile;

    @Argument(
            doc = "Tangent normalization Beta Hats output file",
            shortName = TANGENT_BETA_HATS_SHORT_NAME,
            fullName = TANGENT_BETA_HATS_LONG_NAME,
            optional = true
    )
    protected File betaHatsOutFile;

    @Argument(
            doc = "Pre-tangent normalization counts",
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
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
            final TangentNormalizationResult tangentNormalizationResult = TangentNormalizer.tangentNormalizePcov(pon, readCountCollection);
            outputTangentNormalizationResult(tangentNormalizationResult);
            return "SUCCESS";
        }
    }

    private void outputTangentNormalizationResult(final TangentNormalizationResult tangentNormalizationResult){
        writeTargetFactorNormalizedOutput(tangentNormalizationResult.getTargetFactorNormalizedCounts());
        writePreTangentNormalizationOutput(tangentNormalizationResult.getPreTangentNormalized());
        writeTangentBetaHats(tangentNormalizationResult.getTangentBetaHats(), tangentNormalizationResult.getTargetFactorNormalizedCounts().columnNames());
        writeTangentNormalizedOutput(tangentNormalizationResult.getTangentNormalized());
    }

    /**
     * Writes the pre-tangent-normalization read counts if a file was provided for it.
     *
     * This file is written in the target format originally used in recapseg.  In other words, it is
     *  written as if it was target coverage. If the output file is null, we do not write a file.
     *
     * @param preTangentNormalized the read count collection to write.
     */
    private void writePreTangentNormalizationOutput(final ReadCountCollection preTangentNormalized) {
        // If the output file is null, we do not write a file.
        if (preTangentNormalizationOutFile == null) {
            return;
        }

        try {
            ReadCountCollectionUtils.write(preTangentNormalizationOutFile, preTangentNormalized, "fileFormat = tsv",
                    "commandLine = " + getCommandLine(), "title = Pre tangent normalized coverage profile");
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(preTangentNormalizationOutFile, ex.getMessage());
        }
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
                (i, dataLine) -> dataLine.append(Integer.toString(i)).append(tangentBetaHats.getRow(i)) )) {
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
                final TargetCollection<? extends BEDFeature> result = TargetCollectionUtils.fromBEDFeatureFile(targetFile, bedCodec);
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
     * Even though the tangent normalization code works for multiple samples, our workflows are designed for single samples
     * and so we enforce it here.
     *
     * @param readCountsFile the read-counts file.
     * @param targetCollection the input exon collection. {@code null} indicates that no collection was provided by the user.
     * @param <E>            the element type of {@code exonCollection}, it must be a {@link BEDFeature}.
     * @return never {@code null}.
     * @throws UserException.CouldNotReadInputFile if there was some problem
     *                                             trying to read from {@code inputReadCountsFile}.
     * @throws UserException.BadInput              if there is some format issue with the input file {@code inputReadCountsFile},
     *                                             or it does not contain target names and {@code exonCollection} is {@code null}
     *                                             or the input read counts contain more thanb one sample.
     */
    private <E extends BEDFeature> ReadCountCollection readInputReadCounts(final File readCountsFile,
                                                                           final TargetCollection<E> targetCollection) {
        try {
            final ReadCountCollection result = ReadCountCollectionUtils.parse(readCountsFile, targetCollection, false);
            if (result.columnNames().size() > 1) {
                throw new UserException.BadInput("Only single-sample input read counts are allowed");
            }
            return result;
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
            ReadCountCollectionUtils.write(outFile, tangentNormalized, "fileFormat = tsv",
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
}
