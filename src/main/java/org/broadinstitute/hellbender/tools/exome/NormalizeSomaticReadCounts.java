package org.broadinstitute.hellbender.tools.exome;

import org.apache.logging.log4j.Level;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.HDF5PCACoveragePoN;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.PCACoveragePoN;
import org.broadinstitute.hellbender.tools.pon.coverage.pca.PCATangentNormalizationResult;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Normalizes read counts given the PanelOfNormals (PoN).
 *
 * <p> A note to developers:  If this is extended to use Spark, please be wary that the parallelization in tangent normalization is
 * by case sample, which may not yield benefits for most use cases (which are one sample)  </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 *
 * <h3>Examples</h3>
 * <p>
 *     The command encompasses empirically determined parameters for TCGA project data.
 *     You may obtain better results with different parameters.
 * </p>
 *
 * <p>
 *     The following command is for either whole exome sequencing (WES) or whole genome sequencing (WGS) data.
 * </p>
 *
 * <pre>
 * java -Xmx4g -jar $gatk_jar NormalizeSomaticReadCounts \
 *   --input coverage.tsv \
 *   --panelOfNormals cnv_panel_of_normals.pon \
 *   --tangentNormalized entity_id.tn.tsv \
 *   --factorNormalizedOutput entity_id.fnt.tsv \
 *   --preTangentNormalized entity_id.preTN.tsv \
 *   --betaHatsOutput entity_id.betaHats.tsv
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Normalize PCOV read counts using a panel of normals",
        oneLineSummary = "Normalize proportional coverage (PCOV) read counts using a panel of normals",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class NormalizeSomaticReadCounts extends CommandLineProgram {

    public static final String READ_COUNTS_FILE_FULL_NAME = StandardArgumentDefinitions.INPUT_LONG_NAME;
    public static final String READ_COUNTS_FILE_SHORT_NAME = StandardArgumentDefinitions.INPUT_SHORT_NAME;

    public static final String TANGENT_BETA_HATS_LONG_NAME = "betaHatsOutput";
    public static final String TANGENT_BETA_HATS_SHORT_NAME = "BHO";

    public static final String FACTOR_NORMALIZED_COUNTS_LONG_NAME = "factorNormalizedOutput";
    public static final String FACTOR_NORMALIZED_COUNTS_SHORT_NAME = "FNO";

    @Argument(
            doc = "read counts input file.  This can only contain one sample.",
            shortName = READ_COUNTS_FILE_SHORT_NAME,
            fullName = READ_COUNTS_FILE_FULL_NAME,
            optional = false
    )
    protected File readCountsFile;

    @Argument(
            doc = "target file -- not a BED file.  Should be formatted as a tsv with at least the following header columns: contig, start, stop, name.",
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
            doc = "Tangent normalized counts output",
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected File tangentNormalizationOutFile;

    @Argument(
            doc = "Pre-tangent normalization counts",
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    protected File preTangentNormalizationOutFile;

    @Argument(
            doc = "Tangent normalization Beta Hats output file",
            shortName = TANGENT_BETA_HATS_SHORT_NAME,
            fullName = TANGENT_BETA_HATS_LONG_NAME,
            optional = true
    )
    protected File betaHatsOutFile;

    @Argument(
            doc = "Factor normalized counts output",
            shortName = FACTOR_NORMALIZED_COUNTS_SHORT_NAME,
            fullName = FACTOR_NORMALIZED_COUNTS_LONG_NAME,
            optional = true
    )
    protected File fntOutFile;

    @Override
    protected Object doWork() {
        if (! new HDF5Library().load(null)){ //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }
        IOUtils.canReadFile(ponFile);
        try (final HDF5File hdf5PoNFile = new HDF5File(ponFile)) {
            final PCACoveragePoN pon = new HDF5PCACoveragePoN(hdf5PoNFile, logger);
            final TargetCollection<Target> targetCollection = readTargetCollection(targetFile);
            final ReadCountCollection proportionalCoverageProfile = readInputReadCounts(readCountsFile, targetCollection);
            final PCATangentNormalizationResult tangentNormalizationResult = pon.normalize(proportionalCoverageProfile);;
            tangentNormalizationResult.write(getCommandLine(), tangentNormalizationOutFile, preTangentNormalizationOutFile, betaHatsOutFile, fntOutFile);
            return "SUCCESS";
        }
    }

    /**
     * Reads the target collection from a file.
     * @param targetFile the input target file.
     * @return never {@code null}.
     */
    private TargetCollection<Target> readTargetCollection(final File targetFile) {
        if (targetFile == null) {
            return null;
        } else {
            IOUtils.canReadFile(targetFile);
            logger.log(Level.INFO, String.format("Reading target intervals from exome file '%s' ...", new Object[]{targetFile.getAbsolutePath()}));
            final List<Target> targets = TargetTableReader.readTargetFile(targetFile);
            return new HashedListTargetCollection<>(targets);
        }
    }

    /**
     * Reads the read-counts from the input file using a target collection (if provided) to resolve target names if this
     * are missing.
     *
     * Even though the tangent normalization code works for multiple samples, our workflows are designed for single samples
     * and so we enforce it here.
     *
     * @param readCountsFile the read-counts file.
     * @param targetCollection the input target collection. {@code null} indicates that no collection was provided by the user.
     * @return never {@code null}.
     * @throws UserException.CouldNotReadInputFile if there was some problem
     *                                             trying to read from {@code inputReadCountsFile}.
     * @throws UserException.BadInput              if there is some format issue with the input file {@code inputReadCountsFile},
     *                                             or it does not contain target names and {@code targetCollection} is {@code null}
     *                                             or the input read counts contain more thanb one sample.
     */
    private ReadCountCollection readInputReadCounts(final File readCountsFile,
                                                    final TargetCollection<Target> targetCollection) {
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
}
