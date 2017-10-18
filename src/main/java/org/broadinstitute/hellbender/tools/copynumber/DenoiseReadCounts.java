package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.annotation.GCBiasCorrector;
import org.broadinstitute.hellbender.tools.copynumber.coverage.denoising.svd.HDF5SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.coverage.denoising.svd.SVDDenoisedCopyRatioResult;
import org.broadinstitute.hellbender.tools.copynumber.coverage.denoising.svd.SVDDenoisingUtils;
import org.broadinstitute.hellbender.tools.copynumber.coverage.denoising.svd.SVDReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;

/**
 * Denoises read counts given the panel of normals (PoN) created by {@link CreateReadCountPanelOfNormals} to produce
 * a copy-ratio profile.  Both HDF5 and TSV read-count file input are supported.
 *
 * <h3>Examples</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" DenoiseReadCounts \
 *   --input tumor.readCounts.tsv \
 *   --readCountPanelOfNormals panel_of_normals.hdf5 \
 *   --standardizedCopyRatios tumor.standardizedCR.tsv \
 *   --denoisedCopyRatios tumor.denoisedCR.tsv
 * </pre>
 *
 * The resulting copy-ratio profile is log2 transformed.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Denoise read counts using a panel of normals.",
        oneLineSummary = "Denoise read counts using a panel of normals.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class DenoiseReadCounts extends CommandLineProgram {
    @Argument(
            doc = "Input TSV or HDF5 read-count file containing integer read counts in genomic intervals for a single case sample.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private File inputReadCountFile;

    @Argument(
            doc = "Input HDF5 file containing the panel of normals (output of CreateReadCountPanelOfNormals).",
            fullName = CopyNumberStandardArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.READ_COUNT_PANEL_OF_NORMALS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputPanelOfNormalsFile = null;

    @Argument(
            doc = "Input annotated-interval file containing annotations for GC content in genomic intervals (output of AnnotateIntervals).  " +
                    "Intervals must be identical to and in the same order as those in the input read-count file.  " +
                    "If a panel of normals containing annotations for GC content is provided, this input will be ignored.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME,
            optional = true
    )
    private File annotatedIntervalsFile = null;

    @Argument(
            doc = "Output file for standardized copy-ratio profile.  GC-bias correction will be performed if annotations for GC content are provided.",
            fullName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME
    )
    private File standardizedCopyRatiosFile;

    @Argument(
            doc = "Output file for denoised copy-ratio profile.",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME
    )
    private File denoisedCopyRatiosFile;

    @Argument(
            doc = "Number of eigensamples to use for denoising.  " +
                    "If not specified or if the number of eigensamples available in the panel of normals " +
                    "is smaller than this, all eigensamples will be used.",
            fullName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            shortName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_SHORT_NAME,
            optional = true
    )
    private Integer numEigensamplesRequested = null;

    @Override
    protected Object doWork() {
        if (!new HDF5Library().load(null)) { //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }
        Utils.validateArg(numEigensamplesRequested == null || numEigensamplesRequested > 0,
                "Number of eigensamples to use for denoising must be non-negative.");

        IOUtils.canReadFile(inputReadCountFile);
        logger.info(String.format("Reading read-count file (%s)...", inputReadCountFile));
        final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);

        if (inputPanelOfNormalsFile != null) {  //denoise using panel of normals
            IOUtils.canReadFile(inputPanelOfNormalsFile);
            try (final HDF5File hdf5PanelOfNormalsFile = new HDF5File(inputPanelOfNormalsFile)) {  //HDF5File implements AutoCloseable
                final SVDReadCountPanelOfNormals panelOfNormals = HDF5SVDReadCountPanelOfNormals.read(hdf5PanelOfNormalsFile);

                if (annotatedIntervalsFile != null && panelOfNormals.getOriginalIntervalGCContent() != null) {
                    logger.warn("Panel of normals contains GC-content annotations; ignoring input GC-content annotations...");
                }

                //perform denoising and write result
                final int numEigensamples =
                        numEigensamplesRequested == null ?
                                panelOfNormals.getNumEigensamples() :
                                Math.min(panelOfNormals.getNumEigensamples(), this.numEigensamplesRequested);
                if (numEigensamplesRequested != null && numEigensamples < numEigensamplesRequested) {
                    logger.warn(String.format("%d eigensamples were requested but only %d are available in the panel of normals...",
                            numEigensamplesRequested, numEigensamples));
                }
                final SVDDenoisedCopyRatioResult denoisedCopyRatioResult = panelOfNormals.denoise(readCounts, numEigensamples);

                logger.info("Writing standardized and denoised copy-ratio profiles...");
                denoisedCopyRatioResult.write(standardizedCopyRatiosFile, denoisedCopyRatiosFile);
            }
        } else {    //standardize and perform optional GC-bias correction
            //get GC content (null if not provided)
            final double[] intervalGCContent = GCBiasCorrector.validateIntervalGCContent(readCounts.getIntervals(), annotatedIntervalsFile);

            if (intervalGCContent == null) {
                logger.warn("Neither a panel of normals nor GC-content annotations were provided, so only standardization will be performed...");
            }

            final RealMatrix standardizedCopyRatioValues = SVDDenoisingUtils.preprocessAndStandardizeSample(readCounts.getCounts(), intervalGCContent);

            //construct a result with denoised profile identical to standardized profile
            final SVDDenoisedCopyRatioResult standardizedResult = new SVDDenoisedCopyRatioResult(
                    new SimpleSampleMetadata(readCounts.getSampleName()),
                    readCounts.getIntervals(),
                    standardizedCopyRatioValues,
                    standardizedCopyRatioValues);
            standardizedResult.write(standardizedCopyRatiosFile, denoisedCopyRatiosFile);
        }

        logger.info("Read counts successfully denoised.");

        return "SUCCESS";
    }
}
