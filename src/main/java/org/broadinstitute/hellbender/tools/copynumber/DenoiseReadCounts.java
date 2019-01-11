package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.denoising.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;

/**
 * Denoises read counts to produce denoised copy ratios.
 *
 * <p>
 *     Typically, a panel of normals produced by {@link CreateReadCountPanelOfNormals} is provided as input.
 *     The input counts are then standardized by 1) transforming to fractional coverage,
 *     2) performing optional explicit GC-bias correction (if the panel contains GC-content annotated intervals),
 *     3) filtering intervals to those contained in the panel, 4) dividing by interval medians contained in the panel,
 *     5) dividing by the sample median, and 6) transforming to log2 copy ratio.  The result is then denoised by
 *     subtracting the projection onto the specified number of principal components from the panel.
 * </p>
 *
 * <p>
 *     If no panel is provided, then the input counts are instead standardized by 1) transforming to fractional coverage,
 *     2) performing optional explicit GC-bias correction (if GC-content annotated intervals are provided),
 *     3) dividing by the sample median, and 4) transforming to log2 copy ratio.  No denoising is performed,
 *     so the denoised result is simply taken to be identical to the standardized result.
 * </p>
 *
 * <p>
 *     If performed, explicit GC-bias correction is done by {@link GCBiasCorrector}.
 * </p>
 *
 * <p>
 *     Note that {@code number-of-eigensamples} principal components from the input panel will be used for
 *     denoising; if only fewer are available in the panel, then they will all be used.  This parameter can
 *     thus be used to control the amount of denoising, which will ultimately affect the sensitivity of the analysis.
 * </p>
 *
 * <p>
 *     See comments for {@link CreateReadCountPanelOfNormals} regarding coverage on sex chromosomes.  If sex
 *     chromosomes are not excluded from coverage collection, it is strongly recommended that case samples are
 *     denoised only with panels containing only individuals of the same sex as the case samples.
 * </p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Counts TSV or HDF5 file from {@link CollectReadCounts}.
 *     </li>
 *     <li>
 *         (Optional) Panel-of-normals from {@link CreateReadCountPanelOfNormals}.
 *         If provided, it will be used to standardize and denoise the input counts.  This may include
 *         explicit GC-bias correction if annotated intervals were used to create the panel.
 *     </li>
 *     <li>
 *         (Optional) GC-content annotated-intervals from {@link AnnotateIntervals}.
 *         This can be provided in place of a panel of normals to perform explicit GC-bias correction.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         Standardized-copy-ratios file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CopyRatioCollection.CopyRatioTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 *     <li>
 *         Denoised-copy-ratios file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CopyRatioCollection.CopyRatioTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 *     gatk DenoiseReadCounts \
 *          -I sample.counts.hdf5 \
 *          --count-panel-of-normals panel_of_normals.pon.hdf5 \
 *          --standardized-copy-ratios sample.standardizedCR.tsv \
 *          --denoised-copy-ratios sample.denoisedCR.tsv
 * </pre>
 *
 * <pre>
 *     gatk DenoiseReadCounts \
 *          -I sample.counts.hdf5 \
 *          --annotated-intervals annotated_intervals.tsv \
 *          --standardized-copy-ratios sample.standardizedCR.tsv \
 *          --denoised-copy-ratios sample.denoisedCR.tsv
 * </pre>
 *
 * <pre>
 *     gatk DenoiseReadCounts \
 *          -I sample.counts.hdf5 \
 *          --standardized-copy-ratios sample.standardizedCR.tsv \
 *          --denoised-copy-ratios sample.denoisedCR.tsv
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Denoises read counts to produce denoised copy ratios",
        oneLineSummary = "Denoises read counts to produce denoised copy ratios",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class DenoiseReadCounts extends CommandLineProgram {
    @Argument(
            doc = "Input TSV or HDF5 file containing integer read counts in genomic intervals for a single case sample (output of CollectReadCounts).",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private File inputReadCountFile;

    @Argument(
            doc = "Input HDF5 file containing the panel of normals (output of CreateReadCountPanelOfNormals).",
            fullName = CopyNumberStandardArgument.COUNT_PANEL_OF_NORMALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputPanelOfNormalsFile = null;

    @Argument(
            doc = "Input file containing annotations for GC content in genomic intervals (output of AnnotateIntervals).  " +
                    "Intervals must be identical to and in the same order as those in the input read-counts file.  " +
                    "If a panel of normals is provided, this input will be ignored.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAnnotatedIntervalsFile = null;

    @Argument(
            doc = "Output file for standardized copy ratios.  GC-bias correction will be performed if annotations for GC content are provided.",
            fullName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_LONG_NAME
    )
    private File outputStandardizedCopyRatiosFile;

    @Argument(
            doc = "Output file for denoised copy ratios.",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME
    )
    private File outputDenoisedCopyRatiosFile;

    @Argument(
            doc = "Number of eigensamples to use for denoising.  " +
                    "If not specified or if the number of eigensamples available in the panel of normals " +
                    "is smaller than this, all eigensamples will be used.",
            fullName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private Integer numEigensamplesRequested = null;

    @Override
    protected Object doWork() {
        if (!new HDF5Library().load(null)) { //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        IOUtils.canReadFile(inputReadCountFile);
        logger.info(String.format("Reading read-counts file (%s)...", inputReadCountFile));
        final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);

        if (inputPanelOfNormalsFile != null) {  //denoise using panel of normals
            IOUtils.canReadFile(inputPanelOfNormalsFile);
            try (final HDF5File hdf5PanelOfNormalsFile = new HDF5File(inputPanelOfNormalsFile)) {  //HDF5File implements AutoCloseable
                final SVDReadCountPanelOfNormals panelOfNormals = HDF5SVDReadCountPanelOfNormals.read(hdf5PanelOfNormalsFile);

                if (inputAnnotatedIntervalsFile != null) {
                    logger.warn("Panel of normals was provided; ignoring input GC-content annotations...");
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

                logger.info("Writing standardized and denoised copy ratios...");
                denoisedCopyRatioResult.write(outputStandardizedCopyRatiosFile, outputDenoisedCopyRatiosFile);
            }
        } else {    //standardize and perform optional GC-bias correction
            //get GC content (null if not provided)
            final AnnotatedIntervalCollection annotatedIntervals = CopyNumberArgumentValidationUtils.validateAnnotatedIntervals(
                    inputAnnotatedIntervalsFile, readCounts, logger);
            final double[] intervalGCContent = annotatedIntervals == null
                    ? null
                    : annotatedIntervals.getRecords().stream()
                    .mapToDouble(i -> i.getAnnotationMap().getValue(CopyNumberAnnotations.GC_CONTENT))
                    .toArray();

            if (intervalGCContent == null) {
                logger.warn("Neither a panel of normals nor GC-content annotations were provided, so only standardization will be performed...");
            }

            final RealMatrix standardizedCopyRatioValues = SVDDenoisingUtils.preprocessAndStandardizeSample(readCounts.getCounts(), intervalGCContent);

            //construct a result with denoised result identical to standardized result
            final SVDDenoisedCopyRatioResult standardizedResult = new SVDDenoisedCopyRatioResult(
                    readCounts.getMetadata(),
                    readCounts.getIntervals(),
                    standardizedCopyRatioValues,
                    standardizedCopyRatioValues);
            standardizedResult.write(outputStandardizedCopyRatiosFile, outputDenoisedCopyRatiosFile);
        }

        logger.info("Read counts successfully denoised.");

        return "SUCCESS";
    }
}
