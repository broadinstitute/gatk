package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.stream.DoubleStream;

/**
 * Calls an R script to perform segmentation
 */
public final class RCBSSegmenter {

    private static final String R_SCRIPT = "CBS.R";

    private RCBSSegmenter() {
    }

    /**
     * Create a segmentation file using CBS in R.
     *
     * <p>https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf</p>
     *
     * <p>IMPORTANT:  There is no check that the weights have the same number of entries as the tnFile.</p>
     *
     * @param sampleName Name of the sample being run through the segmenter.  Never {@code null}
     * @param tnFile Tangent-normalized targets file.  Never {@code null}
     * @param outputFile Full path to the outputted segment file.  Never {@code null}
     * @param log whether the tnFile input has already been put into log2CR.  Never {@code null}
     * @param minWidth minimum length for a segment
     * @param weightFile File containing weights for each target (doubles; one per line).  Must be the same length as what is in the tnFile (note that this is not
     *                enforced here).  Typically, 1/var(target) is what is used.  All values must be greater than 0.
     *                Use {@code null} if weighting is not desired.  These values should not be log space.
     */
    public static void writeSegmentFile(final String sampleName, final String tnFile, final String outputFile, final Boolean log, final int minWidth, final File weightFile) {

        String logArg = log ? "TRUE" : "FALSE";
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name="+sampleName, "--targets_file="+tnFile, "--output_file="+outputFile,
                "--log2_input="+logArg, "--min_width="+String.valueOf(minWidth));

        if (weightFile != null) {

            final double[] weights = ParamUtils.readValuesFromFile(weightFile);

            // Check to make sure that no weights are zero.
            if (!DoubleStream.of(weights).allMatch(d -> d > 0 && !Double.isNaN(d) && Double.isFinite(d))) {
                throw new GATKException("A weight for a target was zero or less, which is not allowed.  If you truly want zero, you must remove the target from consideration.");
            }

            // Add the argument to specify weights
            executor.addArgs("--weights_file=" + weightFile.getAbsolutePath());
        }

        executor.exec();
    }

    public static void writeSegmentFile(final String sampleName, final String tnFile, final String outputFile, final Boolean log) {
        writeSegmentFile(sampleName, tnFile, outputFile, log, 2, null);
    }
}