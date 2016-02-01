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

    public enum UndoSplits {

        // Default is NONE
        NONE("none"), PRUNE("prune"), SDUNDO("sdundo");

        private final String value;

        UndoSplits(final String value) {
            this.value = value;
        }

        @Override
        public String toString() {
            return value;
        }
    }

    public enum PMethod {

        // Default is NONE
        HYBRID("hybrid"), PERM("perm");

        private final String value;

        PMethod(final String value) {
            this.value = value;
        }

        @Override
        public String toString() {
            return value;
        }
    }

    private RCBSSegmenter() {
    }

    /**
     * Create a segmentation file using CBS in R.
     *
     * <p>https://www.bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf</p>
     *
     * <p>Please see the above documentation for a more detailed description of the parameters.</p>
     *
     * <p>IMPORTANT:  There is no check that the weights have the same number of entries as the tnFile.</p>
     *
     * Wraps the call:
     * segment(x, weights = NULL, alpha = 0.01, nperm = 10000, p.method =
     *  c("hybrid", "perm"), min.width=2, kmax=25, nmin=200,
     *  eta=0.05, sbdry=NULL, trim = 0.025, undo.splits =
     *  c("none", "prune", "sdundo"), undo.prune=0.05,
     *  undo.SD=3, verbose=1)
     *
     * <p> Note that sbdry and verbose are unavailable through this interface. </p>
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
    public static void writeSegmentFile(final String sampleName, final String tnFile, final String outputFile,
                                        final Boolean log, final File weightFile, final double alpha,
                                        final int nperm, final PMethod pmethod, final int minWidth,
                                        final int kmax, final int nmin, final double eta, final double trim,
                                        final UndoSplits undoSplits, final double undoPrune, final int undoSD
                                        ) {

        String logArg = log ? "TRUE" : "FALSE";
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name="+sampleName, "--targets_file="+tnFile, "--output_file="+outputFile,
                "--log2_input="+logArg, "--min_width="+String.valueOf(minWidth),
                "--alpha=" + String.valueOf(alpha), "--nperm=" + String.valueOf(nperm),
                "--pmethod=" + pmethod.toString(), "--kmax=" + String.valueOf(kmax),
                "--nmin=" + String.valueOf(nmin), "--eta=" + String.valueOf(eta),
                "--trim=" + String.valueOf(trim), "--undosplits=" + undoSplits.toString(),
                "--undoprune=" + String.valueOf(undoPrune), "--undoSD=" + String.valueOf(undoSD));

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

    /**
     *  Write segment file with default parameters
     *
     * @param sampleName Name of the sample being run through the segmenter.  Never {@code null}
     * @param tnFile Tangent-normalized targets file.  Never {@code null}
     * @param outputFile Full path to the outputted segment file.  Never {@code null}
     * @param log whether the tnFile input has already been put into log2CR.  Never {@code null}
     */
    public static void writeSegmentFile(final String sampleName, final String tnFile, final String outputFile, final Boolean log) {
        writeSegmentFile(sampleName, tnFile, outputFile, log, null);
    }

    public static void writeSegmentFile(final String sampleName, final String tnFile, final String outputFile, final Boolean log, final File weightsFile) {
        writeSegmentFile(sampleName, tnFile, outputFile, log, weightsFile, 0.01, 10000, PMethod.HYBRID, 2, 25, 200, 0.05,
                0.025, UndoSplits.NONE, 0.05, 3);
    }
}