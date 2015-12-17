package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

/**
 * Calls an R script to perform segmentation
 */
public final class RCBSSegmenter {
    private static final String R_SCRIPT = "CBS.R";

    private RCBSSegmenter() {
    }

    /**
     * @param sample_name Name of the sample being run through the segmenter
     * @param tnFile Tangent-normalized targets file
     * @param outputFile Full path to the outputted segment file
     */
    public static void writeSegmentFile(final String sample_name, final String tnFile, final String outputFile, final Boolean log) {
        final String logArg = log ? "TRUE" : "FALSE";
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name="+sample_name, "--targets_file="+tnFile, "--output_file="+outputFile,
                "--log2_input="+logArg);
        executor.exec();
    }
}