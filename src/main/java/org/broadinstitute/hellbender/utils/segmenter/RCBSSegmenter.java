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
    public static void writeSegmentFile(String sample_name, String tnFile, String outputFile, Boolean log) {
        String logArg = "FALSE";
        if(log){
            logArg = "TRUE";
        }
        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        /*--args is needed for Rscript to recognize other arguments properly*/
        executor.addArgs("--args", "--sample_name="+sample_name, "--targets_file="+tnFile, "--output_file="+outputFile,
                "--log2_input="+logArg);
        executor.exec();
    }
}