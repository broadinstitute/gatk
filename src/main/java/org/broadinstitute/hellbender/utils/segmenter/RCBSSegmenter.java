package org.broadinstitute.hellbender.utils.segmenter;

import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;
import java.io.File;

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
     * @param minLogValue Any values under this threshold will be set to it
     * @return File created by segmenter
     */
    public static File segment(String sample_name, String tnFile, String outputFile, Float minLogValue) {
        final RScriptExecutor exectutor = new RScriptExecutor();
        exectutor.addScript(new Resource(R_SCRIPT, RCBSSegmenter.class));
        exectutor.addArgs(sample_name, tnFile, outputFile, minLogValue);
        exectutor.exec();
        return new File(outputFile);
    }
}