package org.broadinstitute.hellbender.tools.utilities;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;

/**
 * Command line tool for analyzing and displaying metadata about a SAM/BAM/CRAM/VCF/BCF file or index file.
 */
public abstract class AnalyzeFileBase extends CommandLineProgram {

    @Argument(shortName="target-path", fullName="target-path",
            doc="Path to file to be analyzed",
            optional=false)
    private String targetPath;

    @Argument(shortName="companionIndex", fullName="companionIndex",
            doc="Path of companion index file to be analyzed",
            optional=true)
    private String companionIndex;

    @Argument(shortName="analysisDepth", fullName="analysisDepth",
            doc="depth of analysis to be performed")
    private int analysisDepth;

    protected String fileName;

    public AnalyzeFileBase(String fileName) {
        this.fileName = fileName;
    }

}
