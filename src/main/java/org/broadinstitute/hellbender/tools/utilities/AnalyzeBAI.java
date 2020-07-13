package org.broadinstitute.hellbender.tools.utilities;

import htsjdk.samtools.BAMIndexer;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;

import java.io.File;

/**
 * Analyzer for BAI files.
 */
@ExperimentalFeature
@CommandLineProgramProperties(
        summary = "Print information about a .bai index",
        oneLineSummary = "Print information about a .bai index",
        programGroup = OtherProgramGroup.class
)
public class AnalyzeBAI extends CommandLineProgram {
    @Argument(
            shortName=StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName=StandardArgumentDefinitions.INPUT_LONG_NAME,
            doc="File to be analyzed",
            optional=false)
    private File inputPath; // File due to htsjdk signature

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "File to which to write file analysis information (if not specified, prints to standard output)",
            optional = true)
    private File outputPath;

    /**
     * Run the analyzer for the file.
     * @return
     */
    protected Object doWork() {
        System.out.println(String.format("\nOutput written to %s\n", outputPath));
        BAMIndexer.createAndWriteIndex(inputPath, outputPath, true);

        return 1;
    }

}

