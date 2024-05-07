package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.tools.filediagnostics.CRAMIssue8768Analyzer;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.OtherProgramGroup;

/**
 * A diagnostic tool that analyzes a CRAM input file to look for conditions that trigger issue 8768.
 *
 * Writes a report to a file indicating whether the file appears to have read base corruption caused by 8768.
 *
 * By default, the output report will have a summary of the average mismatch rate for all suspected bad containers
 * and a few presumed good containers in order to determine if there is a large difference in the base mismatch rate.
 * To analyze the base mismatch rate for ALL containers, use the "verbose" option.
 *
 * Works on files ending in .cram.
 *
 * Sample Usage:
 *
 * gatk CRAMIssue8768Detector -I input.cram -O output.txt -R reference.fasta
 */
@ExperimentalFeature
@WorkflowProperties
@CommandLineProgramProperties(
        summary = "Analyze a CRAM file for issue 8768",
        oneLineSummary = "Analyze a CRAM file for issue 8768",
        programGroup = OtherProgramGroup.class
)
public class CRAMIssue8768Detector extends CommandLineProgram {
    // default average mismatch rate threshold above which we consider the file to be corrupt
    private static final double DEFAULT_MISMATCH_RATE_THRESHOLD = 0.05;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input path of file to analyze",
            common = true)
    @WorkflowInput
    public GATKPath inputPath;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output diagnostics text file",
            common = true)
    @WorkflowOutput
    public GATKPath textOutputPath;

    public static final String OUTPUT_TSV__ARG_NAME = "output-tsv";
    @Argument(fullName = OUTPUT_TSV__ARG_NAME,
            shortName = OUTPUT_TSV__ARG_NAME,
            doc = "Output diagnostics tsv file",
            optional = true)
    @WorkflowOutput
    public GATKPath tsvOutputPath;

    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
            doc = "Reference for the CRAM file",
            common = true)
    @WorkflowOutput
    public GATKPath referencePath;

    public static final String MISMATCH_RATE_THRESHOLD_ARG_NAME = "mismatch-rate-threshold";
    @Argument(fullName = MISMATCH_RATE_THRESHOLD_ARG_NAME,
            shortName = MISMATCH_RATE_THRESHOLD_ARG_NAME,
            doc = "Mismatch rate threshold above which we consider the file to be corrupt",
            optional = true)
    public double mismatchRateThreshold = DEFAULT_MISMATCH_RATE_THRESHOLD;

    public static final String VERBOSE_ARG_NAME = "verbose";
    @Argument(fullName = VERBOSE_ARG_NAME,
            shortName= VERBOSE_ARG_NAME,
            doc="Calculate and print the mismatch rate for all containers",
            optional=true)
    public boolean verbose = false;

    public static final String ECHO_ARG_NAME = "echo-to-stdout";
    @Argument(fullName = ECHO_ARG_NAME,
            shortName= ECHO_ARG_NAME,
            doc="Echo text output to stdout",
            optional=true)
    public boolean echoToStdout = false;

    private CRAMIssue8768Analyzer cramAnalyzer;

    @Override
    protected Object doWork() {
        cramAnalyzer = new CRAMIssue8768Analyzer(
                inputPath,
                textOutputPath,
                tsvOutputPath,
                referencePath,
                mismatchRateThreshold,
                verbose,
                echoToStdout);
        cramAnalyzer.doAnalysis();
        return cramAnalyzer.getRetCode();
    }

    @Override
    protected void onShutdown() {
        if ( cramAnalyzer != null ) {
            try {
                cramAnalyzer.close();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
    }
}

