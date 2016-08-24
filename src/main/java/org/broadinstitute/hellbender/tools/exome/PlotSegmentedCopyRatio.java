package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.plotter.CopyRatioSegmentedPlotter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Create plots of copy number variant data. Please note that this tool is only supported for hg19 and b37 references. All other references may fail.",
        oneLineSummary = "Create plots of copy number variant data",
        programGroup = CopyNumberProgramGroup.class
)
public final class PlotSegmentedCopyRatio extends CommandLineProgram {

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    @Argument(
            doc = "Genomic targets file after tangent normalization has been applied, produced by NormalizeSomaticReadCounts: tn",
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected String tangentFile;

    @Argument(
            doc = "Genomic targets before tangent normalization file, produced by NormalizeSomaticReadCounts: preTN",
            shortName = ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.PRE_TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected String preTangentFile;

    @Argument(
            doc = "File of segmented regions of the genome",
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            optional = false
    )
    protected String segmentFile;

    @Argument(
            doc = "Prefix for output image files.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Directory to write plots",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected String outputDir;

    @Argument(
            doc = "If input data has had a log2 transform applied",
            shortName = ExomeStandardArgumentDefinitions.LOG2_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.LOG2_LONG_NAME,
            optional = true
    )
    protected Boolean log = true;

    @Argument(
            doc = "Plot sex chromosomes",
            shortName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_LONG_NAME,
            optional = true
    )
    protected Boolean sexChrs = false;

    @Override
    protected Object doWork() {
        final String sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(new File(tangentFile));
        createPlot(sampleName, tangentFile, preTangentFile, segmentFile, outputDir, outputPrefix, log, sexChrs);
        return "Success";
    }

    private void createPlot(String sampleName, String tangentFile, String preTangentFile, String segmentFile, String outDir, String outPrefix, boolean log, boolean sexChrs){
        CopyRatioSegmentedPlotter.writeSegmentedCopyRatioPlot(sampleName, tangentFile, preTangentFile, segmentFile, outDir, outPrefix, log, sexChrs);
    }
}
