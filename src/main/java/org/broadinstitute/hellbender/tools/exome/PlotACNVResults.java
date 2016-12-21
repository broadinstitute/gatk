package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.plotter.ACNVPlotter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Create plots of allele fraction data used for finding copy number variants. Please note that this tool is only supported for hg19 and b37 references. All other references may fail.",
        oneLineSummary = "Create plots of allele fraction data",
        programGroup = CopyNumberProgramGroup.class
)
public final class PlotACNVResults extends CommandLineProgram {

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    @Argument(
            doc = "File of het SNP positions, ref counts, and alt counts, produced by GetHetCoverage/GetBayesianHetCoverage.",
            fullName =  ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "File of tangent-normalized coverage of targets.",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File coverageFile;

    @Argument(
            doc = "File of segmented regions of the genome, produced by AllelicCNV.",
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File segmentsFile;

    @Argument(
            doc = "Prefix for output image files.",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Directory to write plots.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected String outputDir;

    @Argument(
            doc = "Plot sex chromosomes.",
            fullName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
            optional = true
    )
    protected Boolean useSexChromosomes = false;

    @Override
    protected Object doWork() {
        final String sampleName = SegmentUtils.getSampleNameForCLIsFromSegmentFile(segmentsFile);
        Utils.regularReadableUserFile(snpCountsFile);
        Utils.regularReadableUserFile(coverageFile);
        Utils.regularReadableUserFile(segmentsFile);
        ACNVPlotter.writeSegmentedAlleleFractionPlot(sampleName, snpCountsFile, coverageFile,
                segmentsFile, outputDir, outputPrefix, useSexChromosomes);
        ACNVPlotter.writeSegmentedAlleleFractionPlotPerSeg(sampleName, snpCountsFile, coverageFile,
                segmentsFile, outputDir, outputPrefix, useSexChromosomes);
        return "SUCCESS";
    }
}
