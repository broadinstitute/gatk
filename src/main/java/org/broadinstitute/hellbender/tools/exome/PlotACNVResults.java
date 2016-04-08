package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.utils.plotter.ACNVPlotter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Create plots of allele fraction data used for finding copy number variants.  Please note that this tool is only supported for hg19 and b37 references.  All other references may fail.",
        oneLineSummary = "Create plots of allele fraction data.",
        programGroup = CopyNumberProgramGroup.class
)
public final class PlotACNVResults extends CommandLineProgram {

    @Argument(
            doc = "File of het SNP positions, ref counts, and alt counts, produced by GetHetCoverage.",
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = false
    )
    protected String snpCountsFile;

    @Argument(
            doc = "File of tangent-normalized coverage of targets",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected String coverageFile;

    @Argument(
            doc = "File of segmented regions of the genome, produced by AllelicCNV.",
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            fullName =  ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            optional = false
    )
    protected String segmentsFile;

    @Argument(
            doc = "Directory to write plots",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected String outputDir;

    @Argument(
            doc = "Plot sex chromosomes",
            shortName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.INCLUDE_SEX_CHROMOSOMES_LONG_NAME,
            optional = true
    )
    protected Boolean useSexChromosomes = false;

    @Override
    protected Object doWork() {
        final String sampleName = SegmentUtils.getSampleNameForCLIsFromSegmentFile(new File(segmentsFile));
        ACNVPlotter.writeSegmentedAlleleFractionPlot(sampleName, snpCountsFile, coverageFile,
                segmentsFile, outputDir, useSexChromosomes);
        return "Success";
    }
}
