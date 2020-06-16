package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GATKTool;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.PrintWriter;

@CommandLineProgramProperties(
        summary = "Prints the header from the input SAM/BAM/CRAM file to a textual output file",
        oneLineSummary = "Print the header from a SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public class PrintReadsHeader extends GATKTool {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "file to write the bam header to", optional = false)
    private GATKPath outputFile;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public void traverse() {
        final SAMFileHeader bamHeader = getHeaderForReads();

        try ( final PrintWriter outputWriter = new PrintWriter(outputFile.getOutputStream()) ) {
            final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            codec.encode(outputWriter, bamHeader);
        }

        logger.info("Successfully wrote reads header to destination: " + outputFile);
    }
}
