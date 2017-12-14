package org.broadinstitute.hellbender.tools;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) Emit a single sample name from the bam header into an output file.",
        summary = "If the bam has zero or more than one sample names in the header, this tool will error, by design.\n" +
                "  This tool has not been tested extensively.  Most options supported by the GATK are irrelevant for this tool.",
        programGroup = QCProgramGroup.class
)
@BetaFeature
final public class GetSampleName extends GATKTool{

    @Argument(
            doc = "Output file with only the sample name in it.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outputSampleNameFile;

    @Override
    public void traverse() {
        // Do nothing!
    }

    @Override
    public boolean requiresReads() {return true;}

    @Override
    public void onTraversalStart() {

        // Grab the header info
        if ((getHeaderForReads() == null) || (getHeaderForReads().getReadGroups() == null)) {
            throw new UserException.BadInput("The given input bam has no header or no read groups.  Cannot determine a sample name.");
        }

        final List<String> sampleNames = getHeaderForReads().getReadGroups().stream().map(s -> s.getSample()).distinct().collect(Collectors.toList());
        if (sampleNames.size() > 1) {
            throw new UserException.BadInput("The given input bam has more than one unique sample name: " + StringUtils.join(sampleNames, ", "));
        }

        if (sampleNames.size() == 0) {
            throw new UserException.BadInput("The given bam input has no sample names.");
        }

        try (final FileWriter fileWriter = new FileWriter(outputSampleNameFile, false)) {
            fileWriter.write(sampleNames.get(0));
        } catch (final IOException ioe) {
            throw new UserException("Could not write file.", ioe);
        }
    }
}
