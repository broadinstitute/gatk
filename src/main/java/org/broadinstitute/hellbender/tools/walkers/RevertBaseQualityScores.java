package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.Collections;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Simple tool to revert the quality scores in a SAM/BAM/CRAM file. Copies the scores from the OQ tag to the quality scores.",
        oneLineSummary = "Revert Quality Scores in a SAM/BAM/CRAM file",
        usageExample = "hellbender RevertQualityScores -I input.bam -O output.bam",
        programGroup = ReadDataManipulationProgramGroup.class
)

public class RevertBaseQualityScores extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public String OUTPUT;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(IOUtils.getPath(OUTPUT), true);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        final byte[] originalQuals = ReadUtils.getOriginalBaseQualities(read);
        if ( originalQuals != null ){
            read.setBaseQualities(originalQuals);
        } else {
            throw new UserException("RevertQualityScores can only be applied to SAM/BAM files with original quality scores, "
                    + "caused by read: " + read.getName());
        }
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
