package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Simple tool to revert the quality scores in a SAM/BAM/CRAM file. Copies the scores from the OQ tag to the quality scores.",
        oneLineSummary = "Revert Quality Scores in a SAM/BAM/CRAM file",
        usageExample = "hellbender RevertQualityScores -I input.bam -O output.bam",
        programGroup = ReadProgramGroup.class
)

public class RevertBaseQualityScores extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT, true);
    }

    @Override
    public CountingReadFilter makeReadFilter() {
        return new CountingReadFilter("Allow all reads", ReadFilterLibrary.ALLOW_ALL_READS);
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
