package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowResource;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.RuntimeProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file while hard clipping any soft clipped bases",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file while hard clipping any soft clipped bases",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@RuntimeProperties(memory = "1GB")
public final class HardClipSoftClippedBases extends ReadWalker {

    @WorkflowResource(input=false, output=true, companionResources={StandardArgumentDefinitions.OUTPUT_LONG_NAME + "Index"})
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public GATKPath output;
    private SAMFileGATKReadWriter outputWriter;
        
    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);
    }

    @Override
    public void apply( final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        outputWriter.addRead(ReadClipper.hardClipSoftClippedBases(read));
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
