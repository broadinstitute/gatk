package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.WorkflowProperties;
import org.broadinstitute.barclay.argparser.WorkflowOutput;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.transformers.MisencodedBaseQualityReadTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

@CommandLineProgramProperties(
	summary = "Fixes Illumina base quality scores in a SAM/BAM/CRAM file",
    oneLineSummary = "Fix Illumina base quality scores in a SAM/BAM/CRAM file",
    programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public final class FixMisencodedBaseQualityReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath OUTPUT;

    private SAMFileGATKReadWriter outputWriter;

    private ReadTransformer transform;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT, true);
        transform = new MisencodedBaseQualityReadTransformer();
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(transform.apply(read));
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
