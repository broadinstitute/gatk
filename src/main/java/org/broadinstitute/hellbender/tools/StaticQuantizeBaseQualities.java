package org.broadinstitute.hellbender.tools;

import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
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
import org.broadinstitute.hellbender.transformers.BaseQualityStaticQuantizationTransformer;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

/**
 * Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file while performing static quality score quantization
 */
@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file with static quality score quantization",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file with static quality score quantization",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@RuntimeProperties(memory = "1GB")
public final class StaticQuantizeBaseQualities extends ReadWalker {

    @WorkflowResource(input=false, output=true, companionResources={StandardArgumentDefinitions.OUTPUT_LONG_NAME + "Index"})
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    public GATKPath output;
    private SAMFileGATKReadWriter outputWriter;

    /**
     * Static quantized quals are entirely separate from the quantize_qual option which uses dynamic binning.
     * The two types of binning should not be used together.
     */
    @Argument(fullName="static-quantized-quals", doc = "Use static quantized quality scores, default (10,20,30,40)", optional=true)
    public List<Integer> staticQuantizationQuals = new ArrayList<>();

    /**
     * Round down quantized only works with the static_quantized_quals option, and should not be used with
     * the dynamic binning option provided by quantize_quals.  When roundDown = false, rounding is done in
     * probability space to the nearest bin.  When roundDown = true, the value is rounded to the nearest bin
     * that is smaller than the current bin.
     */
    @Argument(fullName="round-down-quantized", doc = "Round quals down to nearest quantized qual", optional=true)
    public boolean roundDown = false;

    @Override
    public ReadTransformer makePostReadFilterTransformer(){
        return new BaseQualityStaticQuantizationTransformer(staticQuantizationQuals, roundDown);
    }
        
    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);
        if (staticQuantizationQuals.size() == 0) {
            staticQuantizationQuals.add(10);
            staticQuantizationQuals.add(20);
            staticQuantizationQuals.add(30);
            staticQuantizationQuals.add(40);
        }
    }

    @Override
    public void apply( final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
