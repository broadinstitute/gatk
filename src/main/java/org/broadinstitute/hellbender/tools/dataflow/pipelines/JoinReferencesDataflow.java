package org.broadinstitute.hellbender.tools.dataflow.pipelines;


import com.google.cloud.dataflow.sdk.Pipeline;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.dataflow.transforms.JoinReferencesDataflowTransform;


@CommandLineProgramProperties(summary = "Join reads with references", oneLineSummary = "join references", programGroup = DataFlowProgramGroup.class)
public final class JoinReferencesDataflow extends DataflowReadsPipeline {

    static {
        // use Spark serialization for broadcasts to avoid making a copy of the reference (which is ~3GB)
        System.setProperty("dataflow.spark.directBroadcast", "true");
    }

    private static final long serialVersionUID = 1L;
    private Pipeline pipeline;

    @Argument(doc = "local path for the reference file",
            shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME,
            optional = false)
    protected String reference;

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        this.pipeline = pipeline;
        super.setupPipeline(pipeline);
    }

    @Override
    protected PTransformSAM<Long> getTool() {
        ReferenceFileSource refSource = new ReferenceFileSource(reference);
        return new JoinReferencesDataflowTransform(pipeline, refSource);
    }

}



