package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.runners.PipelineRunner;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;


public class DataflowTool extends CommandLineProgram {
    private enum PipelineRunnerType {
        LOCAL, BLOCKING, NONBLOCKING
    }

    @Argument(doc="Run ")
    PipelineRunnerType runner = PipelineRunnerType.LOCAL;

    @Override
    protected Object doWork() {

    }
}
