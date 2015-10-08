package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.markduplicates.AbstractMarkDuplicatesTester;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test MarkDuplicates with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as see AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesTest.
 */
public final class MarkDuplicatesDataflowTester extends AbstractMarkDuplicatesTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesDataflow(); }

    @Override
    protected void addArgs() {
        super.addArgs();
        DataflowCommandLineProgramTest.addDataflowRunnerArgs(args);
    }
}
