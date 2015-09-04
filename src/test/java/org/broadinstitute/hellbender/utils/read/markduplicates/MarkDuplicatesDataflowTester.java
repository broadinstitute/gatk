package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.MarkDuplicatesDataflow;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test MarkDuplicates with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as see AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesTest.
 */
public final class MarkDuplicatesDataflowTester extends AbstractMarkDuplicatesTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesDataflow(); }
}
