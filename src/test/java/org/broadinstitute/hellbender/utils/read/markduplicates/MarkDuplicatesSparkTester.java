package org.broadinstitute.hellbender.utils.read.markduplicates;

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;

/**
 * A tester class for {@link MarkDuplicatesSpark}.
 */
public final class MarkDuplicatesSparkTester extends AbstractMarkDuplicatesTester {

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesSpark(); }
}
