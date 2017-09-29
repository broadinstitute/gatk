package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.picard.sam.markduplicates.MarkDuplicates;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesTester;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test MarkDuplicates with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as see AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesTest.
 */
public final class MarkDuplicatesTester extends AbstractMarkDuplicatesTester {

    public MarkDuplicatesTester() {
        super(DuplicateScoringStrategy.ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicates(); }
}
