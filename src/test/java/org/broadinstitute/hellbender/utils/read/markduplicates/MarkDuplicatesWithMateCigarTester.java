package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.tools.picard.sam.markduplicates.MarkDuplicatesWithMateCigar;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test MarkDuplicatesWithMateCigar with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as see AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesWithMateCigarTest.
 */
public class MarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicatesTester {

    public MarkDuplicatesWithMateCigarTester() {
        // NB: to be equivalent to MarkDuplicates we need to use SUM_OF_BASE_QUALITIES
        super(ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);

        addArg("--MAX_RECORDS_IN_RAM", "1000");
        addArg("--BLOCK_SIZE", "250");
    }

    @Override
    protected PicardCommandLineProgram getProgram() { return new MarkDuplicatesWithMateCigar(); }
}
