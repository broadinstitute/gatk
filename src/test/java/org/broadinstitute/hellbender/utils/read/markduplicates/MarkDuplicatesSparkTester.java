package org.broadinstitute.hellbender.utils.read.markduplicates;

import com.google.common.collect.Lists;
import htsjdk.samtools.DuplicateScoringStrategy;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesTester;

import java.util.Arrays;
import java.util.List;

/**
 * A tester class for {@link MarkDuplicatesSpark}.
 */
public final class MarkDuplicatesSparkTester extends AbstractMarkDuplicatesTester {
    boolean markUnmappedReads = false;

    public MarkDuplicatesSparkTester() {
        this(false);
    }

    public MarkDuplicatesSparkTester(boolean markUnmappedReads) {
        super(DuplicateScoringStrategy.ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
        this.markUnmappedReads = markUnmappedReads;
    }

    @Override
    protected CommandLineProgram getProgram() { return new MarkDuplicatesSpark(); }

    @Override
    protected void addArgs() {
        if (!markUnmappedReads) {
            addArg("--" + MarkDuplicatesSpark.DO_NOT_MARK_UNMAPPED_MATES);
        }
    }
}
