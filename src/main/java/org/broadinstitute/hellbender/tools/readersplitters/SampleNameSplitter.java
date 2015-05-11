package org.broadinstitute.hellbender.tools.readersplitters;

import htsjdk.samtools.SAMReadGroupRecord;

import java.util.function.Function;

/**
 * Splits readers sample names.
 */
public final class SampleNameSplitter extends ReadGroupSplitter<String> {
    @Override
    protected Function<SAMReadGroupRecord, String> getSplitByFunction() {
        return SAMReadGroupRecord::getSample;
    }
}
