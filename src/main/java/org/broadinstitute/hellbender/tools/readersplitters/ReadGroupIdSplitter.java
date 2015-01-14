package org.broadinstitute.hellbender.tools.readersplitters;

import htsjdk.samtools.SAMReadGroupRecord;

import java.util.function.Function;

/**
 * Splits readers read group id.
 */
public class ReadGroupIdSplitter extends ReadGroupSplitter<String> {
    @Override
    protected Function<SAMReadGroupRecord, String> getSplitByFunction() {
        return SAMReadGroupRecord::getId;
    }
}
