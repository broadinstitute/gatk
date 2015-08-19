package org.broadinstitute.hellbender.tools.readersplitters;

import htsjdk.samtools.SAMReadGroupRecord;

import java.util.function.Function;

/**
 * Splits readers by library name.
 */
public final class LibraryNameSplitter extends ReadGroupSplitter<String> {
    @Override
    protected Function<SAMReadGroupRecord, String> getSplitByFunction() {
        return SAMReadGroupRecord::getLibrary;
    }
}
