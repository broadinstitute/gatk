package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;

import java.io.File;

public interface QuerynameSetComparison {
    // Process the case where the queryname set in input1 is not found in input2 (move this to the interface)
    public void processInput1(final ReadPair input1ReadPair);

    public void processInput2(final ReadPair input2ReadPair);

    public void processMatchingQuerynameSets(final ReadPair input1ReadPair, final ReadPair input2ReadPair);

    public void writeSummary(final File outputTable, final SAMFileHeader header);
}
