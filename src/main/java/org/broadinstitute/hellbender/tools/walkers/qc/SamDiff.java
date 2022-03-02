package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.List;

public class SamDiff implements QuerynameSetComparison {
    final SAMFileGATKReadWriter writer1;
    final SAMFileGATKReadWriter writer2;

    public SamDiff(final SAMFileGATKReadWriter writer1, final SAMFileGATKReadWriter writer2){
        this.writer1 = writer1;
        this.writer2 = writer2;
    }

    @Override
    public void processInput1(ReadPair input1ReadPair) {
        final List<GATKRead> reads = input1ReadPair.getReads(false);
        reads.forEach(r -> writer1.addRead(r));
    }

    @Override
    public void processInput2(ReadPair input2ReadPair) {
        final List<GATKRead> reads = input2ReadPair.getReads(false);
        reads.forEach(r -> writer2.addRead(r));
    }

    @Override
    public void processMatchingQuerynameSets(ReadPair input1ReadPair, ReadPair input2ReadPair) {
        // Do nothing. Perhaps count?
    }

    @Override
    public void writeSummary(File outputTable, SAMFileHeader header) {
        // No need to output this.
    }
}
