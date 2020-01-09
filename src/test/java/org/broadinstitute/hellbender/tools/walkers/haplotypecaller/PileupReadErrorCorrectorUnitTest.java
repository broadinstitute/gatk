package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class PileupReadErrorCorrectorUnitTest {

    @Test
    public void TestErrorCorrection() {
        final String refChunk = "GCATAAACATGGCTCACTGC";
        final String refChunkHard = "AGCCTTGAACTCCTGGGCTCAAGTGATCCTCCTGCCTCAGTTTCCCATGTAGCTGGGACCACAGGTGGGGGCTCCACCCCTGGCTGATTTTTTTTTTTTTTTTTTTTTGAGATAGGGT";

        final String sample = "sample";
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rgID");
        rg.setSample(sample);
        rg.setPlatform("illumina");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);

        final int numGoodReads = 500;
        final int numBadReads = 10;
        final int readLength = 15;
        final List<GATKRead> finalizedReadList = new ArrayList<GATKRead>(numGoodReads);
        int offset = 0;
        final byte[] quals = new byte[readLength];

        Arrays.fill(quals,(byte)30);

        for (int k=0; k < numGoodReads; k++) {
            final byte[] bases = Arrays.copyOfRange(refChunk.getBytes(),offset,offset+readLength);
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, bases, quals,readLength+"M");
            read.setName("good_read" + k);
            read.setReadGroup(rg.getReadGroupId());
            finalizedReadList.add(read);
            offset++;
            if (offset >= refChunk.length()-readLength)
                offset = 0;
        }
        offset = 2;
        // coverage profile is now perfectly triangular with "good" bases. Inject now bad bases with errors in them.
        for (int k=0; k < numBadReads; k++) {
            final byte[] bases = finalizedReadList.get(k).getBases().clone();
            bases[offset] = 'N';
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, bases, quals, readLength + "M");
            read.setName("bad_read" + k);
            read.setReadGroup(rg.getReadGroupId());
            finalizedReadList.add(read);
            offset += 7;
            if (offset >= readLength)
                offset = 4; // just some randomly circulating offset for error position
        }

        // now correct all reads
        final ReadErrorCorrector readErrorCorrector = new PileupReadErrorCorrector(3.0, header);
        readErrorCorrector.correctReads(finalizedReadList);

        // check that corrected reads have exactly same content as original reads
        for (int k=0; k < numBadReads; k++) {
            final byte[] badBases = finalizedReadList.get(k).getBases();
            final byte[] originalBases = finalizedReadList.get(k).getBases();
            Assert.assertTrue(Arrays.equals(badBases,originalBases));
        }
    }

}