package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.*;

public class PileupReadErrorCorrectorUnitTest {

    @Test
    public void TestErrorCorrection() {
        Utils.resetRandomGenerator();
        final byte[] ref = "AGCCTTGAACTCCTGGGCTCAAGTGATCCTCCTGCCTCAGTTTCCCATGTAGCTGGGACCACAGGTGGGGGCTCCACCCCTGGCTGATTTTTTTTTTTTTTTTTTTTTGAGATAGGGT".getBytes();

        final String sample = "sample";
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rgID");
        rg.setSample(sample);
        rg.setPlatform("illumina");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);

        final int numGoodReads = 500;
        final int numBadReads = 30;
        final int readLength = 30;
        final int refOffset = 100;
        final List<GATKRead> reads = new ArrayList<GATKRead>(numGoodReads + numBadReads);
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals,(byte)30);

        for (int k=0; k < numGoodReads + numBadReads; k++) {
            final int readStart = Utils.getRandomGenerator().nextInt(ref.length - readLength);
            final byte[] bases = Arrays.copyOfRange(ref, readStart,readStart + readLength);

            if (k >= numGoodReads) {
                final int errorPosition = Utils.getRandomGenerator().nextInt(readLength);
                bases[errorPosition] = BaseUtils.simpleComplement(bases[errorPosition]);
            }

            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, bases, quals,readLength+"M");
            read.setName("read" + k);
            read.setReadGroup(rg.getReadGroupId());
            read.setPosition(read.getContig(), refOffset + readStart);
            reads.add(read);
        }

        Collections.sort(reads, Comparator.comparingInt(GATKRead::getStart));

        // now correct all reads
        final ReadErrorCorrector readErrorCorrector = new PileupReadErrorCorrector(3.0, header);
        final List<GATKRead> correctedReads = readErrorCorrector.correctReads(reads);

        for (final GATKRead read : correctedReads) {
            final int startOnRef = read.getStart() - refOffset;
            for (int n = 0; n < readLength; n++) {
                Assert.assertEquals(read.getBase(n), ref[startOnRef + n]);
            }
        }
    }
}