package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

import static org.testng.Assert.*;

/**
 * Created by David Benjamin on 2/15/17.
 */
public class GATKProtectedVariantContextUtilsUnitTest {
    @Test
    public void testGetPileup() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        final Locatable loc = new SimpleInterval("chr1", 10, 10);
        final int readLength = 3;

        //this read doesn't overlap {@code loc}
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1, readLength);
        read1.setBases(Utils.dupBytes((byte) 'A', readLength));
        read1.setBaseQualities(Utils.dupBytes((byte) 30, readLength));
        read1.setCigar("3M");

        //this read overlaps {@code loc} with a Q30 'A'
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read2", 0, 10, readLength);
        read2.setBases(Utils.dupBytes((byte) 'A', readLength));
        read2.setBaseQualities(Utils.dupBytes((byte) 30, readLength));
        read2.setCigar("3M");

        //this read overlaps {@code loc} with a Q20 'C' (due to the deletion)
        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(header, "read3", 0, 7, readLength);
        read3.setBases(Utils.dupBytes((byte) 'C', readLength));
        read3.setBaseQualities(Utils.dupBytes((byte) 20, readLength));
        read3.setCigar("1M1D2M");

        //this read doesn't overlap {@code loc} due to the deletion
        final GATKRead read4 = ArtificialReadUtils.createArtificialRead(header, "read4", 0, 7, readLength);
        read4.setBases(Utils.dupBytes((byte) 'C', readLength));
        read4.setBaseQualities(Utils.dupBytes((byte) 20, readLength));
        read4.setCigar("1M5D2M");

        final ReadPileup pileup = GATKProtectedVariantContextUtils.getPileup(loc, Arrays.asList(read1, read2, read3, read4));

        // the pileup should contain a Q30 'A' and a Q20 'C'
        final int[] counts = pileup.getBaseCounts();
        Assert.assertEquals(counts, new int[]{1, 1, 0, 0});

    }
}