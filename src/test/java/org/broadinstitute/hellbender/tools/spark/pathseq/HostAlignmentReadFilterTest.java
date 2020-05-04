package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.common.collect.Sets;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Set;

public final class HostAlignmentReadFilterTest {

    private static final int MIN_IDENTITY = 70;
    private final String IGNORED_CONTIG = "chrIgnore";

    @DataProvider(name = "alignmentData")
    public Object[][] getAlignmentData() {
        return new Object[][]{
                {"100M", 0, Boolean.FALSE},
                {"100M", 30, Boolean.FALSE},
                {"100M", 31, Boolean.TRUE},

                {"70M20S", 0, Boolean.FALSE},
                {"69M21S", 0, Boolean.TRUE},
                {"80M20S", 10, Boolean.FALSE},
                {"79M21S", 10, Boolean.TRUE},
                {"80M21S", 11, Boolean.TRUE},
                {"70M30I", 30, Boolean.FALSE},
                {"69M31I", 31, Boolean.TRUE},
                {"80M10D", 10, Boolean.FALSE},
                {"79M11D", 11, Boolean.TRUE},

                {"20S70M", 0, Boolean.FALSE},
                {"21S69M", 0, Boolean.TRUE},
                {"20S80M", 10, Boolean.FALSE},
                {"21S79M", 10, Boolean.TRUE},
                {"21S80M", 11, Boolean.TRUE},
                {"30I70M", 30, Boolean.FALSE},
                {"31I69M", 31, Boolean.TRUE},
                {"10D80M", 10, Boolean.FALSE},
                {"11D79M", 11, Boolean.TRUE},

                {"35M30I35M", 30, Boolean.FALSE},
                {"34M31I35M", 31, Boolean.TRUE},
                {"40M10D40M", 10, Boolean.FALSE},
                {"39M11D40M", 11, Boolean.TRUE}
        };
    }

    @Test(dataProvider = "alignmentData")
    public void testMappedRead(final String cigarString, final int NM, final boolean test_out) {
        final HostAlignmentReadFilter filter = new HostAlignmentReadFilter(MIN_IDENTITY);
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead read_in = ArtificialReadUtils.createArtificialRead(cigar);
        read_in.setAttribute("NM", NM);
        read_in.setPosition("pos", 1);
        final boolean test_i = filter.test(read_in);
        Assert.assertEquals(test_out, test_i);
    }

    @Test
    public void testUnmappedRead() {
        final HostAlignmentReadFilter filter = new HostAlignmentReadFilter(MIN_IDENTITY);
        final byte[] bases = new byte[100];
        final byte[] qual = new byte[100];
        Arrays.fill(bases, (byte) 'A');
        Arrays.fill(qual, (byte) 30);
        final GATKRead read_in = ArtificialReadUtils.createArtificialRead(bases, qual, "*");
        read_in.setIsUnmapped();
        final boolean test_i = filter.test(read_in);
        Assert.assertTrue(test_i);
    }

    @Test
    public void testIgnoredContig() {
        final Set<String> ignoredContigSet = Sets.newHashSet(IGNORED_CONTIG);
        final HostAlignmentReadFilter filter = new HostAlignmentReadFilter(MIN_IDENTITY, ignoredContigSet);
        final byte[] bases = new byte[100];
        final byte[] qual = new byte[100];
        Arrays.fill(bases, (byte) 'A');
        Arrays.fill(qual, (byte) 30);

        // Mapped fully to ignored contig
        final GATKRead readA = ArtificialReadUtils.createArtificialRead(bases, qual, "100M");
        readA.setAttribute("NM", 0);
        readA.setPosition(IGNORED_CONTIG, 1);
        final boolean testA = filter.test(readA);
        Assert.assertTrue(testA);

        // Mapped fully to non-ignored contig
        final GATKRead readB = ArtificialReadUtils.createArtificialRead(bases, qual, "100M");
        readB.setAttribute("NM", 0);
        readB.setPosition("chrNotIgnore", 1);
        final boolean testB = filter.test(readB);
        Assert.assertFalse(testB);

        // Mapped poorly to ignored contig
        final GATKRead readC = ArtificialReadUtils.createArtificialRead(bases, qual, "10M90S");
        readC.setAttribute("NM", 90);
        readC.setPosition(IGNORED_CONTIG, 1);
        final boolean testC = filter.test(readC);
        Assert.assertTrue(testC);

        // Mapped poorly to non-ignored contig
        final GATKRead readD = ArtificialReadUtils.createArtificialRead(bases, qual, "10M90S");
        readD.setAttribute("NM", 90);
        readD.setPosition("chrNotIgnore", 1);
        final boolean testD = filter.test(readD);
        Assert.assertTrue(testD);
    }

    @Test
    public void testMappedReadWithoutNMTag() {
        final String cigarString = "100M";
        final HostAlignmentReadFilter filter = new HostAlignmentReadFilter(MIN_IDENTITY);
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead read_in = ArtificialReadUtils.createArtificialRead(cigar);
        read_in.setPosition("pos", 1);
        final boolean test_i = filter.test(read_in);
        Assert.assertTrue(test_i); //Don't filter out reads if there is no NM tag
    }
}
