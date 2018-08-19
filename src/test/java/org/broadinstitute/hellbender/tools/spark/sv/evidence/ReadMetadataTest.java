package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.util.HashSet;
import java.util.Set;

public class ReadMetadataTest extends GATKBaseTest {
    private final static LibraryStatistics LIBRARY_STATISTICS =
            new LibraryStatistics(IntHistogramTest.genLogNormalSample(400, 175, 10000).getCDF(),
                    60000000000L, 600000000L, 1200000000000L, 3000000000L);

    @Test(groups = "sv")
    void testEverything() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(2, 1, 10000000, 1);
        final String chr1Name = header.getSequenceDictionary().getSequence(0).getSequenceName();
        final String chr2Name = header.getSequenceDictionary().getSequence(1).getSequenceName();
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final Set<Integer> crossContigIgnoreSet = new HashSet<>(3);
        crossContigIgnoreSet.add(1);
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigIgnoreSet, header, LIBRARY_STATISTICS, new ReadMetadata.PartitionBounds[0],
                                1L, 1L, 1);
        Assert.assertEquals(readMetadata.getContigID(chr1Name), 0);
        Assert.assertEquals(readMetadata.getContigID(chr2Name), 1);
        Assert.assertFalse(readMetadata.ignoreCrossContigID(0));
        Assert.assertTrue(readMetadata.ignoreCrossContigID(1));
        Assert.assertThrows(() -> readMetadata.getContigID("not a real name"));
        Assert.assertEquals(readMetadata.getLibraryStatistics(readMetadata.getLibraryName(groupName)),
                LIBRARY_STATISTICS);
        Assert.assertThrows(() -> readMetadata.getLibraryName("not a real name"));

        final File metadataFile = BaseTest.createTempFile("metadata", "");
        ReadMetadata.writeMetadata(readMetadata, metadataFile.toString());
    }

    @Test(groups = "sv")
    void serializationTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>(3);
        crossContigIgnoreSet.add(0);
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigIgnoreSet, header, LIBRARY_STATISTICS,
                        new ReadMetadata.PartitionBounds[0], 1L, 1L, 1);

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, readMetadata);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        final ReadMetadata readMetadata2 = (ReadMetadata)kryo.readClassAndObject(in);
        Assert.assertEquals(readMetadata.getCrossContigIgnoreSet(), readMetadata2.getCrossContigIgnoreSet());
        Assert.assertEquals(readMetadata.getContigNameMap(), readMetadata2.getContigNameMap());
        Assert.assertEquals(readMetadata.getReadGroupToLibraryMap(), readMetadata2.getReadGroupToLibraryMap());
        Assert.assertEquals(readMetadata.getNReads(), readMetadata2.getNReads());
        Assert.assertEquals(readMetadata.getMaxReadsInPartition(), readMetadata2.getMaxReadsInPartition());
        Assert.assertEquals(readMetadata.getCoverage(), readMetadata2.getCoverage());
        Assert.assertEquals(readMetadata.getAllPartitionBounds(), readMetadata2.getAllPartitionBounds());
        Assert.assertEquals(readMetadata.getAllLibraryStatistics().keySet(),
                            readMetadata2.getAllLibraryStatistics().keySet());
        for ( final String readGroupName : readMetadata.getReadGroupToLibraryMap().keySet() ) {
            final LibraryStatistics stats1 = readMetadata.getFragmentLengthStatistics(readGroupName);
            final LibraryStatistics stats2 = readMetadata2.getFragmentLengthStatistics(readGroupName);
            Assert.assertEquals(stats1.getMedian(),stats2.getMedian());
            Assert.assertEquals(stats1.getNegativeMAD(),stats2.getNegativeMAD());
            Assert.assertEquals(stats1.getPositiveMAD(),stats2.getPositiveMAD());
        }
    }
}
