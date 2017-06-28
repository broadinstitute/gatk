package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.HashSet;
import java.util.Set;

public class ReadMetadataTest extends BaseTest {
    @Test(groups = "spark")
    void testEverything() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(2, 1, 10000000, 1);
        final String chr1Name = header.getSequenceDictionary().getSequence(0).getSequenceName();
        final String chr2Name = header.getSequenceDictionary().getSequence(1).getSequenceName();
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final ReadMetadata.LibraryFragmentStatistics statistics = new ReadMetadata.LibraryFragmentStatistics(400, 175, 20);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>(3);
        crossContigIgnoreSet.add(1);
        final ReadMetadata readMetadata = new ReadMetadata(crossContigIgnoreSet, header, statistics, null, 1L, 1L, 1);
        Assert.assertEquals(readMetadata.getContigID(chr1Name), 0);
        Assert.assertEquals(readMetadata.getContigID(chr2Name), 1);
        Assert.assertFalse(readMetadata.ignoreCrossContigID(0));
        Assert.assertTrue(readMetadata.ignoreCrossContigID(1));
        Assert.assertThrows(() -> readMetadata.getContigID("not a real name"));
        Assert.assertEquals(readMetadata.getStatistics(groupName), statistics);
        Assert.assertThrows(() -> readMetadata.getStatistics("not a real name"));
    }

    @Test(groups = "spark")
    void serializationTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final ReadMetadata.LibraryFragmentStatistics statistics = new ReadMetadata.LibraryFragmentStatistics(400, 175, 20);
        final Set<Integer> crossContigIgnoreSet = new HashSet<>(3);
        crossContigIgnoreSet.add(0);
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigIgnoreSet, header, statistics, new ReadMetadata.PartitionBounds[0], 1L, 1L, 1);

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, readMetadata);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        final ReadMetadata readMetadata2 = (ReadMetadata)kryo.readClassAndObject(in);
        Assert.assertEquals(readMetadata, readMetadata2);
    }
}
