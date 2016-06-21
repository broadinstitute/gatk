package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.spark.GATKRegistrator;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Collections;

public class ReadMetadataTest extends BaseTest {
    @Test(groups = "spark")
    void testEverything() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final String chr1Name = header.getSequenceDictionary().getSequence(0).getSequenceName();
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final ReadMetadata.ReadGroupFragmentStatistics statistics = new ReadMetadata.ReadGroupFragmentStatistics(400.f, 75.f);
        final ReadMetadata readMetadata = new ReadMetadata(header, Collections.singletonList(statistics), statistics);
        Assert.assertEquals(readMetadata.getContigID(chr1Name), 0);
        Assert.assertThrows(() -> readMetadata.getContigID("not a real name"));
        Assert.assertEquals(readMetadata.getStatistics(groupName), statistics);
        Assert.assertThrows(() -> readMetadata.getStatistics("not a real name"));
    }

    @Test(groups = "spark")
    void serializationTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final ReadMetadata.ReadGroupFragmentStatistics statistics = new ReadMetadata.ReadGroupFragmentStatistics(400.f, 75.f);
        final ReadMetadata readMetadata = new ReadMetadata(header, Collections.singletonList(statistics), statistics);

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
