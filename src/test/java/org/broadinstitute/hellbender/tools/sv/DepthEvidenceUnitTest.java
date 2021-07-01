package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceBCICodec;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.FeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

public class DepthEvidenceUnitTest {
    private static final SAMSequenceDictionary dict = new SAMSequenceDictionary();
    private static final List<String> samples = new ArrayList<>(3);
    private static final List<DepthEvidence> depths = new ArrayList<>(18);
    static {
        dict.addSequence(new SAMSequenceRecord("21", 46709983));
        dict.addSequence(new SAMSequenceRecord("22", 50818468));
        samples.add("A");
        samples.add("B");
        samples.add("C");
        depths.add(new DepthEvidence("21",25999208,25999308,new int[]{15,23,18}));
        depths.add(new DepthEvidence("21",25999308,25999408,new int[]{22,22,32}));
        depths.add(new DepthEvidence("21",25999408,25999508,new int[]{28,18,19}));
        depths.add(new DepthEvidence("21",25999508,25999608,new int[]{17,20,29}));
        depths.add(new DepthEvidence("21",25999608,25999708,new int[]{22,17,24}));
        depths.add(new DepthEvidence("21",25999708,25999808,new int[]{26,27,26}));
        depths.add(new DepthEvidence("21",25999808,25999908,new int[]{22,20,22}));
        depths.add(new DepthEvidence("21",25999908,26000008,new int[]{24,20,29}));
        depths.add(new DepthEvidence("22",29999964,30000064,new int[]{29,20,33}));
        depths.add(new DepthEvidence("22",30000064,30000164,new int[]{25,18,25}));
        depths.add(new DepthEvidence("22",30000164,30000264,new int[]{21,21,23}));
        depths.add(new DepthEvidence("22",30000264,30000364,new int[]{35,22,32}));
        depths.add(new DepthEvidence("22",30000364,30000464,new int[]{58,52,51}));
        depths.add(new DepthEvidence("22",30000464,30000564,new int[]{35,35,32}));
        depths.add(new DepthEvidence("22",30000564,30000664,new int[]{20,15,16}));
        depths.add(new DepthEvidence("22",30000664,30000764,new int[]{36,22,26}));
        depths.add(new DepthEvidence("22",30000764,30000864,new int[]{23,20,25}));
        depths.add(new DepthEvidence("22",30000864,30000964,new int[]{16,18,30}));
    }

    @Test
    public void testTextRoundTrip() {
        final DepthEvidenceCodec codec = new DepthEvidenceCodec();
        for ( final DepthEvidence de : depths ) {
            Assert.assertEquals(codec.decode(DepthEvidenceCodec.encode(de)), de);
        }
    }

    @Test
    public void testBinaryRoundTrip() {
        final DepthEvidenceBCICodec codec = new DepthEvidenceBCICodec();
        final ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
        final FeaturesHeader header =
                new FeaturesHeader(DepthEvidence.class.getSimpleName(), DepthEvidence.BCI_VERSION, dict, samples);
        final Writer<DepthEvidence> writer =
                new Writer<>("in-memory stream", os, header, codec::encode);
        for ( final DepthEvidence de : depths ) {
            writer.write(de);
        }
        writer.close();
        final ByteArraySeekableStream ss = new ByteArraySeekableStream(os.toByteArray());
        final Reader<DepthEvidence> reader =
                new Reader<>("in-memory stream", ss, codec);
        final List<DepthEvidence> recoveredDepths = new ArrayList<>(18);
        while ( reader.hasNext() ) {
            recoveredDepths.add(reader.readStream());
        }
        Assert.assertEquals(recoveredDepths, depths);
    }
}
