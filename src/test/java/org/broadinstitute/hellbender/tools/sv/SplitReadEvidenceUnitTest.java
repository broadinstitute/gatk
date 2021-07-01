package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceBCICodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.FeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

public class SplitReadEvidenceUnitTest {
    private static final SAMSequenceDictionary dict = new SAMSequenceDictionary();
    private static final List<String> samples = new ArrayList<>(3);
    private static final List<SplitReadEvidence> splits = new ArrayList<>(18);
    static {
        dict.addSequence(new SAMSequenceRecord("21", 46709983));
        dict.addSequence(new SAMSequenceRecord("22", 50818468));
        samples.add("A");
        samples.add("B");
        samples.add("C");
        splits.add(new SplitReadEvidence("B","21",25998421,1,true));
        splits.add(new SplitReadEvidence("B","21",25998456,1,true));
        splits.add(new SplitReadEvidence("B","21",25998529,1,true));
        splits.add(new SplitReadEvidence("B","21",25998668,1,false));
        splits.add(new SplitReadEvidence("B","21",25998825,1,true));
        splits.add(new SplitReadEvidence("C","21",25999220,1,false));
        splits.add(new SplitReadEvidence("B","21",25999674,1,false));
        splits.add(new SplitReadEvidence("B","21",25999957,1,false));
        splits.add(new SplitReadEvidence("C","22",30000009,1,false));
        splits.add(new SplitReadEvidence("B","22",30000115,1,false));
        splits.add(new SplitReadEvidence("B","22",30000125,1,false));
        splits.add(new SplitReadEvidence("A","22",30000133,1,false));
        splits.add(new SplitReadEvidence("B","22",30000137,1,false));
        splits.add(new SplitReadEvidence("B","22",30000189,1,false));
        splits.add(new SplitReadEvidence("B","22",30000192,1,false));
        splits.add(new SplitReadEvidence("A","22",30000196,1,false));
        splits.add(new SplitReadEvidence("B","22",30000335,1,false));
        splits.add(new SplitReadEvidence("A","22",30000338,1,false));
    }

    @Test
    public void testTextRoundTrip() {
        final SplitReadEvidenceCodec codec = new SplitReadEvidenceCodec();
        for ( final SplitReadEvidence sr : splits ) {
            Assert.assertEquals(codec.decode(SplitReadEvidenceCodec.encode(sr)), sr);
        }
    }

    @Test
    public void testBinaryRoundTrip() {
        final SplitReadEvidenceBCICodec codec = new SplitReadEvidenceBCICodec();
        final ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
        final FeaturesHeader header =
                new FeaturesHeader(SplitReadEvidence.class.getSimpleName(), SplitReadEvidence.BCI_VERSION, dict, samples);
        final Writer<SplitReadEvidence> writer =
                new Writer<>("in-memory stream", os, header, codec::encode);
        for ( final SplitReadEvidence sr : splits ) {
            writer.write(sr);
        }
        writer.close();
        final ByteArraySeekableStream ss = new ByteArraySeekableStream(os.toByteArray());
        final Reader<SplitReadEvidence> reader =
                new Reader<>("in-memory stream", ss, codec);
        final List<SplitReadEvidence> recoveredSplitReads = new ArrayList<>(18);
        while ( reader.hasNext() ) {
            recoveredSplitReads.add(reader.readStream());
        }
        Assert.assertEquals(recoveredSplitReads, splits);
    }
}
