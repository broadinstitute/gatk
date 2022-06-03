package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceBCICodec;
import org.broadinstitute.hellbender.utils.codecs.BafEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

public class BafEvidenceUnitTest {
    private static final SAMSequenceDictionary dict = new SAMSequenceDictionary();
    private static final List<String> samples = new ArrayList<>(3);
    private static final List<BafEvidence> bafs = new ArrayList<>(18);
    static {
        dict.addSequence(new SAMSequenceRecord("21", 46709983));
        dict.addSequence(new SAMSequenceRecord("22", 50818468));
        samples.add("A");
        samples.add("B");
        samples.add("C");
        bafs.add(new BafEvidence("B","21",25993443,0.5185185185185185));
        bafs.add(new BafEvidence("A","21",25995506,0.43333333333333335));
        bafs.add(new BafEvidence("A","21",25996383,0.46875));
        bafs.add(new BafEvidence("A","21",25996446,0.6));
        bafs.add(new BafEvidence("A","21",25996679,0.5263157894736842));
        bafs.add(new BafEvidence("A","21",25996772,0.3751515151515151));
        bafs.add(new BafEvidence("A","21",25996907,0.5660377358490566));
        bafs.add(new BafEvidence("A","21",25997040,0.4753753753753754));
        bafs.add(new BafEvidence("A","21",25997094,0.518984337921215));
        bafs.add(new BafEvidence("B","22",30000353,0.5666666666666667));
        bafs.add(new BafEvidence("B","22",30003011,0.5555555555555556));
        bafs.add(new BafEvidence("B","22",30004712,0.5526315789473685));
        bafs.add(new BafEvidence("B","22",30005667,0.5925925925925926));
        bafs.add(new BafEvidence("B","22",30007780,0.36363636363636365));
        bafs.add(new BafEvidence("C","22",30010391,0.514162077104642));
        bafs.add(new BafEvidence("B","22",30012721,0.34782608695652173));
        bafs.add(new BafEvidence("C","22",30012825,0.6266666666666667));
        bafs.add(new BafEvidence("B","22",30016476,0.18181818181818182));
    }

    @Test
    public void testTextRoundTrip() {
        final BafEvidenceCodec codec = new BafEvidenceCodec();
        for ( final BafEvidence bafEvidence : bafs ) {
            // text codec prints just two significant digits to economize on file size
            double roundedValue = Math.round(100*bafEvidence.getValue())/100.;
            final BafEvidence be = new BafEvidence(bafEvidence, roundedValue);
            Assert.assertEquals(codec.decode(BafEvidenceCodec.encode(be)), be);
        }
    }

    @Test
    public void testBinaryRoundTrip() {
        final BafEvidenceBCICodec codec = new BafEvidenceBCICodec();
        final ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
        final SVFeaturesHeader header =
                new SVFeaturesHeader(BafEvidence.class.getSimpleName(), BafEvidence.BCI_VERSION, dict, samples);
        final Writer<BafEvidence> writer =
                new Writer<>("in-memory stream", os, header, codec::encode);
        for ( final BafEvidence be : bafs ) {
            writer.write(be);
        }
        writer.close();
        final ByteArraySeekableStream ss = new ByteArraySeekableStream(os.toByteArray());
        final Reader<BafEvidence> reader =
                new Reader<>("in-memory stream", ss, codec);
        final List<BafEvidence> recoveredBafs = new ArrayList<>(18);
        while ( reader.hasNext() ) {
            recoveredBafs.add(reader.readStream());
        }
        Assert.assertEquals(recoveredBafs, bafs);
    }

    @Test
    public void testValueAlteringConstructor() {
        final BafEvidence bafEvidence = new BafEvidence("sample", "contig", 1234, .4);
        final BafEvidence newValue = new BafEvidence(bafEvidence, .5);
        Assert.assertEquals(newValue.getSample(), bafEvidence.getSample());
        Assert.assertEquals(newValue.getContig(), bafEvidence.getContig());
        Assert.assertEquals(newValue.getStart(), bafEvidence.getStart());
        Assert.assertEquals(bafEvidence.getValue(), .4);
        Assert.assertEquals(newValue.getValue(), .5);
    }
}
