package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceBCICodec;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.FeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

public class DiscordantPairEvidenceUnitTest {
    private static final SAMSequenceDictionary dict = new SAMSequenceDictionary();
    private static final List<String> samples = new ArrayList<>(3);
    private static final List<DiscordantPairEvidence> pairs = new ArrayList<>(18);
    static {
        dict.addSequence(new SAMSequenceRecord("21", 46709983));
        dict.addSequence(new SAMSequenceRecord("22", 50818468));
        dict.addSequence(new SAMSequenceRecord("X", 156040895));
        samples.add("A");
        samples.add("B");
        samples.add("C");
        pairs.add(new DiscordantPairEvidence("C","21",25972981,false,"X",123051950,true));
        pairs.add(new DiscordantPairEvidence("C","21",25974416,false,"X",433179,false));
        pairs.add(new DiscordantPairEvidence("B","21",25974526,true,"X",7510256,false));
        pairs.add(new DiscordantPairEvidence("A","21",25978689,false,"X",23061675,true));
        pairs.add(new DiscordantPairEvidence("C","21",25980279,true,"X",118908694,true));
        pairs.add(new DiscordantPairEvidence("C","21",25991097,false,"X",19552859,true));
        pairs.add(new DiscordantPairEvidence("B","21",25996526,true,"21",25997312,false));
        pairs.add(new DiscordantPairEvidence("B","21",25998677,true,"21",25999518,false));
        pairs.add(new DiscordantPairEvidence("A","21",25999457,true,"21",26000320,false));
        pairs.add(new DiscordantPairEvidence("B","22",30000459,true,"X",112737301,false));
        pairs.add(new DiscordantPairEvidence("A","22",30000461,false,"X",70214634,false));
        pairs.add(new DiscordantPairEvidence("C","22",30000464,true,"X",1113557,false));
        pairs.add(new DiscordantPairEvidence("B","22",30000670,true,"22",30001447,false));
        pairs.add(new DiscordantPairEvidence("C","22",30000827,false,"X",100936820,false));
        pairs.add(new DiscordantPairEvidence("A","22",30003590,false,"X",10293,false));
        pairs.add(new DiscordantPairEvidence("A","22",30006140,true,"22",30007026,false));
        pairs.add(new DiscordantPairEvidence("B","22",30006209,false,"X",116263582,false));
        pairs.add(new DiscordantPairEvidence("C","22",30009296,false,"X",141138844,true));
    }

    @Test
    public void testTextRoundTrip() {
        final DiscordantPairEvidenceCodec codec = new DiscordantPairEvidenceCodec();
        for ( final DiscordantPairEvidence pair : pairs ) {
            Assert.assertEquals(codec.decode(DiscordantPairEvidenceCodec.encode(pair)), pair);
        }
    }

    @Test
    public void testBinaryRoundTrip() {
        final DiscordantPairEvidenceBCICodec codec = new DiscordantPairEvidenceBCICodec();
        final ByteArrayOutputStream os = new ByteArrayOutputStream(1024);
        final FeaturesHeader header = new FeaturesHeader(DiscordantPairEvidence.class.getSimpleName(),
                                            DiscordantPairEvidence.BCI_VERSION, dict, samples);
        final Writer<DiscordantPairEvidence> writer =
                new Writer<>("in-memory stream", os, header, codec::encode);
        for ( final DiscordantPairEvidence be : pairs ) {
            writer.write(be);
        }
        writer.close();
        final ByteArraySeekableStream ss = new ByteArraySeekableStream(os.toByteArray());
        final Reader<DiscordantPairEvidence> reader =
                new Reader<>("in-memory stream", ss, codec);
        final List<DiscordantPairEvidence> recoveredDiscordantPairs = new ArrayList<>(18);
        while ( reader.hasNext() ) {
            recoveredDiscordantPairs.add(reader.readStream());
        }
        Assert.assertEquals(recoveredDiscordantPairs, pairs);
    }
}
