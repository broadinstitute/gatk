package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.sv.SVFeaturesHeader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Collections;

public class BlockCompressIntervalStreamUnitTest extends GATKBaseTest {
    private static int POS_START = 1000;
    private static int POS_INC = 1100;
    private static int POS_END = 10000000;
    private static int FEATURE_LENGTH = 2000;

    private static final SAMSequenceDictionary dict = new SAMSequenceDictionary();
    static {
        dict.addSequence(new SAMSequenceRecord("21", 46709983));
        dict.addSequence(new SAMSequenceRecord("22", 50818468));
    }

    private static void write( final SimpleFeature feature, final Writer<SimpleFeature> writer ) throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getContigIndex(feature.getContig()));
        dos.writeInt(feature.getStart());
        dos.writeInt(feature.getEnd());
    }

    private static class SimpleFeatureCodec implements FeatureCodec<SimpleFeature, Reader<SimpleFeature>> {

        @Override
        public Feature decodeLoc( Reader<SimpleFeature> simpleFeatureReader ) throws IOException {
            return null;
        }

        @Override
        public SimpleFeature decode( Reader<SimpleFeature> reader ) throws IOException {
            final DataInputStream dis = reader.getStream();
            final String contig = reader.getSequenceNames().get(dis.readInt());
            return new SimpleFeature(contig, dis.readInt(), dis.readInt());
        }

        @Override
        public FeatureCodecHeader readHeader( Reader<SimpleFeature> simpleFeatureReader ) throws IOException {
            return null;
        }

        @Override
        public Class<SimpleFeature> getFeatureType() {
            return null;
        }

        @Override
        public Reader<SimpleFeature> makeSourceFromStream( InputStream bufferedInputStream ) {
            return null;
        }

        @Override
        public LocationAware makeIndexableSourceFromStream( InputStream inputStream ) {
            return null;
        }

        @Override
        public boolean isDone( Reader<SimpleFeature> simpleFeatureReader ) {
            return false;
        }

        @Override
        public void close( Reader<SimpleFeature> simpleFeatureReader ) {

        }

        @Override
        public boolean canDecode( String path ) {
            return false;
        }
    }

    @Test
    public void testRoundTrip() {
        final ByteArrayOutputStream os = new ByteArrayOutputStream(200000);
        final SVFeaturesHeader header =
                new SVFeaturesHeader(SimpleFeature.class.getSimpleName(), "1", dict, Collections.singletonList("sample"));
        final Writer<SimpleFeature> writer =
                new Writer<>("in-memory stream", os, header, BlockCompressIntervalStreamUnitTest::write);
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String contig = rec.getSequenceName();
            for ( int start = POS_START; start < POS_END; start += POS_INC ) {
                final SimpleFeature feature = new SimpleFeature(contig, start, start + FEATURE_LENGTH);
                writer.write(feature);
            }
        }
        writer.close();
        final ByteArraySeekableStream ss = new ByteArraySeekableStream(os.toByteArray());
        final SimpleFeatureCodec codec = new SimpleFeatureCodec();
        final BlockCompressedIntervalStream.Reader<SimpleFeature> reader =
                new BlockCompressedIntervalStream.Reader<>("in-memory stream", ss, codec);
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String contig = rec.getSequenceName();
            for ( int start = POS_START; start < POS_END; start += POS_INC ) {
                final SimpleFeature feature = new SimpleFeature(contig, start, start + FEATURE_LENGTH);
                final SimpleFeature recoveredFeature = reader.readStream();
                Assert.assertEquals(feature.getContig(), recoveredFeature.getContig());
                Assert.assertEquals(feature.getStart(), recoveredFeature.getStart());
                Assert.assertEquals(feature.getEnd(), recoveredFeature.getEnd());
            }
        }
        reader.close();
    }

    @Test
    public void testQuery() throws IOException {
        final ByteArrayOutputStream os = new ByteArrayOutputStream(200000);
        final SVFeaturesHeader header =
                new SVFeaturesHeader(SimpleFeature.class.getSimpleName(), "1", dict, Collections.singletonList("sample"));
        final Writer<SimpleFeature> writer =
                new Writer<>("in-memory stream", os, header, BlockCompressIntervalStreamUnitTest::write);

        // write an initial gigantic feature
        final String contig0 = dict.getSequence(0).getSequenceName();
        final SimpleFeature bigFeature = new SimpleFeature(contig0, POS_START - 10, POS_END);
        writer.write(bigFeature);

        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final String contig = rec.getSequenceName();
            for ( int start = POS_START; start < POS_END; start += POS_INC ) {
                final SimpleFeature feature = new SimpleFeature(contig, start, start + FEATURE_LENGTH);
                writer.write(feature);
            }
        }
        writer.close();

        final ByteArraySeekableStream ss = new ByteArraySeekableStream(os.toByteArray());
        final SimpleFeatureCodec codec = new SimpleFeatureCodec();
        final BlockCompressedIntervalStream.Reader<SimpleFeature> reader =
                new BlockCompressedIntervalStream.Reader<>("in-memory stream", ss, codec);
        final CloseableTribbleIterator<SimpleFeature> itr =
                reader.query(contig0, POS_START + 10, POS_START + 20);
        Assert.assertTrue(itr.hasNext());
        final SimpleFeature recoveredBigFeature = itr.next();
        Assert.assertEquals(bigFeature.getContig(), recoveredBigFeature.getContig());
        Assert.assertEquals(bigFeature.getStart(), recoveredBigFeature.getStart());
        Assert.assertEquals(bigFeature.getEnd(), recoveredBigFeature.getEnd());
        Assert.assertTrue(itr.hasNext());
        final SimpleFeature feature = new SimpleFeature(contig0, POS_START, POS_START + FEATURE_LENGTH);
        final SimpleFeature recoveredFeature = itr.next();
        Assert.assertEquals(feature.getContig(), recoveredFeature.getContig());
        Assert.assertEquals(feature.getStart(), recoveredFeature.getStart());
        Assert.assertEquals(feature.getEnd(), recoveredFeature.getEnd());
        Assert.assertFalse(itr.hasNext());
        itr.close();
        reader.close();

        final String contig1 = dict.getSequence(1).getSequenceName();
        // can't requery an in-memory stream -- start fresh
        final ByteArraySeekableStream ss2 = new ByteArraySeekableStream(os.toByteArray());
        final BlockCompressedIntervalStream.Reader<SimpleFeature> reader2 =
                new BlockCompressedIntervalStream.Reader<>("in-memory stream", ss2, codec);
        final CloseableTribbleIterator<SimpleFeature> itr2 =
                reader2.query(contig1, POS_START + 10, POS_START + 20);
        final SimpleFeature feature1 = new SimpleFeature(contig1, POS_START, POS_START + FEATURE_LENGTH);
        final SimpleFeature recoveredFeature1 = itr2.next();
        Assert.assertEquals(feature1.getContig(), recoveredFeature1.getContig());
        Assert.assertEquals(feature1.getStart(), recoveredFeature1.getStart());
        Assert.assertEquals(feature1.getEnd(), recoveredFeature1.getEnd());
        Assert.assertFalse(itr2.hasNext());
        itr2.close();
        reader2.close();
    }
}
