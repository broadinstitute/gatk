package org.broadinstitute.hellbender.utils.io;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;
import org.broadinstitute.hellbender.tools.sv.SVIOUtils;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;

public class FeatureOutputStreamUnitTest extends GATKBaseTest {

    private final SAMSequenceDictionary dictionary = ReferenceUtils.loadFastaDictionary(new File(FULL_HG38_DICT));

    @DataProvider
    public Object[][] featureOutputStreamData() {
        return new Object[][] {
                {
                    Lists.newArrayList(
                                new SplitReadEvidence("sample1", "chr1", 4783443, 3, true),
                                new SplitReadEvidence("sample2", "chr1", 4783443, 2, true),
                                new SplitReadEvidence("sample3", "chr1", 4783443, 1, true),
                                new SplitReadEvidence("sample1", "chr1", 6373883, 1, false),
                                new SplitReadEvidence("sample1", "chr1", 8398393, 7, true)
                        ),
                        SplitReadEvidenceCodec.FORMAT_SUFFIX,
                        new SplitReadEvidenceCodec(),
                        null
                },
                {
                        Lists.newArrayList(
                                new DepthEvidence("chr1", 4000, 4099, new int[]{0, 5, 0}),
                                new DepthEvidence("chr1", 4100, 4199, new int[]{0, 3, 0}),
                                new DepthEvidence("chr1", 4200, 4300, new int[]{0, 1, 1})
                        ),
                        DepthEvidenceCodec.FORMAT_SUFFIX,
                        new DepthEvidenceCodec(),
                        "#Chr\tStart\tEnd\tsample1\tsample2\tsample3"
                },
        };
    }

    @Test(dataProvider = "featureOutputStreamData")
    public void testTabixIndexedStream(final ArrayList<Feature> featureList, final String extension,
                                       final FeatureCodec codec, final String expectedHeader) throws IOException {
        final File tempDir = IOUtils.createTempDir(TabixIndexedFeatureOutputStream.class.getSimpleName());
        final Path outFilePath = Paths.get(tempDir.toString(), getClass().getSimpleName() + extension + ".gz");
        final FeatureOutputStream stream = new TabixIndexedFeatureOutputStream(
                new GATKPath(outFilePath.toString()),
                codec,
                SVIOUtils::encodeSVEvidenceFeature,
                dictionary,
                4
        );
        testWithStream(stream, featureList, outFilePath, codec, expectedHeader, true);
    }

    @Test(dataProvider = "featureOutputStreamData")
    public void testUncompressedStream(final ArrayList<Feature> featureList, final String extension,
                                       final FeatureCodec codec, final String expectedHeader) throws IOException {
        final File tempDir = IOUtils.createTempDir(TabixIndexedFeatureOutputStream.class.getSimpleName());
        final Path outFilePath = Paths.get(tempDir.toString(), getClass().getSimpleName() + extension);
        final FeatureOutputStream stream = new UncompressedFeatureOutputStream(
                new GATKPath(outFilePath.toString()),
                SVIOUtils::encodeSVEvidenceFeature
        );
        testWithStream(stream, featureList, outFilePath, codec, expectedHeader, false);
    }

    private void testWithStream(final FeatureOutputStream stream, final ArrayList<Feature> featureList,
                                final Path outFilePath, final FeatureCodec codec,
                                final String expectedHeader, final boolean indexed) throws IOException {
        final Path outIndexPath = Paths.get(outFilePath + FileExtensions.TABIX_INDEX);
        if (expectedHeader != null) {
            stream.writeHeader(expectedHeader);
        }
        featureList.stream().forEachOrdered(s -> stream.add(s));
        stream.close();
        if (indexed) {
            Assert.assertTrue(Files.exists(outIndexPath), "Index does not exist");
        }
        final InputStream inputStream = IOUtil.isBlockCompressed(outFilePath) ?
                new BlockCompressedInputStream(outFilePath.toFile()) : new FileInputStream(outFilePath.toFile());
        final LineIterator reader = new AsciiLineReaderIterator(AsciiLineReader.from(inputStream));

        //Check header
        if (expectedHeader != null) {
            final Object actualHeader = codec.readHeader(reader).getHeaderValue();
            if (actualHeader instanceof String) {
                Assert.assertEquals((String) actualHeader, expectedHeader);
            }
        } else {
            Assert.assertNull(codec.readHeader(reader).getHeaderValue());
        }

        //Check records
        final Iterator<Feature> expectedIterator = featureList.iterator();
        while (reader.hasNext()) {
            Assert.assertTrue(expectedIterator.hasNext(), "Actual file had more lines than expected");
            final Feature feature = codec.decode(reader);
            Assert.assertEquals(feature, expectedIterator.next());
        }
    }
}