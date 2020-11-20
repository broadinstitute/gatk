package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

public final class ShardingVCFWriterUnitTest extends GATKBaseTest {

    private static final String testVcfBaseFilename = "test";
    private static final String testVariantId = "test_var";
    private static final String contigName = "00";
    private final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord(contigName, 1000000)));
    private final Path vcfBasePath = Paths.get(createTempDir(getTestedClassName()).getAbsolutePath(), testVcfBaseFilename);
    private final String firstShardPath = ShardingVCFWriter.getShardFilename(vcfBasePath, 0);

    private ShardingVCFWriter createTestWriter(final Path path) {
        return new ShardingVCFWriter(path, 100, dictionary, false, Options.INDEX_ON_THE_FLY);
    }

    private VariantContext createTestVariant(final String id) {
        return new VariantContextBuilder(null, contigName, 1, 1, Collections.singletonList(Allele.REF_N)).id(id).make();
    }

    private void asssertHeaderDictionariesEqual(final VCFHeader header, final String path) {
        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(path);
        final VCFHeader headerTest = vcf.getKey();
        Assert.assertEquals(header.getSequenceDictionary(), headerTest.getSequenceDictionary(),
                "Header sequence dictionaries do not match");
    }

    private void assertVariantCount(final int expectedCount, final String path, final boolean hasHeader) {
        final int count;
        if (hasHeader) {
            final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(path);
            count = vcf.getValue().size();
        } else {
            try (final BufferedReader reader = new BufferedReader(IOUtils.makeReaderMaybeGzipped(Paths.get(path)))) {
                final List<String> lines = reader.lines().collect(Collectors.toList());
                lines.stream().forEach(line -> Assert.assertFalse(line.startsWith("#"), "Unexpected header line"));
                count = lines.size();
            } catch (final IOException e) {
                throw new TestException(e);
            }
        }
        Assert.assertEquals(count, expectedCount, "Incorrect variant count");
    }

    @Test
    public void testWriteHeader() {
        final ShardingVCFWriter writer = createTestWriter(vcfBasePath);
        final VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(dictionary);
        writer.writeHeader(header);
        writer.close();
        asssertHeaderDictionariesEqual(header, firstShardPath);
    }

    @Test(expectedExceptions = RuntimeIOException.class)
    public void testDoubleClose() {
        final ShardingVCFWriter writer = createTestWriter(vcfBasePath);
        writer.close();
        writer.close();
    }

    @Test
    public void testCheckNoError() {
        final VariantContextWriter baseWriter = new VariantContextWriterBuilder()
                .setOutputPath(Paths.get(firstShardPath))
                .setReferenceDictionary(dictionary)
                .setCreateMD5(false)
                .setOption(Options.INDEX_ON_THE_FLY)
                .build();
        final ShardingVCFWriter writer = createTestWriter(vcfBasePath);
        Assert.assertFalse(writer.checkError());
        Assert.assertEquals(baseWriter.checkError(), writer.checkError());
    }

    @Test
    public void testAdd() {
        final ShardingVCFWriter writer = createTestWriter(vcfBasePath);
        final VariantContext variant = createTestVariant(testVariantId);
        final VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(dictionary);
        writer.writeHeader(header);
        writer.add(variant);
        writer.close();

        final Pair<VCFHeader, List<VariantContext>> vcf = VariantContextTestUtils.readEntireVCFIntoMemory(firstShardPath);
        final Iterator<VariantContext> iter = vcf.getValue().iterator();
        Assert.assertTrue(iter.hasNext(), "No variants could be read");
        final VariantContext variantTest = iter.next();
        Assert.assertEquals(variant.getID(), variantTest.getID(), "Test variant does not match");
    }

    @DataProvider(name = "ShardingDataProvider")
    public Object[][] getShardingData() {
        // Files, expected to support serial iteration (false if any input is a .sam)
        return new Object[][] {
                { 10, 25, 3, 5, true },
                { 10, 25, 3, 5, false },
                { 10, 30, 3, 10, true },
                { 10, 11, 2, 1, true },
                { 10, 0, 1, 0, true }
        };
    }

    @Test(dataProvider = "ShardingDataProvider")
    public void testSharding(final int shardSize, final int numVariants, final int numExpectedShards,
                             int numExpectedVariantsLastShard, final boolean writeHeader) {
        //Initialize our writer
        final Path shardBasePath = Paths.get(createTempDir(getTestedClassName() + "_ShardTest").getAbsolutePath(), testVcfBaseFilename);
        final ShardingVCFWriter writer = new ShardingVCFWriter(shardBasePath, shardSize, dictionary, false, Options.INDEX_ON_THE_FLY);
        final VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(dictionary);
        if (writeHeader) {
            writer.writeHeader(header);
        } else {
            writer.setHeader(header);
        }

        //Write test variants
        for (int i = 0; i < numVariants; i++) {
            final VariantContext variant = createTestVariant(testVariantId + "_" + i);
            writer.add(variant);
        }
        writer.close();

        //Check that correct number of shards are written and each is the correct size
        for (int i = 0; i < numExpectedShards; i++) {
            final String shardPath = ShardingVCFWriter.getShardFilename(shardBasePath, i);
            final int expectedCount = i < numExpectedShards - 1 ? shardSize : numExpectedVariantsLastShard;
            if (writeHeader) {
                asssertHeaderDictionariesEqual(header, shardPath);
            }
            assertVariantCount(expectedCount, shardPath, writeHeader);
        }
    }

    @Test
    public void testSetHeader() {
        final ShardingVCFWriter writer = createTestWriter(vcfBasePath);
        final VariantContext variant = createTestVariant(testVariantId);
        final VCFHeader header = new VCFHeader();
        header.setSequenceDictionary(dictionary);
        writer.setHeader(header);
        writer.add(variant);
        writer.close();
        assertVariantCount(1, firstShardPath, false);
    }
}