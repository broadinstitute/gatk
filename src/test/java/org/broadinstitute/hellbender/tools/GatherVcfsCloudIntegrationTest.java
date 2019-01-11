package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class GatherVcfsCloudIntegrationTest extends CommandLineProgramTest{

    private static final File TINY_HG_19_VCF = new File(hg19_chr1_1M_exampleVCF);
    private static final File LARGE_VCF = new File(largeFileTestDir, "gvcfs/combined.gatk3.7_30_ga4f720357.expected.vcf");

    @DataProvider
    public Object[][] getGatherTypes(){
        return new Object[][] {
                { GatherVcfsCloud.GatherType.BLOCK },
                { GatherVcfsCloud.GatherType.AUTOMATIC},
                { GatherVcfsCloud.GatherType.CONVENTIONAL}
        };
    }

    @Test(groups = "bucket", dataProvider = "getGatherTypes")
    public void testGatherOverNIO(final GatherVcfsCloud.GatherType gatherType) throws IOException {
        assertGatherProducesCorrectVariants(gatherType,
                new File(largeFileTestDir, "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf"),
                getTestFile("cloud-inputs.args"));
    }

    @Test(dataProvider = "getGatherTypes")
    public void testLocalGathers(final GatherVcfsCloud.GatherType gatherType) throws IOException {
        assertGatherProducesCorrectVariants(gatherType,
                getTestFile("gzipped.vcf.gz"),
                getTestFile("bgzipped_shards.args"));
    }

    private void assertGatherProducesCorrectVariants(GatherVcfsCloud.GatherType gatherType, File expected, File inputs) throws IOException {
        assertGatherProducesCorrectVariants(gatherType, expected, Collections.singletonList(inputs));
    }

    @SuppressWarnings({"unchecked"})
    private void assertGatherProducesCorrectVariants(GatherVcfsCloud.GatherType gatherType, File expected, List<File> inputs) throws IOException {
        final File output = createTempFile("gathered", ".vcf.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        inputs.forEach(args::addInput);
        args.addOutput(output)
                .addArgument(GatherVcfsCloud.GATHER_TYPE_LONG_NAME, gatherType.toString());
        runCommandLine(args);

        try (final AbstractFeatureReader<VariantContext, ?> actualReader = AbstractFeatureReader.getFeatureReader(output.getAbsolutePath(), null, new VCFCodec(), false, Function.identity(), Function.identity());
             final AbstractFeatureReader<VariantContext, ?> expectedReader = AbstractFeatureReader.getFeatureReader(expected.getAbsolutePath(), null, new VCFCodec(), false, Function.identity(), Function.identity())) {

            final List<VariantContext> actualVariants = StreamSupport.stream(Spliterators.spliteratorUnknownSize(actualReader.iterator(),Spliterator.ORDERED), false).collect(Collectors.toList());
            final List<VariantContext> expectedVariants = StreamSupport.stream(Spliterators.spliteratorUnknownSize(expectedReader.iterator(),Spliterator.ORDERED), false).collect(Collectors.toList());
            VariantContextTestUtils.assertEqualVariants(actualVariants, expectedVariants);

            Assert.assertEquals(((VCFHeader) actualReader.getHeader()).getMetaDataInInputOrder(),
                    ((VCFHeader) expectedReader.getHeader()).getMetaDataInInputOrder());
        }
    }

    @DataProvider
    public Object[][] getHeaderOnlyCases(){
        final File onlyHeader = getTestFile("onlyHeader.vcf.gz");
        final File shard0 = getTestFile("shard_0.vcfs.gz");
        return new Object[][]{
                {Collections.singletonList(onlyHeader), GatherVcfsCloud.GatherType.BLOCK, onlyHeader},
                {Collections.singletonList(onlyHeader), GatherVcfsCloud.GatherType.CONVENTIONAL, onlyHeader},
                {Arrays.asList(onlyHeader, onlyHeader), GatherVcfsCloud.GatherType.BLOCK, onlyHeader},
                {Arrays.asList(onlyHeader, onlyHeader), GatherVcfsCloud.GatherType.CONVENTIONAL, onlyHeader},
                {Arrays.asList(onlyHeader, shard0), GatherVcfsCloud.GatherType.BLOCK, shard0},
                {Arrays.asList(onlyHeader, shard0), GatherVcfsCloud.GatherType.CONVENTIONAL, shard0},
                {Arrays.asList(shard0, onlyHeader), GatherVcfsCloud.GatherType.BLOCK, shard0},
                {Arrays.asList(shard0, onlyHeader), GatherVcfsCloud.GatherType.CONVENTIONAL, shard0},
        };
    }

    @Test(dataProvider = "getHeaderOnlyCases")
    public void testVcfsWithOnlyAHeader(List<File> inputs, final GatherVcfsCloud.GatherType gatherType, File expected) throws IOException {
        assertGatherProducesCorrectVariants(gatherType, expected, inputs);
    }

    @DataProvider
    public Object[][] getVcfsToShard(){
        final File db_snpVCF = new File(dbsnp_138_b37_20_21_vcf);
        return new Object[][] {
                //the shard number is very silly, vcf.size / shards is used as the partition size,
                //since it's integer division you won't get the exact number of shards you requested, but larger
                //numbers will produce more shards
                {TINY_HG_19_VCF, 1},
                {TINY_HG_19_VCF, 2},
                {TINY_HG_19_VCF, 10},
                {db_snpVCF, 1},
                {db_snpVCF, 2},
                {db_snpVCF, 10},
                {db_snpVCF, 100},
                {db_snpVCF, 1000},
                {LARGE_VCF, 1},
                {LARGE_VCF, 2},
                {LARGE_VCF, 5},
                {LARGE_VCF, 10},
                {LARGE_VCF, 100},
                {LARGE_VCF, 1000},
        };
    }

    @Test(dataProvider = "getVcfsToShard", enabled = false)
    public void testBlockGather(final File vcf, final int numShards) throws IOException {
        try (final FeatureDataSource<VariantContext> input = new FeatureDataSource<>(vcf)) {
            final ArrayList<VariantContext> expected = Lists.newArrayList(input);
            final List<File> shards = scatterVariants(expected, (VCFHeader) input.getHeader(), numShards,
                    createTempDir("vcfshards"));

            final File output = createTempFile("testBlockGather_gathered", ".vcf.gz");
            final ArgumentsBuilder args = new ArgumentsBuilder();
            shards.forEach(args::addInput);

            args.addOutput(output)
                    .addArgument(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, "0")
                    .addArgument(GatherVcfsCloud.GATHER_TYPE_LONG_NAME, GatherVcfsCloud.GatherType.BLOCK.toString())
                    .addBooleanArgument(GatherVcfsCloud.IGNORE_SAFETY_CHECKS_LONG_NAME, true) // much faster this way
                    .addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, false);
            runCommandLine(args);

            try (final FeatureDataSource<VariantContext> outputDataSource = new FeatureDataSource<>(output)) {
                final List<VariantContext> actual = Lists.newArrayList(outputDataSource);
                VariantContextTestUtils.assertEqualVariants(actual, expected);

                Assert.assertEquals(((VCFHeader) outputDataSource.getHeader()).getMetaDataInInputOrder(),
                        ((VCFHeader) input.getHeader()).getMetaDataInInputOrder());
            }

            //test that the files have the same number of header lines as well as that the headers are identical
            //this catches a case we saw where header blocks were being inserted into the middle of the file
            assertFilesHaveSameNumberOfHeaderLines(vcf, output);
        }
    }

    private static void assertFilesHaveSameNumberOfHeaderLines(File vcf, File output) throws IOException {
        //test that the files have the same number of header lines as well as that the headers are identical
        //this catches a case we saw where header blocks were being inserted into the middle of the file
        try (final XReadLines actualReader = new XReadLines(output);
             final XReadLines expectedReader = new XReadLines(vcf)) {
            final List<String> actualHeaderLines = getHeaderLines(actualReader);
            final List<String> expectedHeaderLines = getHeaderLines(expectedReader);
            Assert.assertEquals(actualHeaderLines.size(), expectedHeaderLines.size(),
                    "actual: " + String.join("\n", actualHeaderLines) +
                            "\nexpected: " + String.join("\n", expectedHeaderLines));
        }
    }

    private static List<String> getHeaderLines(final XReadLines actualReader) {
        return actualReader.readLines().stream()
                .filter(line -> line.startsWith("#"))
                .collect(Collectors.toList());
    }

    /**
     * Write variants into evenly split vcf.gz files and return the list of files
     *
     * note that shards is used in a very silly way and may not do what you're expecting
     */
    private static List<File> scatterVariants(final List<VariantContext> variants, final VCFHeader header, final int shards, final File shardDir){
        final List<List<VariantContext>> partitions = Lists.partition(variants, variants.size() / shards);
        final List<File> list = new ArrayList<>();
        for (int i = 0; i < partitions.size(); i++) {
            final File file = writeShard(partitions.get(i), header, shardDir, i);
            try(final FeatureDataSource<VariantContext> features = new FeatureDataSource<>(file)) {
                VariantContextTestUtils.assertEqualVariants(Lists.newArrayList(features),       partitions.get(i));

            }
            list.add(file);
        }

        return list;
    }

    private static File writeShard(final List<VariantContext> variants, final VCFHeader header, final File dir, final int index){
        final File shard = new File( dir, + index + ".vcf.gz");
        try(final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(shard, header.getSequenceDictionary(), false)){
            writer.writeHeader(header);
            variants.forEach(writer::add);
        }
        return shard;
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBlockCopyFailsIfNonBGZipped(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(createTempFile("output","vcf.gz"))
                .addArgument(GatherVcfsCloud.GATHER_TYPE_LONG_NAME, GatherVcfsCloud.GatherType.BLOCK.toString())
                .addInput(TINY_HG_19_VCF);
        runCommandLine(args);
    }

}