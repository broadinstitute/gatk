package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import com.intel.gkl.compression.IntelDeflaterFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class GatherVcfsIntegrationTest extends CommandLineProgramTest{

    private static final File TINY_HG_19_VCF = new File(hg19_chr1_1M_exampleVCF);
    private static final File LARGE_VCF = new File(largeFileTestDir, "gvcfs/combined.gatk3.7_30_ga4f720357.expected.vcf");


    @DataProvider
    public Object[][] getGatherTypes(){
        return new Object[][] {
            { GatherVcfs.GatherType.BLOCK },
            { GatherVcfs.GatherType.AUTOMATIC},
            {GatherVcfs.GatherType.CONVENTIONAL}
        };
    }

    @Test(groups = "bucket", dataProvider = "getGatherTypes")
    public void testGatherOverNIO(final GatherVcfs.GatherType gatherType){
        final File output = createTempFile("gathered", ".vcf.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput( getTestFile("cloud-inputs.list"))
            .addOutput(output)
            .addArgument(GatherVcfs.GATHER_TYPE_LONG_NAME, gatherType.toString());
        runCommandLine(args);

        final List<VariantContext> actual = Lists.newArrayList(new FeatureDataSource<VariantContext>(output));
        final List<VariantContext> expected = Lists.newArrayList(new FeatureDataSource<VariantContext>(
                new File(largeFileTestDir, "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf")));
        VariantContextTestUtils.assertEqualVariants(actual, expected);
    }

    @Test(dataProvider = "getGatherTypes")
    public void testGatherOverNIO(final GatherVcfs.GatherType gatherType){
        final File output = createTempFile("gathered", ".vcf.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput( getTestFile("bgzipped_shards.list"))
                .addOutput(output)
                .addArgument(GatherVcfs.GATHER_TYPE_LONG_NAME, gatherType.toString());
        runCommandLine(args);

        final List<VariantContext> actual = Lists.newArrayList(new FeatureDataSource<VariantContext>(output));
        final List<VariantContext> expected = Lists.newArrayList(new FeatureDataSource<VariantContext>(
                new File(, "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf")));
        VariantContextTestUtils.assertEqualVariants(actual, expected);
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
    
    @Test(dataProvider = "getVcfsToShard")
    public void testBlockGather(final File vcf, final int numShards) throws IOException {
        try (final FeatureDataSource<VariantContext> input = new FeatureDataSource<>(vcf)) {
            final ArrayList<VariantContext> expected = Lists.newArrayList(input);
            final List<File> shards = scatterVariants(expected, (VCFHeader) input.getHeader(), numShards,
                                                      createTempDir("vcfshards"));

            final File output = new File("failAgain.vcf.gz");//createTempFile("gathered", ".vcf.gz");
            final ArgumentsBuilder args = new ArgumentsBuilder();
            shards.forEach(args::addInput);

            args.addOutput(output)
                    .addArgument(StandardArgumentDefinitions.CLOUD_PREFETCH_BUFFER_LONG_NAME, "0")
                    .addArgument("gatherType", GatherVcfs.GatherType.BLOCK.toString())
                    .addBooleanArgument(GatherVcfs.IGNORE_SAFETY_CHECKS_LONG_NAME, true) // much faster this way
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
            assertSameNumberOfHeaderLines(vcf, output);
        }
    }

    private static void assertSameNumberOfHeaderLines(File vcf, File output) throws IOException {
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
                .addArgument(GatherVcfs.GATHER_TYPE_LONG_NAME, GatherVcfs.GatherType.BLOCK.toString())
                .addInput(TINY_HG_19_VCF);
        runCommandLine(args);
    }

}