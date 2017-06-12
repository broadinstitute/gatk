package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class GatherVcfsIntegrationTest extends CommandLineProgramTest{

    private static final File VCF = new File(largeFileTestDir, "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf");
    private static final File SMALL_GNOMAD_VCF = new File(largeFileTestDir, "very-small-gnomad.vcf");
    private static final File TINY_HG_19_VCF = new File(hg19_chr1_1M_exampleVCF);

    @Test(groups = "bucket")
    public void testConventionalGatherOverNIO(){
        final File output = createTempFile("gathered", ".vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput( getTestFile("cloud-inputs.list"))
            .addOutput(output)
            .addBooleanArgument("useConventionalGather", true);
        runCommandLine(args);

        final List<VariantContext> actual = Lists.newArrayList(new FeatureDataSource<VariantContext>(output));
        final List<VariantContext> expected = Lists.newArrayList(new FeatureDataSource<VariantContext>(
                VCF));
        VariantContextTestUtils.assertEqualVariants(actual, expected);
    }

    @DataProvider
    public Object[][] getVcfsToShard(){
        final File DB_SNP_VCF = new File(dbsnp_138_b37_20_21_vcf);
        return new Object[][] {
                {TINY_HG_19_VCF, 1},
                {TINY_HG_19_VCF, 2},
                {TINY_HG_19_VCF, 10},
                {DB_SNP_VCF, 1},
                {DB_SNP_VCF, 2},
                {DB_SNP_VCF, 10},
                {DB_SNP_VCF, 1000}
        };
    }

    @Test(dataProvider = "getVcfsToShard")
    public void testBlockGather(final File vcf, final int numShards) throws IOException {
        final FeatureDataSource<VariantContext> input = new FeatureDataSource<>( vcf );
        final ArrayList<VariantContext> expected = Lists.newArrayList(input);
        final List<File> shards = scatterVariants(expected, (VCFHeader) input.getHeader(), numShards,
                                                  createTempDir("vcfshards"));

        final File output = createTempFile("gathered", ".vcf.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        shards.forEach(args::addInput);

        args.addOutput(output)
                .addArgument("gatherType", GatherVcfs.GatherType.BLOCK.toString())
                .addBooleanArgument(GatherVcfs.IGNORE_SAFETY_CHECKS_LONG_NAME, true) // nomad file is missing a sequence dictionary
                .addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, false);
        runCommandLine(args);

        final FeatureDataSource<VariantContext> outputDataSource = new FeatureDataSource<>(output);
        final List<VariantContext> actual = Lists.newArrayList(outputDataSource);
        VariantContextTestUtils.assertEqualVariants(actual, expected);

        Assert.assertEquals(((VCFHeader) outputDataSource.getHeader()).getMetaDataInInputOrder(),
                            ((VCFHeader) input.getHeader()).getMetaDataInInputOrder());

        //test that the files have the same number of header lines as well as that the headers are identical
        //this catches a case we saw where header blocks were being inserted into the middle of the file
        try (final XReadLines actualReader = new XReadLines(output);
             final XReadLines expectedReader = new XReadLines(vcf)){
            final List<String> actualHeaderLines = getHeaderLines(actualReader);
            final List<String> expectedHeaderLines = getHeaderLines(expectedReader);
            Assert.assertEquals(actualHeaderLines.size() , expectedHeaderLines.size(),
                                "actual: " + String.join("\n",actualHeaderLines )+
                                        "\nexpected: " +String.join("\n", expectedHeaderLines ));
        }
    }

    private static List<String> getHeaderLines(final XReadLines actualReader) {
        return actualReader.readLines().stream().filter(line -> line.startsWith("#")).collect(Collectors.toList());
    }

    /**
     * Write variants into evenly split vcf.gz files and return the list of files
     */
    private static List<File> scatterVariants(List<VariantContext> variants, VCFHeader header, int shards, File shardDir){
        final List<List<VariantContext>> partitions = Lists.partition(variants, variants.size() / shards);
        return IntStream.range(0, partitions.size())
                .mapToObj(i -> writeShard(partitions.get(i), header, shardDir, i))
                .collect(Collectors.toList());
    }

    private static File writeShard(List<VariantContext> variants, VCFHeader header, File dir, int index){
        final File shard = new File( dir, + index + ".vcf.gz");
        try(final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(shard, header.getSequenceDictionary(), false)){
            writer.writeHeader(header);
            variants.forEach(writer::add);
        }
        return shard;
    }
    
}