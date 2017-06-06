package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class GatherVcfsIntegrationTest extends CommandLineProgramTest{

    private static final File VCF = new File(largeFileTestDir, "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf");

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
        return new Object[][] {
                {new File(hg19_chr1_1M_exampleVCF), 1},
                {new File(hg19_chr1_1M_exampleVCF), 2},
                {new File(hg19_chr1_1M_exampleVCF), 10},
                {new File(largeFileTestDir, "very-small-gnomad.vcf"), 10}
        };
    }

    @Test(dataProvider = "getVcfsToShard")
    public void testBlockGather(final File vcf, final int numShards){
        final FeatureDataSource<VariantContext> input = new FeatureDataSource<>( vcf );
        final ArrayList<VariantContext> expected = Lists.newArrayList(input);
        final List<File> shards = scatterVariants(expected, (VCFHeader) input.getHeader(), numShards);

        final File output = createTempFile("gathered", "vcf.gz");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        shards.forEach(args::addInput);

        args.addOutput(output)
                .addBooleanArgument(GatherVcfs.IGNORE_SAFETY_CHECKS_LONG_NAME, true) // nomad file is missing a sequence dictionary
                .addBooleanArgument(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, false);
        runCommandLine(args);

        final List<VariantContext> actual = Lists.newArrayList(new FeatureDataSource<VariantContext>(output));
        VariantContextTestUtils.assertEqualVariants(actual, expected);
    }

    /**
     * Write variants into evenly split vcf.gz files and return the list of files
     */
    private static List<File> scatterVariants(List<VariantContext> variants, VCFHeader header, int shards){
        final List<List<VariantContext>> partitions = Lists.partition(variants, shards);
        final File shardDir = createTempDir("vcfshards");
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