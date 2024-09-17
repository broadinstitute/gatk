package org.broadinstitute.hellbender.engine;

import com.google.errorprone.annotations.Var;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.writers.IntervalFilteringVcfWriter;
import org.broadinstitute.hellbender.utils.variant.writers.ShardingVCFWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GatkToolIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Test
    public void testSitesOnlyMode() {
        File out = createTempFile("GTStrippedOutput", "vcf");
        String[] args = new String[] {
                "-V",  TEST_DIRECTORY + "vcf_with_genotypes.vcf",
                "--" + StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
                "-O",
                out.getAbsolutePath()};
        runCommandLine(Arrays.asList(args), SelectVariants.class.getSimpleName());

        // Assert that the genotype field has been stripped from the file
        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(out.getAbsolutePath());

        Assert.assertFalse(results.getLeft().hasGenotypingData());
        for (VariantContext v: results.getRight()) {
            Assert.assertFalse(v.hasGenotypes());
        }
    }

    @Test (expectedExceptions = java.lang.IllegalArgumentException.class)
    // test asserting that if the reference dictionary exists but is not valid we get a more helpful exception than a null pointer exception
    public void testBrokenReferenceDictionaryErrorMessage() throws IOException {
        File out = createTempFile("GTStrippedOutput", "vcf");

        Path refCopy = Files.copy(IOUtils.getPath(hg19_chr1_1M_Reference), createTempPath("reference", ".fasta"), StandardCopyOption.REPLACE_EXISTING);
        Path indexCopy = Files.copy(ReferenceSequenceFileFactory.getFastaIndexFileName(IOUtils.getPath(hg19_chr1_1M_Reference)), ReferenceSequenceFileFactory.getFastaIndexFileName(refCopy));
        File emptyDict = new File(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(refCopy).toString());
        IOUtils.deleteOnExit(indexCopy);

        emptyDict.createNewFile();
        IOUtils.deleteOnExit(emptyDict.toPath());

        String[] args = new String[] {
                "-R", refCopy.toString(),
                "-I", TEST_DIRECTORY + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam",
                "-O", out.getAbsolutePath()
        };

        runCommandLine(Arrays.asList(args), Mutect2.class.getSimpleName());
    }

    @Test
    public void testSharding() {
        final String outDir = createTempDir("GTShardedOutput").getAbsolutePath();
        final String fileBase = "test";
        final String out = Paths.get(outDir, fileBase + FileExtensions.COMPRESSED_VCF).toString();
        final String[] args = new String[] {
                "-V",  TEST_DIRECTORY + "example_variants_withSequenceDict.vcf",
                "-R", hg19MiniReference,
                "--" + StandardArgumentDefinitions.MAX_VARIANTS_PER_SHARD_LONG_NAME, "10",
                "-O", out};
        runCommandLine(Arrays.asList(args), SelectVariants.class.getSimpleName());

        // 12 total records in the test input should create 2 vcf shards with 10 and 2 records
        final Path basePath = Paths.get(outDir, fileBase).toAbsolutePath();
        final String firstShard = ShardingVCFWriter.getShardFilename(basePath, 0);
        final String secondShard = ShardingVCFWriter.getShardFilename(basePath, 1);
        final Pair<VCFHeader, List<VariantContext>> firstResults = VariantContextTestUtils.readEntireVCFIntoMemory(firstShard);
        final Pair<VCFHeader, List<VariantContext>> secondResults = VariantContextTestUtils.readEntireVCFIntoMemory(secondShard);
        Assert.assertEquals(firstResults.getValue().size(), 10, "First shard has wrong number of records");
        Assert.assertEquals(secondResults.getValue().size(), 2, "Second shard has wrong number of records");
        Assert.assertTrue(Files.exists(Paths.get(firstShard + FileExtensions.COMPRESSED_VCF_INDEX)));
        Assert.assertTrue(Files.exists(Paths.get(secondShard + FileExtensions.COMPRESSED_VCF_INDEX)));
    }

    @CommandLineProgramProperties(summary = "testTool which emits specific variants",
            oneLineSummary = "Test tool",
            programGroup = TestProgramGroup.class)
    public static class VariantEmitter extends GATKTool{
        static final SimpleInterval INT1 = new SimpleInterval("1",10, 15);
        static final SimpleInterval INT2 = new SimpleInterval("1",100, 105);
        static final SimpleInterval INT3 = new SimpleInterval("1",1000, 1005);
        static final SimpleInterval INT4 = new SimpleInterval("1",10000, 10005);
        static final SimpleInterval INT5 = new SimpleInterval("2",20, 25);
        static final SimpleInterval INT6 = new SimpleInterval("2",200, 205);
        static final SimpleInterval INT7 = new SimpleInterval("2",2000, 2005);
        static final SimpleInterval INT8 = new SimpleInterval("2",20000, 20005);
        static final SimpleInterval INT9_INS = new SimpleInterval("3",20000, 20000);

        static final List<Locatable> INTERVALS = List.of(INT1, INT2, INT3, INT4, INT5, INT6, INT7, INT8, INT9_INS);

        @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
        File output;

        @Override
        public void traverse() {
            //nope
        }
        @Argument(fullName = StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME,
                doc = "Restrict the output variants to ones that match the specified intervals according to the specified matching mode.",
                optional = true)
        @Advanced
        public IntervalFilteringVcfWriter.Mode outputVariantIntervalFilteringMode  = IntervalFilteringVcfWriter.Mode.ANYWHERE;

        @Override
        public IntervalFilteringVcfWriter.Mode getVariantFilteringOutputModeIfApplicable() {
            if (outputVariantIntervalFilteringMode != null) {
                return outputVariantIntervalFilteringMode;
            } else {
                return super.getVariantFilteringOutputModeIfApplicable();
            }
        }

        @Override
        public void onTraversalStart() {
            try(final VariantContextWriter vcfWriter = createVCFWriter(output)){
                vcfWriter.writeHeader(new VCFHeader());
                for(final Locatable interval : INTERVALS){
                    final VariantContextBuilder vcb = new VariantContextBuilder();
                    if (interval.getEnd()==interval.getStart()) {
                        vcb.alleles("A", "AAAAAAAAAA");
                        vcfWriter.add(vcb.chr(interval.getContig()).start(interval.getStart()).computeEndFromAlleles(Arrays.asList(Allele.create("A", true), Allele.create("AAAAAAAAAA")) ,interval.getStart()).make());
                    } else {
                        vcb.alleles("AAAAAA", "A").chr("1");
                        vcfWriter.add(vcb.loc(interval.getContig(),interval.getStart(), interval.getEnd()).make());
                    }
                }
            }
        }
    }

    @DataProvider
    public Object[][] getIntervalsAndOverlapMode(){
        final SimpleInterval chr1Interval = new SimpleInterval("1", 101, 10001);
        final SimpleInterval chr2Interval = new SimpleInterval("2", 201, 20001);
        final SimpleInterval chr1IntervalLeft = new SimpleInterval("1", 99, 102);
        final SimpleInterval chr1IntervalRight = new SimpleInterval("1", 103, 110);
        final SimpleInterval chr1IntervalNonAbutting = new SimpleInterval("1", 104, 110);

        final SimpleInterval chr1Interval99 = new SimpleInterval("1", 99, 110);
        final SimpleInterval chr1Interval100 = new SimpleInterval("1", 100, 110);
        final SimpleInterval chr1Interval101 = new SimpleInterval("1", 101, 110);

        final SimpleInterval chr3Interval99 = new SimpleInterval("3", 19999, 19999);
        final SimpleInterval chr3Interval100 = new SimpleInterval("3", 20000, 20000);
        final SimpleInterval chr3Interval101 = new SimpleInterval("3", 20001, 20001);

        return new Object[][]{
                {Arrays.asList(chr1Interval, chr2Interval), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.ANYWHERE, VariantEmitter.INTERVALS },
                {Arrays.asList(chr1Interval, chr2Interval), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.OVERLAPS, List.of(VariantEmitter.INT2, VariantEmitter.INT3, VariantEmitter.INT4, VariantEmitter.INT6, VariantEmitter.INT7, VariantEmitter.INT8)},
                {Arrays.asList(chr1Interval, chr2Interval), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, List.of(VariantEmitter.INT3, VariantEmitter.INT4, VariantEmitter.INT7,  VariantEmitter.INT8)},
                {Arrays.asList(chr1Interval, chr2Interval), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.ENDS_IN, List.of(VariantEmitter.INT2, VariantEmitter.INT3, VariantEmitter.INT6, VariantEmitter.INT7)},
                {Arrays.asList(chr1Interval, chr2Interval), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.CONTAINED, List.of(VariantEmitter.INT3, VariantEmitter.INT7)},

                // Tests specifically aimed at documenting how the --interval-merging-rule argument works in conjunction with the interval filtering
                {Arrays.asList(chr1IntervalLeft, chr1IntervalRight), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.CONTAINED, List.of(VariantEmitter.INT2)}, // Default is to merge all
                {Arrays.asList(chr1IntervalLeft, chr1IntervalRight), Arrays.asList("--interval-merging-rule", "ALL"), IntervalFilteringVcfWriter.Mode.CONTAINED, List.of(VariantEmitter.INT2)},
                {Arrays.asList(chr1IntervalLeft, chr1IntervalRight), Arrays.asList("--interval-merging-rule", "OVERLAPPING_ONLY"), IntervalFilteringVcfWriter.Mode.CONTAINED, new ArrayList<>()},
                {Arrays.asList(chr1IntervalLeft, chr1IntervalNonAbutting), Arrays.asList("--interval-merging-rule", "OVERLAPPING_ONLY"), IntervalFilteringVcfWriter.Mode.CONTAINED, new ArrayList<>()},
                {Arrays.asList(chr1IntervalLeft, chr1IntervalNonAbutting), Arrays.asList("--interval-merging-rule", "OVERLAPPING_ONLY", "--interval-padding", "10"), IntervalFilteringVcfWriter.Mode.CONTAINED, List.of(VariantEmitter.INT2)},

                // Demonstrating the exact behavior of the starts_in/ends_in modes
                {Arrays.asList(chr1Interval99), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, List.of(VariantEmitter.INT2)},
                {Arrays.asList(chr1Interval100), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, List.of(VariantEmitter.INT2)}, //deletion where left base is at 100
                {Arrays.asList(chr1Interval101), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, new ArrayList<>()},

                // Deomstrating the behavior for starts/ends_in with an insertion
                {Arrays.asList(chr3Interval99), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, new ArrayList<>()},
                {Arrays.asList(chr3Interval100), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, List.of(VariantEmitter.INT9_INS)},
                {Arrays.asList(chr3Interval101), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.STARTS_IN, new ArrayList<>()},
                {Arrays.asList(chr3Interval99), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.ENDS_IN, new ArrayList<>()},
                {Arrays.asList(chr3Interval100), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.ENDS_IN, List.of(VariantEmitter.INT9_INS)},
                {Arrays.asList(chr3Interval101), new ArrayList<>(), IntervalFilteringVcfWriter.Mode.ENDS_IN, new ArrayList<>()},
        };
    }

    @Test(dataProvider = "getIntervalsAndOverlapMode")
    public void testVcfOutputFilterMode(List<? extends Locatable> intervals, List<String> extraArguments, IntervalFilteringVcfWriter.Mode mode, List<Locatable> expected){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final File out = createTempFile("out", ".vcf");
        args.addOutput(out);
        intervals.forEach(args::addInterval);
        args.addReference(b37Reference);
        extraArguments.forEach(args::addRaw);
        if( mode != null) {
            args.add(StandardArgumentDefinitions.VARIANT_OUTPUT_INTERVAL_FILTERING_MODE_LONG_NAME, mode);
        }

        runCommandLine(args, VariantEmitter.class.getSimpleName());
        final Pair<VCFHeader, List<VariantContext>> vcfHeaderListPair = VariantContextTestUtils.readEntireVCFIntoMemory(out.toString());

        final List<VariantContext> actual = vcfHeaderListPair.getRight();
        Assert.assertEquals(actual.size(), expected.size());
        BaseTest.assertCondition(actual, expected, (left, right) -> {
            Assert.assertEquals(left.getContig(), right.getContig());
            Assert.assertEquals(left.getStart(), right.getStart());
            Assert.assertEquals(left.getEnd(), right.getEnd());
        } );

    }
}
