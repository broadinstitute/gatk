package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.VCFHeaderReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import org.apache.commons.codec.digest.DigestUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBArgumentCollection;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GenotypeGVCFsIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test --tests org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFsIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    private static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();
    private static final String BASE_PAIR_EXPECTED = "gvcf.basepairResolution.gatk3.7_30_ga4f720357.output.vcf";
    private static final String b38_reference_20_21 = largeFileTestDir + "Homo_sapiens_assembly38.20.21.fasta";
    private static final String BASE_PAIR_GVCF = "gvcf.basepairResolution.gvcf";
    private static final double TLOD_THRESHOLD = 4.0;

    private static final File CEUTRIO_20_21_GATK3_4_G_VCF = new File(largeFileTestDir, "gvcfs/CEUTrio.20.21.gatk3.4.g.vcf");
    private static final String CEUTRIO_20_21_EXPECTED_VCF = "CEUTrio.20.21.gatk3.7_30_ga4f720357.expected.vcf";
    private static final File NA12878_HG37 = new File(toolsTestDir + "GenomicsDBImport/expected.testGVCFMode.gatk4.g.vcf");
    private static final List<String> ATTRIBUTES_WITH_JITTER = Arrays.asList(
            "AS_QD",
            "QD",//TODO QD and AS_QD have cap values and anything that reaches that is randomized.  It's difficult to reproduce the same random numbers across gatk3 -> 4
            "FS");//TODO There's some bug in either gatk3 or gatk4 fisherstrand that's making them not agree still, I'm not sure which is correct
    private static final List<String> ATTRIBUTES_TO_IGNORE = Arrays.asList("FS","RAW_MQ","RGQ","MQ"); //MQ data format and key have changed since GATK3

    private static final String ALLELE_SPECIFIC_DIRECTORY = toolsTestDir + "walkers/annotator/allelespecific";

    private static <T> void assertForEachElementInLists(final List<T> actual, final List<T> expected, final BiConsumer<T, T> assertion) {
        Assert.assertEquals(actual.size(), expected.size(), "different number of elements in lists:\n"
                + actual.stream().map(Object::toString).collect(Collectors.joining("\n","actual:\n","\n"))
            +  expected.stream().map(Object::toString).collect(Collectors.joining("\n","expected:\n","\n")));
        for (int i = 0; i < actual.size(); i++) {
            assertion.accept(actual.get(i), expected.get(i));
        }
    }

    @DataProvider(name = "gvcfsToGenotype")
    public Object[][] gvcfsToGenotype() {
        return new Object[][]{
                //combine not supported yet, see https://github.com/broadinstitute/gatk/issues/2429 and https://github.com/broadinstitute/gatk/issues/2584
                //{"combine.single.sample.pipeline.1.vcf", null, Arrays.asList("-V", getTestFile("combine.single.sample.pipeline.2.vcf").toString() , "-V", getTestFile("combine.single.sample.pipeline.3.vcf").toString()), b37_reference_20_21},
                {getTestFile("leadingDeletion.g.vcf"), getTestFile("leadingDeletionRestrictToStartExpected.vcf"), Arrays.asList("-L", "20:69512-69513", "--"+GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME), b37_reference_20_21},
                {getTestFile("leadingDeletion.g.vcf"), getTestFile("leadingDeletionExpected.vcf"), Arrays.asList("-L", "20:69512-69513"), b37_reference_20_21},
                {getTestFile(BASE_PAIR_GVCF), getTestFile( BASE_PAIR_EXPECTED), NO_EXTRA_ARGS, b37_reference_20_21}, //base pair level gvcf
                {getTestFile("testUpdatePGT.gvcf"), getTestFile( "testUpdatePGT.gatk3.7_30_ga4f720357.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},   //testUpdatePGT
                {getTestFile("gvcfExample1.vcf"), getTestFile( "gvcfExample1.gatk3.7_30_ga4f720357.expected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //single sample vcf
                {getTestFile("gvcfExample1.vcf"), getTestFile( "gvcfExample1.gatk3.7_30_ga4f720357.expected.vcf"), Arrays.asList("-L", "20"), b37_reference_20_21}, //single sample vcf with -L
                {getTestFile("combined_genotype_gvcf_exception.original.vcf"), getTestFile( "combined_genotype_gvcf_exception.gatk3.7_30_ga4f720357.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //test that an input vcf with 0/0 already in GT field is overwritten
                {getTestFile("combined_genotype_gvcf_exception.nocall.vcf"), getTestFile( "combined_genotype_gvcf_exception.gatk3.7_30_ga4f720357.output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},  //same test as above but with ./.
                {getTestFile(BASE_PAIR_GVCF), getTestFile("ndaTest.gatk3.7_30_ga4f720357.expected.vcf"), Collections.singletonList("--annotate-with-num-discovered-alleles"), b37_reference_20_21},  //annotating with the number of alleles discovered option
                {getTestFile(BASE_PAIR_GVCF), getTestFile("maxAltAllelesTest.gatk3.7_30_ga4f720357.expected.vcf"), Arrays.asList("--max-alternate-alleles", "1"), b37_reference_20_21 }, //restricting the max number of alt alleles
                {getTestFile(BASE_PAIR_GVCF), getTestFile("standardConfTest.gatk3.7_30_ga4f720357.expected.vcf"), Arrays.asList("-stand-call-conf", "300"),b37_reference_20_21}, //set minimum calling threshold
                {getTestFile("spanningDel.combined.g.vcf"), getTestFile( "spanningDel.combined.gatk3.7_30_ga4f720357.expected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                {getTestFile("spanningDel.delOnly.g.vcf"), getTestFile( "spanningDel.delOnly.gatk3.7_30_ga4f720357.expected.vcf"), NO_EXTRA_ARGS, b37_reference_20_21},
                {getTestFile("spanningDel.depr.delOnly.g.vcf"), getTestFile( "spanningDel.depr.delOnly.gatk3.7_30_ga4f720357.expected.vcf" ), NO_EXTRA_ARGS, b37_reference_20_21},
                {getTestFile("ad-bug-input.vcf"), getTestFile( "ad-bug-gatk3.7_30_ga4f720357-output.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}, //Bad AD Propagation Haploid Bug
                {CEUTRIO_20_21_GATK3_4_G_VCF, getTestFile(CEUTRIO_20_21_EXPECTED_VCF), Arrays.asList("--dbsnp", largeFileTestDir + "dbsnp_138.b37.20.21.vcf"), b37_reference_20_21},
                {getTestFile("CEUTrio.20.21.missingIndel.g.vcf"), getTestFile( "CEUTrio.20.21.missingIndel.gatk3.7_30_ga4f720357.expected.vcf"), Arrays.asList("--dbsnp", "src/test/resources/large/dbsnp_138.b37.20.21.vcf"), b37_reference_20_21},
                {new File(largeFileTestDir + "gvcfs/gatk3.7_30_ga4f720357.24_sample.21.g.vcf"), new File( largeFileTestDir + "gvcfs/gatk3.7_30_ga4f720357.24_sample.21.expected.vcf"), NO_EXTRA_ARGS, b38_reference_20_21},
                {getTestFile("chr21.bad.pl.g.vcf"), getTestFile( "chr21.bad.pl.gatk3.7_30_ga4f720357.expected.vcf"), Arrays.asList("-L", "chr21:28341770-28341790"), b38_reference_20_21},
                {new File(largeFileTestDir + "gvcfs/combined.gatk3.7_30_ga4f720357.g.vcf.gz"),  new File(largeFileTestDir + "gvcfs/combined.gatk3.7_30_ga4f720357.expected.vcf"), NO_EXTRA_ARGS, b38_reference_20_21},
                //Tests for Allele-Specific Annotations
                {new File(ALLELE_SPECIFIC_DIRECTORY, "NA12878.AS.chr20snippet.g.vcf"), getTestFile( "AS_Annotations.gatk3.7_30_ga4f720357.expected.vcf"), Arrays.asList( "-A", "ClippingRankSumTest", "-G", "AS_StandardAnnotation", "-G", "StandardAnnotation"), b37_reference_20_21},
                {new File(ALLELE_SPECIFIC_DIRECTORY, "NA12878.AS.chr20snippet.g.vcf"), getTestFile( "AS_Annotations.keepRawCombined.expected.vcf"), Arrays.asList( "-A", "ClippingRankSumTest", "-G", "AS_StandardAnnotation", "-G", "StandardAnnotation", "-keep-combined"), b37_reference_20_21},
                //input GVCF doesn't have AS annotations, just new RAW_MQandDP format; expected output annotations include AS_QD because that can be calculated from genotypes; definitely should have MQ and retain F1R2:F2R1 for OxoG in FORMAT
                {getTestFile( "withOxoGReadCounts.g.vcf"), getTestFile( "withOxoGReadCounts.vcf"), Arrays.asList("-G", "AS_StandardAnnotation", "-G", "StandardAnnotation"), b37_reference_20_21},
                {getTestFile( "multiSamples.g.vcf"), getTestFile( "multiSamples.GATK3expected.g.vcf"), Arrays.asList( "-A", "ClippingRankSumTest", "-G", "AS_StandardAnnotation", "-G", "StandardAnnotation"), b37_reference_20_21},
                {getTestFile( "testAlleleSpecificAnnotations.CombineGVCF.output.g.vcf"), getTestFile( "testAlleleSpecificAnnotations.CombineGVCF.expected.g.vcf"), Arrays.asList( "-A", "ClippingRankSumTest", "-G", "AS_StandardAnnotation", "-G", "StandardAnnotation"), b37_reference_20_21},

                // all sites/--include-non-variant-sites tests
                // The results from these tests differ from GATK3 in the following ways:
                //  - sites where the only alternate allele is a spanning deletion are emitted by GATK3, but not emitted by GATK4
                //  - LowQual variants are not emitted by GATK3, but are emitted by GATK4
                //  - GATK3 added `AN` annotations to non-variant sites, but GATK3 does not
                {getTestFile(BASE_PAIR_GVCF), getTestFile( "expected/gvcf.basepairResolution.includeNonVariantSites.vcf"), Collections.singletonList("--" + GenotypeGVCFs.ALL_SITES_LONG_NAME), b37_reference_20_21 },
                {getTestFile( "combine.single.sample.pipeline.1.vcf"),
                        getTestFile( "expected/combine.single.sample.pipeline.1.include_nonvariant.vcf"),
                        Arrays.asList( " --" + GenotypeGVCFs.ALL_SITES_LONG_NAME + " -L 20:10,030,000-10,033,000 -L 20:10,386,000-10,386,500 "),
                        b37_reference_20_21},
                {getTestFile( "combine.single.sample.pipeline.2.vcf"),
                        getTestFile( "expected/combine.single.sample.pipeline.2.include_nonvariant.vcf"),
                        Arrays.asList( " --" + GenotypeGVCFs.ALL_SITES_LONG_NAME + " -L 20:10,030,000-10,033,000 -L 20:10,386,000-10,386,500 "),
                        b37_reference_20_21},
                {getTestFile( "combine.single.sample.pipeline.3.vcf"),
                        getTestFile( "expected/combine.single.sample.pipeline.3.include_nonvariant.vcf"),
                        Arrays.asList( " --" + GenotypeGVCFs.ALL_SITES_LONG_NAME + " -L 20:10,030,000-10,033,000 -L 20:10,386,000-10,386,500 "),
                        b37_reference_20_21},
                // combined, with intervals
                {getTestFile( "combined.single.sample.pipeline.gatk3.vcf"),
                        getTestFile( "expected/combined.single.sample.pipeline.include_nonvariant.vcf"),
                        Arrays.asList( " --" + GenotypeGVCFs.ALL_SITES_LONG_NAME + " -L 20:10,030,000-10,033,000 -L 20:10,386,000-10,386,500 "),
                        b37_reference_20_21},
                // test site 10096905 - 10096907 to force coverage around a spanning deletion only site, and 20:10624924-1062492 to
                // force coverage around a multi-allelic variant that includes a spanning deletion
                {getTestFile( "combined.single.sample.pipeline.gatk3.vcf"),
                        getTestFile( "expected/testSpanningDeletion.vcf"),
                        Arrays.asList( " --" + GenotypeGVCFs.ALL_SITES_LONG_NAME + " -L 20:10,096,905-10,096,907 -L 20:10624924-10624926"),
                        b37_reference_20_21},
                // test site 20:10,012,730-10,012,740 to force coverage around LowQual site
                {getTestFile( "combined.single.sample.pipeline.gatk3.vcf"),
                        getTestFile( "expected/includeLowQualSites.vcf"),
                        Arrays.asList( " --" + GenotypeGVCFs.ALL_SITES_LONG_NAME + " -L 20:10,012,730-10,012,740"),
                        b37_reference_20_21}
        };
    }

    @DataProvider(name = "singleSampleGVCFWithNewMQFormat")
    public Object[][] singleSampleGVCFWithNewMQFormat() {
        return new Object[][]{
                {NA12878_HG37, getTestFile("newMQcalc.singleSample.genotyped.vcf"), new SimpleInterval("20", 1, 11_000_000), b37_reference_20_21}
        };
    }

    @DataProvider(name = "GVCFsWithNewMQFormat")
    public Object[][] GVCFsWithNewMQFormat() {
        return new Object[][]{
                {NA12878_HG37, getTestFile("newMQcalc.singleSample.genotyped.vcf"), new SimpleInterval("20", 1, 11_000_000), b37_reference_20_21},
                {new File(getTestDataDir() + "/walkers/CombineGVCFs/newMQcalc.combined.g.vcf"), getTestFile("newMQcalc.combined.genotyped.vcf"), new SimpleInterval("20", 1, 11_000_000), b37_reference_20_21}
        };
    }

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    /*
    This test is useful for testing changes in GATK4 versus different versions of GATK3.
    To use, set GATK3_PATH to point to a particular version of gatk, then enable this test and run.

    It will cache the gatk3 outputs in a folder called gatk3results, it does it's best to avoid reusing bad results by
    comparing the md5 of the gatk3 path, input file path, reference, and commandline, but it doesn't know about internal changes to files
    You have to manually delete the cache if you make changes inside the input files.

    The expected outputs based on gatk3's results will be put in a folder called expectedResults.
    These are overwritten during each test, the names are based on the name of the existing expected output file.

    This method should be removed after GenotypeGVCFs has been completely validated against GATK3.
     */
    @Test(dataProvider = "gvcfsToGenotype", enabled = false)
    public void compareToGATK3(File input, File outputFile, List<String> extraArgs, String reference) throws IOException, NoSuchAlgorithmException {
        final String GATK3_PATH = "gatk3.7-30-ga4f720357.jar";
        final String params = GATK3_PATH + input.getAbsolutePath() + extraArgs.stream().collect(Collectors.joining()) + reference;
        final String md5 = DigestUtils.md5Hex(params);
        final File gatk3ResultsDir = new File("gatk3results");
        if(! gatk3ResultsDir.exists()){
            Assert.assertTrue(gatk3ResultsDir.mkdir());
        }
        final File gatk3Result = new File(gatk3ResultsDir, md5 + ".vcf");
        if (!gatk3Result.exists()) {
            List<String> gatk3Command = new ArrayList<>(
                    Arrays.asList("java", "-jar", GATK3_PATH, "-T", "GenotypeGVCFs"));
            gatk3Command.add("-V");
            gatk3Command.add(input.getAbsolutePath());
            gatk3Command.add("-o");
            gatk3Command.add(gatk3Result.getAbsolutePath());
            gatk3Command.add("-R");
            gatk3Command.add(reference);
            gatk3Command.addAll(extraArgs);

            runProcess(new ProcessController(), gatk3Command.toArray(new String[gatk3Command.size()]));
        } else {
            System.out.println("Found precomputed gatk3Result");
        }
        final Path expectedResultsDir = Paths.get("expectedResults");
        if ( !Files.exists(expectedResultsDir)) {
            Files.createDirectory(expectedResultsDir);
        }
        Files.copy(gatk3Result.toPath(), expectedResultsDir.resolve(outputFile.getName()), StandardCopyOption.REPLACE_EXISTING);

        assertGenotypesMatch(input, gatk3Result, extraArgs, reference);
        assertVariantContextsMatch(input, gatk3Result, extraArgs, reference);
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testGenotypesOnly(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        final List<String> extra = new ArrayList<>(extraArgs);
        assertGenotypesMatch(input, expected, extra, reference);
    }

    @DataProvider
    //this is different from the above data provider because we can currently only load a single interval into a genomics db in a sane way
    //so we need to provide a list of intervals and then look at each one
    public Object[][] getGVCFsForGenomicsDB(){
        return new Object[][]{
                {getTestFile(BASE_PAIR_GVCF), getTestFile(BASE_PAIR_EXPECTED), new SimpleInterval("20", 1, 11_000_000), b37_reference_20_21},
                {CEUTRIO_20_21_GATK3_4_G_VCF, getTestFile("CEUTrio.20.gatk3.7_30_ga4f720357.expected.vcf"), new SimpleInterval("20", 1, 11_000_000), b37_reference_20_21},
                {CEUTRIO_20_21_GATK3_4_G_VCF, getTestFile("CEUTrio.21.gatk3.7_30_ga4f720357.expected.vcf"), new SimpleInterval("21", 1, 11_000_000), b37_reference_20_21}
        };
    }

    @DataProvider
    public Object[][] getGVCFsForGenomicsDBOverMultipleIntervals() {
        LinkedList<SimpleInterval> intervals = new LinkedList<SimpleInterval>();
        //[ 10000117, 10020107 ]
        int base = 10000117;
        for (int i = 0; i < 1000; ++i)
            intervals.add(new SimpleInterval("20", base + 20 * i, base + 20 * i + 10)); //intervals of size 10 separated by 10

        return new Object[][]{
                {CEUTRIO_20_21_GATK3_4_G_VCF, getTestFile("NA12878.mergedIntervals.vcf"), intervals, b37_reference_20_21}
        };
    }

    //this only tests single-sample
    @Test(dataProvider = "getGVCFsForGenomicsDB", timeOut = 1000000)
    public void assertMatchingGenotypesFromGenomicsDB(File input, File expected, Locatable interval, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(input, interval);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);
        runGenotypeGVCFSAndAssertSomething(genomicsDBUri, expected, NO_EXTRA_ARGS, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes, reference);

        // The default option with GenomicsDB input uses VCFCodec for decoding, test BCFCodec explicitly
        final List<String> args = new ArrayList<String>();
        args.add("--"+GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME);
        runGenotypeGVCFSAndAssertSomething(genomicsDBUri, expected, args, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes, reference);
    }

    private void runAndCheckGenomicsDBOutput(final ArgumentsBuilder args, final File expected, final File output) {
        Utils.resetRandomGenerator();
        runCommandLine(args);

        // Note that if this isn't working it will take *FOREVER*
        // runs in 0.06 minutes with no input intervals specfied
        final List<VariantContext> expectedVC = VariantContextTestUtils.getVariantContexts(expected);
        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes);
    }

    @Test(dataProvider = "getGVCFsForGenomicsDBOverMultipleIntervals")
    public void testGenotypeGVCFsMultiIntervalGDBQuery(File input, File expected, List<Locatable> intervals, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(input, intervals, true);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .add("V", genomicsDBUri);
        args.addOutput(output);
        intervals.forEach(args::addInterval);
        args.addRaw("--" + GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME);
        args.addRaw("--only-output-calls-starting-in-intervals");  //note that this will restrict calls to just the specified intervals

        runAndCheckGenomicsDBOutput(args, expected, output);

        // The default option with GenomicsDB input uses VCFCodec for decoding, test BCFCodec explicitly
        args.addRaw("--"+GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME);
        runAndCheckGenomicsDBOutput(args, expected, output);
    }

    //this tests single-sample with new MQ format
    @Test (dataProvider = "singleSampleGVCFWithNewMQFormat")
    public void assertMatchingAnnotationsFromGenomicsDB_newMQformat(File input, File expected, Locatable interval, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(input, interval);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expected.getAbsolutePath())));
        final List<String> attributesToFilter = Stream.concat(ATTRIBUTES_WITH_JITTER.stream(), ATTRIBUTES_TO_IGNORE.stream()).collect(Collectors.toList());
        runGenotypeGVCFSAndAssertSomething(genomicsDBUri, expected, NO_EXTRA_ARGS, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, ATTRIBUTES_WITH_JITTER, header), reference);

        // The default option with GenomicsDB input uses VCFCodec for decoding, test BCFCodec explicitly
        final List<String> args = new ArrayList<String>();
        args.add("--"+GenomicsDBArgumentCollection.USE_BCF_CODEC_LONG_NAME);
        runGenotypeGVCFSAndAssertSomething(genomicsDBUri, expected, args, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, ATTRIBUTES_WITH_JITTER, header), reference);
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testEntireVariantContext(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        assertVariantContextsMatch(input, expected, extraArgs, reference);
    }

    private void assertVariantContextsMatch(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        try {
            final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expected.getAbsolutePath())));
            runGenotypeGVCFSAndAssertSomething(input, expected, extraArgs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, ATTRIBUTES_WITH_JITTER, header), reference);
        } catch (java.io.IOException e) {
            throw new AssertionError("There was a problem reading your expected input file");
        }
    }

    private void assertGenotypesMatch(File input, File expected, List<String> additionalArguments, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, additionalArguments, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes,
                reference);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void assertDeprecatedMQThrowsUserException() {
        final File output = createTempFile("genotypegvcf", ".vcf");
        // This old gatk3 output file contains the old MQ format
        final File inputWithOldArgument = getTestFile( "combined.single.sample.pipeline.gatk3.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .add("V", inputWithOldArgument.getAbsolutePath())
                .addOutput(output);

        // This is expected to fail because RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT is not specified to allow old format MQ calculations.
        runCommandLine(args);
    }

    //this test is separate because all the others use old data and ignore the MQ annotations
    @Test(dataProvider = "GVCFsWithNewMQFormat")
    public void assertNewMQWorks(File input, File expected, Locatable interval, String reference) throws IOException {
        final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expected.getAbsolutePath())));
        runGenotypeGVCFSAndAssertSomething(input, expected, NO_EXTRA_ARGS, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, ATTRIBUTES_WITH_JITTER, header), reference);
    }

    private void runGenotypeGVCFSAndAssertSomething(File input, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input.getAbsolutePath(), expected, additionalArguments, assertion, reference
        );
    }

    private void runGenotypeGVCFSAndAssertSomething(String input, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        final File output = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected : createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .add("V", input)
                .addFlag(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false)
                .addOutput(output);

        additionalArguments.forEach(args::addRaw);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        if (! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            final List<VariantContext> expectedVC = VariantContextTestUtils.getVariantContexts(expected);
            final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output);
            try {
                assertForEachElementInLists(actualVC, expectedVC, assertion);
            } catch (final AssertionError error) {
                throw error;
            }
        }
    }

    @Test
    public void testIndexIsCreated(){
        final File output = createTempFile("test", ".vcf");
        final File index = new File(output.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX);
        Assert.assertFalse(index.exists());
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(getTestFile(BASE_PAIR_GVCF))
                .addOutput(output)
                .addReference(new File(b37_reference_20_21))
                .add(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, "true");

        runCommandLine(args);
        Assert.assertTrue(index.exists());
    }

    @Test
    public void testIntervalsAndOnlyOutputCallsStartingInIntervalsAreMutuallyRequired(){
        ArgumentsBuilder args =   new ArgumentsBuilder()
                .addVCF(getTestFile("leadingDeletion.g.vcf"))
                .addReference(new File(b37_reference_20_21))
                .addOutput( createTempFile("tmp",".vcf"))
                .add(GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME, true);

        Assert.assertThrows(CommandLineException.MissingArgument.class, () -> runCommandLine(args));
        args.add("L", "20:69512-69513");
        runCommandLine(args);
    }

    @Test
    public void testGenotypingForSomaticGVCFs() throws IOException{
        final List<Integer> starPositions = new ArrayList<>();
        starPositions.add(319);
        starPositions.add(321);
        starPositions.add(322);
        starPositions.add(323);

        final File output = createTempFile("tmp", ".vcf");
        ArgumentsBuilder args =   new ArgumentsBuilder()
                .addVCF(getTestFile("threeSamples.MT.g.vcf"))
                .addReference(new File(b37Reference))
                .addOutput(output)
                .add(CombineGVCFs.SOMATIC_INPUT_LONG_NAME, true);
        runCommandLine(args);

        //compared with the combined GVCF, this output should have called GTs and no alts with LODs less than TLOD_THRESHOLD
        //uncalled alleles should be removed

        final List<VariantContext> results = VariantContextTestUtils.getVariantContexts(output);

        //qualitative match
        for (final VariantContext vc : results) {
            Assert.assertTrue(!vc.getAlleles().contains(Allele.NON_REF_ALLELE));
            Assert.assertTrue(vc.getAlternateAlleles().size() >= 1);
            Assert.assertTrue(vc.filtersWereApplied());  //filtering should happen during combine, but make sure filters aren't dropped

            for (final Genotype g : vc.getGenotypes()) {
                Assert.assertTrue(g.isCalled());
                //homRef sites are sometimes filtered because they had low quality, filtered alleles that GGVCFs genotyped out
                //ideally if the site wasn't homRef and had a good allele it would be PASS, but htsjdk will only output genotype filters if at least one genotype is filtered
            }

            //MT:302 has an alphabet soup of alleles in the GVCF -- make sure the ones we keep are good
            if (vc.getStart() == 302) {
                Assert.assertEquals(vc.getNAlleles(), 6);
                double[] sample0LODs = VariantContextGetters.getAttributeAsDoubleArray(vc.getGenotype(0), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
                double[] sample1LODs = VariantContextGetters.getAttributeAsDoubleArray(vc.getGenotype(1), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
                double[] sample2LODs = VariantContextGetters.getAttributeAsDoubleArray(vc.getGenotype(2), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
                for (int i = 0; i < vc.getNAlleles() - 1; i++) {
                    Assert.assertTrue(sample0LODs[i] > TLOD_THRESHOLD || sample1LODs[i] > TLOD_THRESHOLD || sample2LODs[i] > TLOD_THRESHOLD);
                }
            }

            //make sure phasing is retained
            if (vc.getStart() == 317 || vc.getStart() == 320) {
                VariantContextTestUtils.assertGenotypeIsPhasedWithAttributes(vc.getGenotype(2));
            }

            //make sure *-only variants are dropped
            Assert.assertFalse(starPositions.contains(vc.getStart()));

            //sample 0 is also phased
            if (vc.getStart() == 4713 || vc.getStart() == 4720) {
                VariantContextTestUtils.assertGenotypeIsPhasedWithAttributes(vc.getGenotype(0));
            }

            //MT:4769 in combined GVCF has uncalled alleles
            if (vc.getStart() == 4769) {
                Assert.assertEquals(vc.getNAlleles(), 2);
                Assert.assertTrue(vc.getAlternateAlleles().get(0).basesMatch("G"));
            }
        }

        //exact match
        final File expectedFile = getTestFile("threeSamples.MT.vcf");
        final List<VariantContext> expected = VariantContextTestUtils.getVariantContexts(expectedFile);
        final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expectedFile.getAbsolutePath())));
        assertForEachElementInLists(results, expected,
                (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, ATTRIBUTES_WITH_JITTER, header));

    }

    @Test
    public void testGenotypingForSomaticGVCFs_subsetAlts() {
        final File output = createTempFile("tmp", ".vcf");
        ArgumentsBuilder args =   new ArgumentsBuilder()
                .addVCF(new File(getToolTestDataDir() + "../CombineGVCFs/twoSamples.MT.g.vcf"))
                .addReference(new File(b37Reference))
                .addOutput(output)
                .add(CombineGVCFs.SOMATIC_INPUT_LONG_NAME, true)
                .add("max-alternate-alleles", "2")
                .addInterval(new SimpleInterval("MT:73"));
        runCommandLine(args);

        List<VariantContext> results = VariantContextTestUtils.getVariantContexts(output);
        //MT:302 originally had 5 alts
        for (final VariantContext vc : results) {
            Assert.assertTrue(vc.getNAlleles() <= 3);  //NAlleles includes ref
        }
    }

    @Test
    public void testRawAndFinalizedAlleleSpecificAnnotationsThoroughly() {
        final File output = createTempFile("tmp", ".vcf");
        ArgumentsBuilder args =   new ArgumentsBuilder()
                .addVCF(new File(ALLELE_SPECIFIC_DIRECTORY, "NA12878.AS.chr20snippet.g.vcf"))
                .addReference(new File(b37Reference))
                .addOutput(output)
                .add("keep-combined", true)
                .add("A", "ClippingRankSumTest")
                .add("G", "AS_StandardAnnotation")
                .add("G", "StandardAnnotation")
                .add("allow-old-rms-mapping-quality-annotation-data", true);
        runCommandLine(args);

        List<VariantContext> results = VariantContextTestUtils.getVariantContexts(output);
        //there are only about 25 VCs here so we can read them all into memory
        for (final VariantContext vc : results) {
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY,
                    VCFHeaderLineCount.A, false);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY,
                    VCFHeaderLineCount.A, false);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY,
                    VCFHeaderLineCount.A, false);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_SB_TABLE_KEY,
                    VCFHeaderLineCount.R);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_FISHER_STRAND_KEY,
                    VCFHeaderLineCount.A, false);
            VariantContextTestUtils.assertAlleleSpecificAnnotationLengthsCorrect(vc, GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY,
                    VCFHeaderLineCount.A, false);
        }

    }

    @Test
    public void testForceOutput() {
        final File input = getTestFile( "combine.single.sample.pipeline.1.vcf");
        final File output1 = createTempFile("output", ".vcf");

        final ArgumentsBuilder argsWithoutForceCalling = new ArgumentsBuilder()
                .addReference(b37Reference)
                .addVCF(input)
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20")
                .addInterval(new SimpleInterval("20", 10000000, 10010000))
                .add(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT, true)
                .addOutput(output1);

        Utils.resetRandomGenerator();
        runCommandLine(argsWithoutForceCalling);

        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output1);

        // every site has output
        Assert.assertEquals(actualVC.size(), 10001);

        final File output2 = createTempFile("output", ".vcf");

        final ArgumentsBuilder argsWithForceCalling = new ArgumentsBuilder()
                .addReference(b37Reference)
                .addVCF(input)
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20:10000100")
                .addInterval(new SimpleInterval("20", 10000000, 10010000))
                .add(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT, true)
                .addOutput(output2);

        Utils.resetRandomGenerator();
        runCommandLine(argsWithForceCalling);

        final List<VariantContext> actualVC2 = VariantContextTestUtils.getVariantContexts(output2);

        // one requested site and one variant site have output
        Assert.assertEquals(actualVC2.size(), 2);
        Assert.assertEquals(actualVC2.get(0).getStart(), 10000100);
        Assert.assertTrue(actualVC2.get(0).isMonomorphicInSamples());
        Assert.assertTrue(actualVC2.get(1).isPolymorphicInSamples());
    }

    /**
     * This tests whether NON_REF alleles are properly removed, including multi-allelic sites
     */
    @Test
    public void testForceOutputNonRef() {
        final File input = new File(getToolTestDataDir() + "../CombineGVCFs/NA12878.AS.chr20snippet.g.vcf");

        // No sites should be output
        final File output1 = createTempFile("output", ".vcf");
        final ArgumentsBuilder argsWithoutForceSpecificSites = new ArgumentsBuilder()
                .addReference(b37Reference)
                .addVCF(input)
                .add(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT, true)
                .addOutput(output1);

        Utils.resetRandomGenerator();
        runCommandLine(argsWithoutForceSpecificSites);

        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output1);

        Assert.assertEquals(actualVC.size(), 24);

        // No sites should output
        final File output2 = createTempFile("output2", ".vcf");
        final ArgumentsBuilder argsWithSpecificSites = new ArgumentsBuilder()
                .addReference(b37Reference)
                .addVCF(input)
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20:10433049")
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20:10433197")
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20:10433312")
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20:10684106")
                .add(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT, true)
                .addOutput(output2);

        Utils.resetRandomGenerator();
        runCommandLine(argsWithSpecificSites);

        final List<VariantContext> actualVC2 = VariantContextTestUtils.getVariantContexts(output2);

        Assert.assertEquals(actualVC2.size(), 28);
        actualVC2.forEach(vc -> {
            Assert.assertTrue(!vc.getAlleles().contains(Allele.NON_REF_ALLELE));
        });

        for (VariantContext vc : actualVC2) {
            //If non-used alleles are pruned, this will be true
            Assert.assertEquals(!vc.isPolymorphicInSamples(), vc.getAlleles().size() == 1);
            Assert.assertEquals(vc.getAlleles(), vc.subContextFromSamples(vc.getSampleNames(), true).getAlleles());

            if (vc.getStart() == 10433049) {
                Assert.assertEquals(vc.getAlleles(), Arrays.asList(Allele.REF_C));
            }
            else if (vc.getStart() == 10433197) {
                Assert.assertEquals(vc.getAlleles(), Arrays.asList(Allele.REF_C));
            }
            else if (vc.getStart() == 10433312) {
                Assert.assertEquals(vc.getAlleles(), Arrays.asList(Allele.create("CAAAAAAA", true)));
            }
            else if (vc.getStart() == 10684106) {
                Assert.assertEquals(vc.getAlleles(), Arrays.asList(Allele.create("CCTTTCTTTCTTT", true)));
            }
        }
    }

    @Test
    public void testForceOutputWithSpanningDeletion() {
        final File input = getTestFile("leadingDeletion.g.vcf");
        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(b37Reference)
                .addVCF(input)
                .add(GenotypeGVCFs.FORCE_OUTPUT_INTERVALS_NAME, "20")
                .addInterval(new SimpleInterval("20", 69511, 69515))
                .addOutput(output);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output);

        Assert.assertEquals(actualVC.size(), 5);
        for (final int n : new int[] {1, 2, 3}) {
            Assert.assertTrue(actualVC.get(n).getAlternateAlleles().stream().anyMatch(a -> a == Allele.SPAN_DEL));
        }
    }

    @Test
    public void testWithReblockedGVCF() {
        final File reblockedGVCF = new File("src/test/resources/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFs/twoReblocked.g.vcf");
        final File output = createTempFile("reblockedAndGenotyped", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(b37_reference_20_21)
                .addVCF(reblockedGVCF)
                .addOutput(output);
        /*runCommandLine(args);

        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output);
        Assert.assertTrue(actualVC.stream().anyMatch(vc -> vc.getGenotype(1).isHomRef() && !vc.getGenotype(1).hasPL()));  //second sample has a bunch of 0/0s
*/
        final File bigCombinedReblockedGVCF = new File("src/test/resources/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFs/combineReblocked.g.vcf");
        final File cohortOutput = createTempFile("biggerCohort.rb", ".vcf");

        final ArgumentsBuilder args2 = new ArgumentsBuilder();
        args2.addReference(hg38Reference)
                .addVCF(bigCombinedReblockedGVCF)
                .addOutput(cohortOutput);
        runCommandLine(args2);

        final List<VariantContext> outputVCs = VariantContextTestUtils.getVariantContexts(cohortOutput);
        final VariantContext vc0 = outputVCs.get(0);
        Assert.assertTrue(vc0.getAttributeAsDouble(GATKVCFConstants.EXCESS_HET_KEY, 1000.0) < 10.0);
        Assert.assertTrue(vc0.hasAttribute(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY));  //will get dropped if homRefs aren't counted
        Assert.assertEquals(vc0.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, 0), 362);
        Assert.assertEquals(vc0.getAlternateAlleles().size(), 1);  //had another low quality alt
        Assert.assertEquals(vc0.getAlternateAllele(0).getBaseString(), "G");
        Assert.assertTrue(vc0.getGenotypes().stream().allMatch(g -> g.isCalled() && g.hasGQ() && g.hasDP()));

        final VariantContext vc1 = outputVCs.get(1);
        Assert.assertTrue(vc1.getAttributeAsDouble(GATKVCFConstants.EXCESS_HET_KEY, 1000.0) < 10.0); //will be ~72 if homRefs aren't counted
        Assert.assertTrue(vc1.hasAttribute(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY));
        Assert.assertEquals(vc1.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, 0), 362);
        Assert.assertEquals(vc1.getAlternateAlleles().size(), 3);
        Assert.assertTrue(vc1.isIndel());
        Assert.assertTrue(vc0.getGenotypes().stream().allMatch(g -> g.isCalled() && g.hasGQ() && g.hasDP()));
    }
}
