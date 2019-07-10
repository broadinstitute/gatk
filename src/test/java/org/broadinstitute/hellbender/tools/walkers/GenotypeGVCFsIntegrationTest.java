package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.VCFHeaderReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.codec.digest.DigestUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.GenomicsDBTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
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
    private static final List<String> ATTRIBUTES_TO_IGNORE = Arrays.asList("AS_QD","QD","FS","RAW_MQ","RGQ","MQ"); //MQ data format and key have changed since GATK3

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

    @DataProvider(name = "gvcfWithPPs")
    public Object[][] gvcfWithPPs() {
        return new Object[][] {
                {getTestFile("../../GenomicsDBImport/expected.testGVCFMode.gatk4.posteriors.g.vcf"),
                        getTestFile( "expected.posteriors.genotyped.vcf"), NO_EXTRA_ARGS, b37_reference_20_21}
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
        assertGenotypesMatch(input, expected, extraArgs, reference);
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
    public void assertMatchingGenotypesFromTileDB(File input, File expected, Locatable interval, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(input, interval);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);
        runGenotypeGVCFSAndAssertSomething(genomicsDBUri, expected, NO_EXTRA_ARGS, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes, reference);
    }

    @Test(dataProvider = "getGVCFsForGenomicsDBOverMultipleIntervals")
    public void testGenotypeGVCFsMultiIntervalGDBQuery(File input, File expected, List<Locatable> intervals, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(input, intervals, true);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .addArgument("V", genomicsDBUri)
                .addOutput(output);
        intervals.forEach(args::addInterval);
        args.add("--" + GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME);
        args.add("--only-output-calls-starting-in-intervals");  //note that this will restrict calls to just the specified intervals

        Utils.resetRandomGenerator();
        runCommandLine(args);

        //Note that if this isn't working it will take *FOREVER*
        // runs in 0.06 minutes with no input intervals specfied
        final List<VariantContext> expectedVC = VariantContextTestUtils.getVariantContexts(expected);
        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes);

    }

    //this tests single-sample with new MQ format
    @Test (dataProvider = "singleSampleGVCFWithNewMQFormat")
    public void assertMatchingAnnotationsFromGenomicsDB_newMQformat(File input, File expected, Locatable interval, String reference) throws IOException {
        final File tempGenomicsDB = GenomicsDBTestUtils.createTempGenomicsDB(input, interval);
        final String genomicsDBUri = GenomicsDBTestUtils.makeGenomicsDBUri(tempGenomicsDB);

        final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expected.getAbsolutePath())));
        final List<String> attributesToIgnore = Stream.concat(ATTRIBUTES_WITH_JITTER.stream(), ATTRIBUTES_TO_IGNORE.stream()).collect(Collectors.toList());
        runGenotypeGVCFSAndAssertSomething(genomicsDBUri, expected, NO_EXTRA_ARGS, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, attributesToIgnore, header), reference);
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testEntireVariantContext(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        assertVariantContextsMatch(input, expected, extraArgs, reference);
    }

    private void assertVariantContextsMatch(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        try {
            final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expected.getAbsolutePath())));
            runGenotypeGVCFSAndAssertSomething(input, expected, extraArgs, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, header), reference);
        } catch (java.io.IOException e) {
            throw new AssertionError("There was a problem reading your expected input file");
        }
    }

    private void assertGenotypesMatch(File input, File expected, List<String> additionalArguments, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, additionalArguments, VariantContextTestUtils::assertVariantContextsHaveSameGenotypes,
                reference);
    }

    @Test(dataProvider = "gvcfWithPPs")
    public void assertPPsAreStripped(File input, File expected, List<String> extraArgs, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input, expected, extraArgs, VariantContextTestUtils::assertGenotypePosteriorsAttributeWasRemoved, reference);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void assertDeprecatedMQThrowsUserException() {
        final File output = createTempFile("genotypegvcf", ".vcf");
        // This old gatk3 output file contains the old MQ format
        final File inputWithOldArgument = getTestFile( "combined.single.sample.pipeline.gatk3.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addArgument("V", inputWithOldArgument.getAbsolutePath())
                .addOutput(output);

        // This is expected to fail because RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT is not specified to allow old format MQ calculations.
        runCommandLine(args);
    }

    //this test is separate because all the others use old data and ignore the MQ annotations
    @Test(dataProvider = "GVCFsWithNewMQFormat")
    public void assertNewMQWorks(File input, File expected, Locatable interval, String reference) throws IOException {
        final VCFHeader header = VCFHeaderReader.readHeaderFrom(new SeekablePathStream(IOUtils.getPath(expected.getAbsolutePath())));
        runGenotypeGVCFSAndAssertSomething(input, expected, NO_EXTRA_ARGS, (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_WITH_JITTER, header), reference);
    }

    private void runGenotypeGVCFSAndAssertSomething(File input, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        runGenotypeGVCFSAndAssertSomething(input.getAbsolutePath(), expected, additionalArguments, assertion, reference
        );
    }

    private void runGenotypeGVCFSAndAssertSomething(String input, File expected, List<String> additionalArguments, BiConsumer<VariantContext, VariantContext> assertion, String reference) throws IOException {
        final File output = createTempFile("genotypegvcf", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(reference))
                .addArgument("V", input)
                .addArgument(RMSMappingQuality.RMS_MAPPING_QUALITY_OLD_BEHAVIOR_OVERRIDE_ARGUMENT)
                .addOutput(output);

        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = VariantContextTestUtils.getVariantContexts(expected);
        final List<VariantContext> actualVC = VariantContextTestUtils.getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, assertion);
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
                .addArgument(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, "true");

        runCommandLine(args);
        Assert.assertTrue(index.exists());
    }

    @Test
    public void testIntervalsAndOnlyOutputCallsStartingInIntervalsAreMutuallyRequired(){
        ArgumentsBuilder args =   new ArgumentsBuilder()
                .addVCF(getTestFile("leadingDeletion.g.vcf"))
                .addReference(new File(b37_reference_20_21))
                .addOutput( createTempFile("tmp",".vcf"))
                .addBooleanArgument(GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME, true);

        Assert.assertThrows(CommandLineException.MissingArgument.class, () -> runCommandLine(args));
        args.addArgument("L", "20:69512-69513");
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
                .addBooleanArgument(CombineGVCFs.SOMATIC_INPUT_LONG_NAME, true);
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
                double[] sample0LODs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc.getGenotype(0), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
                double[] sample1LODs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc.getGenotype(1), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
                double[] sample2LODs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc.getGenotype(2), GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, () -> null, 0.0);
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
                (a, e) -> VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, ATTRIBUTES_TO_IGNORE, header));

    }

    @Test
    public void testGenotypingForSomaticGVCFs_subsetAlts() {
        final File output = createTempFile("tmp", ".vcf");
        ArgumentsBuilder args =   new ArgumentsBuilder()
                .addVCF(new File(getToolTestDataDir() + "../CombineGVCFs/twoSamples.MT.g.vcf"))
                .addReference(new File(b37Reference))
                .addOutput(output)
                .addBooleanArgument(CombineGVCFs.SOMATIC_INPUT_LONG_NAME, true)
                .addArgument("max-alternate-alleles", "2")
                .addInterval(new SimpleInterval("MT:73"));
        runCommandLine(args);

        List<VariantContext> results = VariantContextTestUtils.getVariantContexts(output);
        //MT:302 originally had 5 alts
        for (final VariantContext vc : results) {
            Assert.assertTrue(vc.getNAlleles() <= 3);  //NAlleles includes ref
        }
    }
}
