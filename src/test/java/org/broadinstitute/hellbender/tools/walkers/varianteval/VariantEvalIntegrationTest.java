package org.broadinstitute.hellbender.tools.walkers.varianteval;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class VariantEvalIntegrationTest extends CommandLineProgramTest {

    private static String fundamentalTestVCF = "FundamentalsTest.annotated.db.subset.snps_and_indels.vcf";
    private static String fundamentalTestSNPsVCF = "FundamentalsTest.annotated.db.subset.final.vcf";
    private static String fundamentalTestSNPsWithMLEVCF = "FundamentalsTest.annotated.db.subset.final.withMLE.vcf";
    private static String fundamentalTestSNPsSplit1of2VCF = "FundamentalsTest.annotated.db.subset.final.split_1_of_2.vcf";
    private static String fundamentalTestSNPsSplit2of2VCF = "FundamentalsTest.annotated.db.subset.final.split_2_of_2.vcf";
    private static String fundamentalTestSNPsOneSampleVCF = "FundamentalsTest.annotated.db.subset.final.NA12045.vcf";
    private static String snpEffVcf = largeFileTestDir + "snpEff2.0.5.AFR.unfiltered.VariantAnnotator.output.vcf";

    private String getExpectedFile(String testName) {
        return getToolTestDataDir() + "expected/" + testName + ".expected.txt";
    }
    
    protected String getTestFilePath(String fileName) {
        return super.getTestFile(fileName).getPath();
    }

    @Test
    public void testFunctionClassWithSnpeff() throws IOException {
        String name = "testFunctionClassWithSnpeff";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                                        " -R " + b37Reference +
                                        " --dbsnp " + dbsnp_138_b37_1_65M_vcf +
                                        " --eval " + snpEffVcf +
                                        " -no-ev" +
                                        " -EV TiTvVariantEvaluator" +
                                        " -no-st" +
                                        " -ST FunctionalClass" +
                                        " -L " + snpEffVcf +
                                        " -O %s"
                                , Arrays.asList(getExpectedFile(name)));

            spec.executeTest(name, this);
    }

    @Test
    public void testStratifySamplesAndExcludeMonomorphicSites() throws IOException {
        String vcf = getTestFilePath("/CEU.trio.callsForVE.vcf");
        String name = "testStratifySamplesAndExcludeMonomorphicSites";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                                    " -R " + b37Reference +
                                    " --dbsnp " + dbsnp_138_b37_1_65M_vcf+
                                    " --eval " + vcf +
                                    " -no-ev" +
                                    " -EV TiTvVariantEvaluator" +
                                    " -ST Sample" +
                                    " -L " + vcf +
                                    " -O %s"
                            , Arrays.asList(getExpectedFile(name))
                              );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndels() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndels";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                                        " -R " + b37_reference_20_21 +
                                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                                        " -no-ev" +
                                        " -EV CountVariants" +
                                        " -no-st" +
                                        " -L " + getTestFilePath(fundamentalTestVCF) +
                                        " -O %s", Arrays.asList(getExpectedFile(name))
                              );
            spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNovelty() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithNovelty";
        getTestFilePath(fundamentalTestVCF);
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST Novelty" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @DataProvider(name = "testContigStratWithUserSuppliedIntervalsData")
    public Object[][] testContigStratWithUserSuppliedIntervalsData() {
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{"1:1-1480226", "testContigStratWithUserSuppliedIntervals"});
        tests.add(new Object[]{"1", "testContigStratWithUserSuppliedIntervals2"});
        tests.add(new Object[]{null, "testContigStratWithUserSuppliedIntervals3"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "testContigStratWithUserSuppliedIntervalsData")
    public void testContigStratWithUserSuppliedIntervals(String intervalString, String expectedOutputFile) throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -R " + b37Reference +
                        " --eval " + snpEffVcf +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST Contig" +
                        (intervalString == null ? "" : " -L " + intervalString) +
                        " -O %s",
                Arrays.asList(getExpectedFile(expectedOutputFile))
        );
        spec.executeTest(expectedOutputFile, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithNoveltyAndFilter";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST Novelty" +
                        " -ST Filter" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithCpG() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithCpG";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST CpG" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClasses() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithFunctionalClass";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST FunctionalClass" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithDegeneracy";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST Degeneracy" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithSample() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithSample";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST Sample" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithJexlExpression";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST JexlExpression" +
                        " -select 'DP < 20'" +
                        " -select-name DepthSelect" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions() throws IOException {
        String name = "testFundamentalsCountVariantsSNPsAndIndelsWithMultipleJexlExpressions";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST JexlExpression" +
                        " -select 'DP < 20'" +
                        " -select-name DepthLt20" +
                        " -select 'DP > 20'" +
                        " -select-name DepthGt20" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testFundamentalsCountVariantsNoCompRod() throws IOException {
        String name = "testFundamentalsCountVariantsNoCompRod";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --eval " + getTestFilePath(fundamentalTestVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -L " + getTestFilePath(fundamentalTestVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testSelect1() throws IOException {
        String name = "testSelect1";
        String vcf1 = getTestFilePath("yri.trio.gatk_glftrio.intersection.annotated.filtered.chr1.500.noheader.vcf");
        String vcf2 = getTestFilePath("yri.trio.gatk.ug.head.vcf");

        String tests = " -R " + b37Reference +
                " --dbsnp " + dbsnp_138_b37_1_65M_vcf +
                " --eval " + vcf1 +
                " --comp:comp_genotypes " + vcf2;
        IntegrationTestSpec spec = new IntegrationTestSpec(withSelect(tests, "DP < 50", "DP50") + " " + " -ST CpG -O %s", Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name, this);
    }

    @Test
    public void testVEMendelianViolationEvaluator() throws IOException {
        String vcfFile = getTestFilePath("MendelianViolationEval.vcf");
        String pedFile = "MendelianViolationEval.ped";

        String name = "testVEMendelianViolationEvaluator";
        IntegrationTestSpec spec = new IntegrationTestSpec(" -R " + b37Reference + " --eval " + vcfFile + " -ped "+ getTestFilePath(pedFile) +" -no-ev -EV MendelianViolationEvaluator -L 1:10109-10315 -O %s -mvq 0 -no-st",
                Arrays.asList(getExpectedFile(name)));

        spec.executeTest(name, this);
    }

    @Test
    public void testMVEvalFamilyStrat() throws IOException {
        String name = "testMVEvalFamilyStrat";
        String vcfFile = "PhaseByTransmission.IntegrationTest.TP.vcf";
        String pedFile = "PhaseByTransmission.IntegrationTest.goodFamilies.ped";

        IntegrationTestSpec spec = new IntegrationTestSpec(" -R " + b37Reference + " -ped " + getTestFilePath(pedFile) + " -eval " + getTestFilePath(vcfFile) + " -no-ev -no-st -ST Family -EV MendelianViolationEvaluator -O %s",
                Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name, this);
    }


    private static String withSelect(String cmd, String select, String name) {
        return String.format("%s -select '%s' -select-name %s", cmd, select, name);
    }

    @Test
    public void testCompOverlap() throws IOException {
        String eval = getTestFilePath("pacbio.ts.recalibrated.vcf");
        String comp = largeFileTestDir + "genotypes_r27_nr.b37_fwd.subset.vcf";

        String name = "testCompOverlap";
        String extraArgs = " -R " + b37_reference_20_21 + " -L " + getTestFilePath("pacbio.hg19.intervals") + " --comp:comphapmap " + comp + " --eval " + eval + " -no-ev -EV CompOverlap -sn NA12878 -no-st -ST Novelty -O %s";
        IntegrationTestSpec spec = new IntegrationTestSpec(extraArgs,Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name,this);
    }

    @Test
    public void testEvalTrackWithoutGenotypes() throws IOException {
        String vcf = largeFileTestDir + "ALL.20100201.chr20.subset.bi.sites.vcf";
        String name = "testEvalTrackWithoutGenotypes";
        String extraArgs = " -R " + b37_reference_20_21 +
                           " -L 20" +
                           " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                           " --eval:evalBI " + vcf +
                           " -no-st -ST Novelty -O %s";
        IntegrationTestSpec spec = new IntegrationTestSpec(extraArgs,Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name,this);
    }

    @Test
    public void testEvalTrackWithoutGenotypesWithSampleFields() throws IOException {
        String name = "testEvalTrackWithoutGenotypesWithSampleFields";
        String vcf = getTestFilePath("noGenotypes.vcf");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37Reference +
                        " -eval " + vcf +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))); //There is no md5 because we only care that this completes without an exception.
        spec.executeTest(name, this);

    }

    @Test
    public void testEvalTrackWithoutGenotypesWithSampleFieldsWrongRef() throws IOException {
        String name = "testEvalTrackWithoutGenotypesWithSampleFieldsWrongRef";
        String vcf = getTestFilePath("noGenotypes.vcf");

        //The VCF has chr 1, but FASTA has 20/21 only
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " -eval " + vcf +
                        " -O %s",
                1, UserException.class);
        spec.executeTest(name, this);
    }

    @Test
    public void testMultipleEvalTracksWithoutGenotypes() throws IOException {
        String name = "testMultipleEvalTracksWithoutGenotypes";
        String vcf1 = largeFileTestDir + "ALL.20100201.chr20.subset.bi.sites.vcf";
        String vcf2 = largeFileTestDir + "ALL.20100201.chr20.subset.bc.sites.vcf";

        String extraArgs =
                " -R " + b37_reference_20_21 +
                " -L 20" +
                " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                " --eval:evalBI " + vcf1 +
                " --eval:evalBC " + vcf2 +
                " -no-st -ST Novelty -O %s";
        IntegrationTestSpec spec = new IntegrationTestSpec(extraArgs,Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name,this);
    }

    @Test
    public void testMultipleCompTracks() throws IOException {
        String vcf = largeFileTestDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf";
        String eval = getTestFilePath("NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.subset.vcf");

        String extraArgs = " -R " + b37_reference_20_21 +
                           " --comp " + vcf +
                           " --eval " + eval +
                           " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                           " -L 20:10000000-10100000" +
                           " -no-st -no-ev -ST Novelty -EV CompOverlap" +
                           " -O %s";

        String name = "testMultipleCompTracks";
        IntegrationTestSpec spec = new IntegrationTestSpec(extraArgs,Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name,this);
    }

    @Test
    public void testPerSampleAndSubsettedSampleHaveSameResults() throws IOException {
        String name = "testPerSampleAndSubsettedSampleHaveSameResults-subset";
        IntegrationTestSpec spec = new IntegrationTestSpec(
              " -R " + b37_reference_20_21 +
                    " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                    " --eval " + getTestFilePath(fundamentalTestSNPsVCF) +
                    " -no-ev" +
                    " -EV CompOverlap" +
                    " -sn NA12045" +
                    " -no-st" +
                    " -L " + getTestFilePath(fundamentalTestSNPsVCF) +
                    " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);

        String name2 = "testPerSampleAndSubsettedSampleHaveSameResults-Onesample";
        IntegrationTestSpec spec2 = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestSNPsOneSampleVCF) +
                        " -no-ev" +
                        " -EV CompOverlap" +
                        " -no-st" +
                        " -L " + getTestFilePath(fundamentalTestSNPsOneSampleVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name2))
        );
        spec2.executeTest(name2, this);
    }


    @Test
    public void testAlleleCountStrat() throws IOException {
        String name = "testAlleleCountStrat";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestSNPsVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST AlleleCount" +
                        " -L " + getTestFilePath(fundamentalTestSNPsVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testAlleleCountStratWithMLE() throws IOException {
        String name = "testAlleleCountStratWithMLE";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestSNPsWithMLEVCF) +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST AlleleCount" +
                        " -L " + getTestFilePath(fundamentalTestSNPsWithMLEVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testMultipleEvalTracksAlleleCountWithMerge() throws IOException {
        String name = "testMultipleEvalTracksAlleleCountWithMerge";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestSNPsSplit1of2VCF) +
                        " --eval " + getTestFilePath(fundamentalTestSNPsSplit2of2VCF) +
                        " --merge-evals" +
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST AlleleCount" +
                        " -L " + getTestFilePath(fundamentalTestSNPsVCF) +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testMultipleEvalTracksAlleleCountWithoutMerge() throws IOException {
        String name = "testMultipleEvalTracksAlleleCountWithoutMerge";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " --dbsnp " + dbsnp_138_b37_20_21_vcf +
                        " --eval " + getTestFilePath(fundamentalTestSNPsSplit1of2VCF) +
                        " --eval " + getTestFilePath(fundamentalTestSNPsSplit2of2VCF) +
                        //"--merge-evals" + No merge with AC strat ==> error
                        " -no-ev" +
                        " -EV CountVariants" +
                        " -no-st" +
                        " -ST AlleleCount" +
                        " -L " + getTestFilePath(fundamentalTestSNPsVCF) +
                        " -O %s",
                1,
                CommandLineException.class
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testIntervalStrat() throws IOException {
        String name = "testIntervalStrat";
        String vcf = getTestFilePath("/withSymbolic.b37.vcf");
        String bed = getTestFilePath("/overlapTest.bed");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                                  " -R " + b37_reference_20_21 +
                                        " -eval " + vcf +
                                        " -no-ev" +
                                        " -EV CountVariants" +
                                        " -no-st" +
                                        " -strat-intervals " + bed +
                                        " -ST IntervalStratification" +
                                        " -L 20" +
                                        " -O %s",
                                Arrays.asList(getExpectedFile(name))
                              );
        spec.executeTest(name, this);
    }

    @Test
    public void testModernVCFWithLargeIndels() throws IOException {
        String name = "testModernVCFWithLargeIndels";
        String vcf = largeFileTestDir + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.chr20.vcf";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                                  " -R " + b37_reference_20_21 +
                                        " -eval " + vcf +
                                        " -L 20" +
                                        " -D " + dbsnp_138_b37_20_21_vcf +
                                        " -O %s",
                                Arrays.asList(getExpectedFile(name))
                              );
        spec.executeTest(name, this);
    }

    @Test
    public void testStandardIndelEval() throws IOException {
        String name = "testStandardIndelEval";
        String vcf = largeFileTestDir + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.chr20.vcf";
        String gold = largeFileTestDir + "/Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " -eval " + vcf +
                        " -L 20" +
                        " -no-st -ST Sample -ST OneBPIndel -ST TandemRepeat" +
                        " -no-ev -EV IndelSummary -EV IndelLengthHistogram" +
                        " -gold " + gold +
                        " -D " + dbsnp_138_b37_20_21_vcf +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))
        );
        spec.executeTest(name, this);
    }

    @Test
    public void testBadACValue() throws IOException {
        String name = "testBadACValue";
        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " -eval " + getTestFilePath("vcfexample.withBadAC.vcf") +
                        " -no-st -ST AlleleCount" +
                        " -no-ev -EV VariantSummary" +
                        " -O %s",
                1,
                UserException.class);
        spec.executeTest(name, this);
    }


    @Test()
    public void testIncompatibleEvalAndStrat() throws IOException {
        String name = "testIncompatibleEvalAndStrat";
        String vcf = largeFileTestDir + "/NA12878.HiSeq.WGS.b37_decoy.indel.recalibrated.chr20.vcf";

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " -eval " + vcf +
                        " -L 20 -no-st -ST AlleleCount -no-ev -EV VariantSummary -O %s",
                1,
                CommandLineException.BadArgumentValue.class);
        spec.executeTest(name, this);
    }

    public void testIncludingAC0(boolean includeAC0) throws IOException {
        String name = "testIncludingAC0 keep ac 0 = " + includeAC0;
        String vcf = getTestFilePath("ac0.vcf");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " -eval " + vcf +
                        " -L 20:81006 -no-st -no-ev -EV VariantSummary -O %s" + (includeAC0 ? " -keep-ac0" : ""),
                Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name, this);
    }

    @Test public void testWithAC0() throws IOException { testIncludingAC0(true); }
    @Test public void testWithoutAC0() throws IOException { testIncludingAC0(false); }

    //
    // Test validation report is doing the right thing with sites only and genotypes files
    // where the validation comp has more genotypes than eval
    //
    @Test(dataProvider = "testValidationReportData")
    public void testValidationReport(final String suffix, final String eval, final String comp) throws IOException {
        String name = "testValidationReportData-" + suffix;

        IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + b37_reference_20_21 +
                        " -eval " + eval +
                        " -comp " + comp +
                        " -L 20:10,000,000-10,000,010 -no-st -no-ev -EV ValidationReport -O %s",
                Arrays.asList(getExpectedFile(name)));
        spec.executeTest(name, this);
    }

    @DataProvider(name = "testValidationReportData")
    public Object[][] testValidationReportData() {
        final String compGenotypes = getTestFilePath("validationReportComp.vcf");
        final String compSites = getTestFilePath("validationReportComp.noGenotypes.vcf");
        final String evalGenotypes = getTestFilePath("validationReportEval.vcf");
        final String evalSites = getTestFilePath("validationReportEval.noGenotypes.vcf");

        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{"sites-sites", evalSites, compSites});
        tests.add(new Object[]{"sites-genotypes", evalSites, compGenotypes});
        tests.add(new Object[]{"genotypes-sites", evalGenotypes, compSites});
        tests.add(new Object[]{"genotypes-genotypes", evalGenotypes, compGenotypes});
        return tests.toArray(new Object[][]{});
    }

    @Test
    public void testPrintMissingComp() throws IOException {
        String name = "testPrintMissingComp";
        String vcf1 = getTestFilePath("validationReportEval.noGenotypes.vcf");
        String vcf2 = getTestFilePath("validationReportComp.noGenotypes.vcf");

        IntegrationTestSpec spec = new IntegrationTestSpec(
                        " -R " + b37_reference_20_21 +
                        " -eval " + vcf1 +
                        " --comp " + vcf2 +
                        " -L 20" +
                        " -EV PrintMissingComp" +
                        " -O %s",
                Arrays.asList(getExpectedFile(name))); // sato: make sure it doesn't throw a null pointer exception.
        spec.executeTest(name, this);

    }

    @Test
    public void testOutputFileCreation() throws IOException {
        // duplicate of "testPrintMissingComp", this time without using IntegrationTestSpec in order to force the
        // tool to create the output file directly (tests fix for https://github.com/broadinstitute/gatk/issues/5674)

        final String testName = "testOutputFileCreation";
        final File tmpDir = createTempDir(testName);
        final File outputFile = new File(tmpDir, testName + ".txt");

        ArgumentsBuilder argBuilder= new ArgumentsBuilder();
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.add("eval", getTestFilePath("validationReportEval.noGenotypes.vcf"));
        argBuilder.add(StandardArgumentDefinitions.COMPARISON_SHORT_NAME, getTestFilePath("validationReportComp.noGenotypes.vcf"));
        argBuilder.addInterval(new SimpleInterval("20"));
        argBuilder.add("EV", "PrintMissingComp");
        argBuilder.addOutput(outputFile);

        runCommandLine(argBuilder);

        IntegrationTestSpec.assertEqualTextFiles(
                outputFile,
                new File(getToolTestDataDir() + "expected/" + "testPrintMissingComp" + ".expected.txt"));
    }

    private void testForCrashWithGivenEvaluator(final String countVariants) {
        final String vcf = getTestFilePath("/CEU.trio.callsForVE.vcf");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput( createTempFile("out",".stuff"))
                .add("eval", vcf)
                .addFlag("do-not-use-all-standard-modules")
                .add("EV", countVariants);

        runCommandLine(args);
    }
}
