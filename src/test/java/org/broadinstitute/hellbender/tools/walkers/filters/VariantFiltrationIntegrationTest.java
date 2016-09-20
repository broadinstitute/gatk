package org.broadinstitute.hellbender.tools.walkers.filters;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class VariantFiltrationIntegrationTest extends CommandLineProgramTest {

    public String baseTestString(final String vcf, final String options) {
        final String file = getToolTestDataDir() + vcf;
        return "--variant " + file + " " + options + " -O %s" + " -R" + hg19_chr1_1M_Reference;
    }

    @Test
    public void testNoAction() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", ""),//"-L 1:10,020,000-10,021,000"
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testNoAction.vcf")
        );

        spec.executeTest("testNoAction", this);
    }

    @Test
    public void testClusteredSnps() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
        baseTestString("vcfexample2.vcf", " -window 10 "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testClusteredSnps.vcf")
        );

        spec.executeTest("testClusteredSnps", this);
    }

    @DataProvider(name="masks")
    public Object[][] masks() {
        return new String[][]{
                {"foo", "--mask " + getToolTestDataDir() + "vcfexample2.vcf", "testVariantFiltration_testMask1.vcf"},
                {"foo", "--mask VCF:" + getToolTestDataDir() + "vcfMask.vcf", "testVariantFiltration_testMask2.vcf"},
                {"foo", "-maskExtend 10 --mask VCF:" + getToolTestDataDir() + "vcfMask.vcf", "testVariantFiltration_testMask3.vcf"},
        };
    }

    @Test(dataProvider = "masks")
    public void testMask(final String maskName, final String mask, final String expected) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -maskName " + maskName + " " + mask),
                Arrays.asList(getToolTestDataDir() + "expected/" + expected)
        );

        spec.executeTest("testMask", this);
    }

    @Test
    public void testMaskReversed() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -maskName outsideGoodSites -filterNotInMask --mask BED:" + getToolTestDataDir() + "goodMask.bed"),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testMaskReversed.vcf")
        );

        spec.executeTest("testMaskReversed", this);
    }

    @Test
    public void testIllegalFilterName() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName 'foo < foo' "),
                1,
                UserException.class
        );

        spec.executeTest("testIllegalFilterName", this);
    }

    @Test
    public void testFilter1() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testFilter1.vcf")
        );

        spec.executeTest("testFilter1", this);
    }

    @Test
    public void testFilter2() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testFilter2.vcf")
        );

        spec.executeTest("testFilter2", this);
    }

    @Test
    public void testFilterWithSeparateNames() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testFilterWithSeparateNames.vcf")
        );

        spec.executeTest("testFilterWithSeparateNames", this);
    }

    @Test
    public void testInvertFilter() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --invertFilterExpression "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testInvertFilter.vcf")
        );

        spec.executeTest("testInvertFilter", this);
    }

    @Test
    public void testInvertJexlFilter() throws IOException {
        //Note: the "invert" in the name refers to the logic being the opposite of testFilterWithSeparateNames (and same as testInvertFilter)
        //Note: Output differs from testInvertFilter because FILTER description uses the -genotypeFilterExpression argument
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " --filterName ABF -filter 'AlleleBalance >= 0.7' --filterName FSF -filter 'FisherStrand != 1.4' "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testInvertJexlFilter.vcf")
        );

        spec.executeTest("testInvertJexlFilter", this);
    }

    @Test
    public void testGenotypeFilters1() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -G_filter 'GQ == 0.60' -G_filterName foo "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testGenotypeFilters1.vcf")
        );

        spec.executeTest("testGenotypeFilters1", this);
    }

    @Test
    public void testGenotypeFilters2() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("vcfexample2.vcf", " -G_filter 'isHomVar == 1' -G_filterName foo "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testGenotypeFilters2.vcf")
        );

        spec.executeTest("testGenotypeFilters2", this);
    }

    @Test
    public void testDeletions() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("twoDeletions.vcf", " --filterExpression 'QUAL < 100' --filterName foo "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testDeletions.vcf")
        );

        spec.executeTest("testDeletions", this);
    }

    @Test
    public void testUnfilteredBecomesFilteredAndPass() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("unfilteredForFiltering.vcf", " --filterExpression 'FS > 60.0' --filterName SNP_FS "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testUnfilteredBecomesFilteredAndPass.vcf")
        );

        spec.executeTest("testUnfilteredBecomesFilteredAndPass", this);
    }

    @Test
    public void testFilteringDPfromINFO() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("filteringDepthInFormat.vcf", " --filterExpression 'DP < 8' --filterName lowDP "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testFilteringDPfromINFO.vcf")
        );

        spec.executeTest("testFilteringDPfromINFO", this);
    }

    @Test
    public void testFilteringDPfromFORMAT() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("filteringDepthInFormat.vcf", " -genotypeFilterExpression 'DP < 8' --genotypeFilterName lowDP "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testFilteringDPfromFORMAT.vcf")
        );

        spec.executeTest("testFilteringDPfromFORMAT", this);
    }

    @Test
    public void testInvertGenotypeFilterExpression() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("filteringDepthInFormat.vcf", " --genotypeFilterExpression 'DP < 8' --genotypeFilterName highDP --invertGenotypeFilterExpression "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testInvertGenotypeFilterExpression.vcf")
        );

        spec.executeTest("testInvertGenotypeFilterExpression", this);
    }

    @Test
    public void testInvertJexlGenotypeFilterExpression() throws IOException {
        //Note: the "invert" in the name refers to the logic being the opposite of testFilteringDPfromFORMAT (and same as testInvertGenotypeFilterExpression_
        //Note: Output differs from testInvertGenotypeFilterExpression because FILTER description uses the -genotypeFilterExpression argument
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("filteringDepthInFormat.vcf", " --genotypeFilterExpression 'DP >= 8' --genotypeFilterName highDP "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testInvertJexlGenotypeFilterExpression.vcf")
        );

        spec.executeTest("testInvertJexlGenotypeFilterExpression", this);
    }

    @Test
    public void testSetFilteredGtoNocall() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("filteringDepthInFormat.vcf", " --genotypeFilterExpression 'DP < 8' --genotypeFilterName lowDP --setFilteredGtToNocall "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testSetFilteredGtoNocall.vcf")
        );

        spec.executeTest("testSetFilteredGtoNocall", this);
    }

    @Test
    public void testSetFilteredGtoNocallUpdateInfo()  throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("variantFiltrationInfoField.vcf", " -G_filter 'GQ < 20' -G_filterName lowDP -G_filter 'DP < 10' -G_filterName lowGQ --setFilteredGtToNocall "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testSetFilteredGtoNocallUpdateInfo.vcf")
        );

        spec.executeTest("testSetFilteredGtoNocallUpdateInfo", this);
    }

    @Test
    public void testSetVcfFilteredGtoNocall()  throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString("filteredSamples.vcf", " --setFilteredGtToNocall "),
                Arrays.asList(getToolTestDataDir() + "expected/" + "testVariantFiltration_testSetVcfFilteredGtoNocall.vcf")
        );

        spec.executeTest("testSetVcfFilteredGtoNocall", this);
    }
}
