package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.TreeMap;

import static org.broadinstitute.hellbender.tools.variantdb.nextgen.FilterSensitivityTools.*;
import static org.testng.Assert.*;

public class FilterSensitivityToolsTest {

    // for testing inputs
    private Double definedDoubleInput = 0.0;
    private Double undefinedDoubleInput = null;
    private String definedStringInput = "I'm defined!";
    private String undefinedStringInput = null;

    private Map<Double, Double> testTrancheMap = new TreeMap<>();

    @BeforeMethod
    public void setUp() {
        // deliberately construct out of order to test ordering of TreeMap
        testTrancheMap.put(95.0, 2.5);
        testTrancheMap.put(99.0, -6.0);
        testTrancheMap.put(92.0, 10.0);
        testTrancheMap.put(100.0, -42.);
        testTrancheMap.put(90.0, 12.5);
    }

    // tests for validateFilteringCutoffs

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testValidateFilteringCutoffsOnlySnpTruthInput() {
        Double snpTruthInput = definedDoubleInput;
        Double indelTruthInput = undefinedDoubleInput;
        Double snpVqslodInput = undefinedDoubleInput;
        Double indelVqslodInput = undefinedDoubleInput;
        String trancheTableInput = definedStringInput;

        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testValidateFilteringCutoffsOnlyIndelTruthInput() {
        Double snpTruthInput = undefinedDoubleInput;
        Double indelTruthInput = definedDoubleInput;
        Double snpVqslodInput = undefinedDoubleInput;
        Double indelVqslodInput = undefinedDoubleInput;
        String trancheTableInput = definedStringInput;

        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testValidateFilteringCutoffsOnlySnpVqslodInput() {
        Double snpTruthInput = undefinedDoubleInput;
        Double indelTruthInput = undefinedDoubleInput;
        Double snpVqslodInput = definedDoubleInput;
        Double indelVqslodInput = undefinedDoubleInput;
        String trancheTableInput = definedStringInput;

        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testValidateFilteringCutoffsOnlyIndelVqslodInput() {
        Double snpTruthInput = undefinedDoubleInput;
        Double indelTruthInput = undefinedDoubleInput;
        Double snpVqslodInput = undefinedDoubleInput;
        Double indelVqslodInput = definedDoubleInput;
        String trancheTableInput = definedStringInput;

        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testValidateFilteringCutoffsBothInputs() {
        Double snpTruthInput = definedDoubleInput;
        Double indelTruthInput = definedDoubleInput;
        Double snpVqslodInput = undefinedDoubleInput;
        Double indelVqslodInput = definedDoubleInput;
        String trancheTableInput = definedStringInput;

        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testValidateFilteringCutoffsNoTranches() {
        Double snpTruthInput = undefinedDoubleInput;
        Double indelTruthInput = undefinedDoubleInput;
        Double snpVqslodInput = undefinedDoubleInput;
        Double indelVqslodInput = undefinedDoubleInput;
        String trancheTableInput = undefinedStringInput;

        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    @Test
    public void testValidateFilteringCutoffsVqslodNoTranches() {
        Double snpTruthInput = undefinedDoubleInput;
        Double indelTruthInput = undefinedDoubleInput;
        Double snpVqslodInput = definedDoubleInput;
        Double indelVqslodInput = definedDoubleInput;
        String trancheTableInput = undefinedStringInput;

        // this should pass fine - don't need tranches if vqslod cutoff is defined
        validateFilteringCutoffs(snpTruthInput, indelTruthInput, snpVqslodInput, indelVqslodInput, trancheTableInput);
    }

    // tests for getVqslodThreshold

    @Test
    public void testGetVqslodThresholdExactMatch() {
        Double testSensitivityThresh = 92.0;
        Double expectedVqslod = 10.0;
        assertEquals(getVqslodThreshold(testTrancheMap, testSensitivityThresh, GATKVCFConstants.SNP), expectedVqslod);
    }

    @Test
    public void testGetVqslodThresholdBetweenTranches() {
        Double testSensitivityThresh = 94.0;
        Double expectedVqslod = 2.5; // should use next highest tranche, i.e. 95.0
        assertEquals(getVqslodThreshold(testTrancheMap, testSensitivityThresh, GATKVCFConstants.SNP), expectedVqslod);
    }

    @Test
    public void testGetVqslodThresholdLowerThanAll() {
        Double testSensitivityThresh = 3.0;
        Double expectedVqslod = 12.5; // should use next highest tranche, i.e. 90.0
        assertEquals(getVqslodThreshold(testTrancheMap, testSensitivityThresh, GATKVCFConstants.SNP), expectedVqslod);
    }

    @Test
    public void testGetVqslodThresholdNullSNP() {
        Double testSensitivityThresh = null; // should default to 99.7 for snps
        Double expectedVqslod = -42.0;
        assertEquals(getVqslodThreshold(testTrancheMap, testSensitivityThresh, GATKVCFConstants.SNP), expectedVqslod);
    }

    @Test
    public void testGetVqslodThresholdNullINDEL() {
        Double testSensitivityThresh = null; // should default to 99.0 for indels
        Double expectedVqslod = -6.0;
        assertEquals(getVqslodThreshold(testTrancheMap, testSensitivityThresh, GATKVCFConstants.INDEL), expectedVqslod);
    }

    // tests for getVqslodThresholdFromTranches

    @Test
    public void testGetVqslodThresholdFromTranchesExactMatch() {
        Double testSensitivityThresh = 92.0;
        Double expectedVqslod = 10.0;
        assertEquals(getVqslodThresholdFromTranches(testTrancheMap, testSensitivityThresh), expectedVqslod);
    }

    @Test
    public void testGetVqslodThresholdFromTranchesBetweenTranches() {
        Double testSensitivityThresh = 94.0;
        Double expectedVqslod = 2.5; // should use next highest tranche, i.e. 95.0
        assertEquals(getVqslodThresholdFromTranches(testTrancheMap, testSensitivityThresh), expectedVqslod);
    }

    @Test
    public void testGetVqslodThresholdFromTranchesLowerThanAll() {
        Double testSensitivityThresh = 3.0;
        Double expectedVqslod = 12.5; // should use next highest tranche, i.e. 90.0
        assertEquals(getVqslodThresholdFromTranches(testTrancheMap, testSensitivityThresh), expectedVqslod);
    }

    @Test(expectedExceptions = UserException.class)
    public void testGetVqslodThresholdFromTranchesHigherThanAll() {
        Double testSensitivityThresh = 200.0;
        getVqslodThresholdFromTranches(testTrancheMap, testSensitivityThresh);
    }


    // tests for getVqsLodHeader

    @Test
    public void testGetVqsLodHeaderSNP() {
        Double vqsLodThreshold = 0.0;
        String model = GATKVCFConstants.SNP;
        VCFFilterHeaderLine expectedHeader = new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed SNP model VQSLOD cutoff of 0.0");

        assertEquals(getVqsLodHeader(vqsLodThreshold, GATKVCFConstants.SNP), expectedHeader);
    }

    @Test
    public void testGetVqsLodHeaderINDEL() {
        Double vqsLodThreshold = 0.0;
        String model = GATKVCFConstants.INDEL;
        VCFFilterHeaderLine expectedHeader = new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed INDEL model VQSLOD cutoff of 0.0");

        assertEquals(getVqsLodHeader(vqsLodThreshold, GATKVCFConstants.INDEL), expectedHeader);
    }

    // tests for getTruthSensitivityHeader

    @Test
    public void testGetTruthSensitivityHeaderSNP() {
        Double vqsLodThreshold = 0.0;
        Double truthSensitivityThreshold = 90.0;
        String model = GATKVCFConstants.SNP;
        VCFFilterHeaderLine expectedHeader = new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed SNP model sensitivity cutoff (90.0), corresponding with VQSLOD cutoff of 0.0");

        assertEquals(getTruthSensitivityHeader(truthSensitivityThreshold, vqsLodThreshold, GATKVCFConstants.SNP), expectedHeader);
    }

    @Test
    public void testGetTruthSensitivityHeaderINDEL() {
        Double vqsLodThreshold = 0.0;
        Double truthSensitivityThreshold = 90.0;
        String model = GATKVCFConstants.INDEL;
        VCFFilterHeaderLine expectedHeader = new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed INDEL model sensitivity cutoff (90.0), corresponding with VQSLOD cutoff of 0.0");

        assertEquals(getTruthSensitivityHeader(truthSensitivityThreshold, vqsLodThreshold, GATKVCFConstants.INDEL), expectedHeader);
    }

}
