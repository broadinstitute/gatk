package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.TreeMap;

import static org.broadinstitute.hellbender.tools.variantdb.nextgen.FilterSensitivityTools.*;
import static org.testng.Assert.*;

public class FilterSensitivityToolsTest {

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




    // TODO how to test the exceptions in this method?
//    @Test
//    public void testValidateFilteringCutoffs() {
//        String definedInput = "I'm defined!";
//        String undefinedInput = null;
//
//        assertThrows(new testValidateFilteringCutoffs(definedInput, undefinedInput, undefinedInput, undefinedInput, definedInput));
//    }

    @Test
    public void testGetTrancheMaps() {

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

    // tests for getVqsLodHeader

    @Test
    public void testGetVqsLodHeader() {
        Double vqsLodSNPThreshold = 90.0;
        String model = GATKVCFConstants.SNP;
        VCFFilterHeaderLine expectedHeader = new VCFFilterHeaderLine(GATKVCFConstants.VQSR_FAILURE_PREFIX + model,
                "Site failed SNP model VQSLOD cutoff of 90.0");

        assertEquals(getVqsLodHeader(vqsLodSNPThreshold, GATKVCFConstants.SNP), expectedHeader);
    }

    @Test
    public void testGetTruthSensitivityHeader() {
    }

}
