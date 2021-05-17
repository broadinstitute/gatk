package org.broadinstitute.hellbender.tools.gvs.filtering;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class ExtractFeaturesEngineTest extends GATKBaseTest {

    @Test
    public void testExtractFeaturesEngineQueryLabelCorrect() {
        List<String> labelStringList = new ArrayList<String>();
        labelStringList.add("labelkey=labelvalue");
        Map<String, String> labelMap = ExtractFeaturesEngine.createQueryLabels(labelStringList);
        Assert.assertEquals(labelMap.get("gvs_tool_name"), "extract-features");
        Assert.assertEquals(labelMap.get("gvs_query_name"), "extract-features");
    }

    @Test(expectedExceptions = { UserException.class })
    public void testExtractFeaturesEngineQueryLabelTooMany() {
        List<String> labelStringList = new ArrayList<>();
        for(int i=0; i<63; i++) { // up to 64 are allowed, and 2 are static, so only 62 additional can be passed in
            labelStringList.add("labelkey"+i+"=labelvalue"+i);
        }
        ExtractFeaturesEngine.createQueryLabels(labelStringList);
    }

    @Test(expectedExceptions = { UserException.class })
    public void testExtractFeaturesEngineQueryLabelKeyTooLong() {
        List<String> labelStringList = new ArrayList<>();
        labelStringList.add("labelkey_that_is_quite_long_and_needs_to_be_shortened_in_order_to_be_acceptable_to_BQ=labelvalue");
        ExtractFeaturesEngine.createQueryLabels(labelStringList);
    }

    @Test(expectedExceptions = { UserException.class })
    public void testExtractFeaturesEngineQueryLabelValueTooLong() {
        List<String> labelStringList = new ArrayList<>();
        labelStringList.add("labelkey=labelvalue_that_is_quite_long_and_needs_to_be_shortened_in_order_to_be_acceptable_to_BQ");
        ExtractFeaturesEngine.createQueryLabels(labelStringList);
    }

    @Test(expectedExceptions = { UserException.class })
    public void testExtractFeaturesEngineQueryLabelSpaces() {
        List<String> labelStringList = new ArrayList<>();
        labelStringList.add("label key=label value");
        ExtractFeaturesEngine.createQueryLabels(labelStringList);
    }

    @Test(expectedExceptions = { UserException.class })
    public void testExtractFeaturesEngineQueryLabelCaps() {
        List<String> labelStringList = new ArrayList<>();
        labelStringList.add("labelKey=labelValue");
        ExtractFeaturesEngine.createQueryLabels(labelStringList);
    }

    @Test(expectedExceptions = { UserException.class })
    public void testExtractFeaturesEngineQueryLabelBadDelimiter() {
        List<String> labelStringList = new ArrayList<>();
        labelStringList.add("labelkey:labelvalue");
        ExtractFeaturesEngine.createQueryLabels(labelStringList);
    }
}
