package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class TextMDCodecTest {

    @DataProvider(name = "mdStringProvider")
    public Object[][] getMdStrings() {
        return new Object[][]{
                { "151", Arrays.asList(new TextMDCodec.MatchMDElement(151)) },
                { "75^AC76", Arrays.asList(new TextMDCodec.MatchMDElement(75), new TextMDCodec.DeletionMDElement(2), new TextMDCodec.MatchMDElement(76)) },
                // these were seen in the wild so need to support them
                { "0AT151", Arrays.asList(new TextMDCodec.MatchMDElement(0), new TextMDCodec.MismatchMDElement(), new TextMDCodec.MismatchMDElement(), new TextMDCodec.MatchMDElement(151)) }
        };
    }

    @Test(dataProvider = "mdStringProvider")
    public void testParseMDString(String mdString, List<TextMDCodec.MDElement> expectedElements) throws Exception {
        final List<TextMDCodec.MDElement> actualElements = TextMDCodec.parseMDString(mdString);
        Assert.assertEquals(actualElements.size(), expectedElements.size());
        for (int i = 0; i < actualElements.size(); i++) {
            Assert.assertEquals(actualElements.get(i), expectedElements.get(i));
        }

    }

}