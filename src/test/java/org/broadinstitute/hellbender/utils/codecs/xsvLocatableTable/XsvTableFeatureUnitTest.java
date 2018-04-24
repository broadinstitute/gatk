package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * A Class to hold unit tests for {@link XsvTableFeature}.
 * Created by jonn on 12/18/17.
 */
public class XsvTableFeatureUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final String TEST_RESOURCE_DIR = publicTestDir + "org/broadinstitute/hellbender/utils/codecs/xsvLocatableTable" + File.separator;

    private static final List<String> file1Headers = Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_chr", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_start", "XSV_LOCATABLE_TEST_NAME_end", "XSV_LOCATABLE_TEST_NAME_Bond");
    private static final List<String> file1Line1 = Arrays.asList("Blofeld", "chr19", "test_val_chr19", "8959519", "9092018", "Connery");

    final XsvTableFeature exampleFile1Feature = new XsvTableFeature(1, 3,4, file1Headers, file1Line1, "XSV_LOCATABLE_TEST_NAME");

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    //==================================================================================================================
    // Tests:

    @Test
    public void testGetStart() {
        Assert.assertEquals(exampleFile1Feature.getStart(), 8959519);
    }

    @Test
    public void testGetEnd() {
        Assert.assertEquals(exampleFile1Feature.getEnd(), 9092018);
    }

    @Test
    public void testGetContig() {
        Assert.assertEquals(exampleFile1Feature.getContig(), "chr19");
    }

    @Test
    public void testGetDataSourceName() {
        Assert.assertEquals(exampleFile1Feature.getDataSourceName(), "XSV_LOCATABLE_TEST_NAME");
    }

    @Test
    public void testGet_string_based() {
        Assert.assertEquals(exampleFile1Feature.get("XSV_LOCATABLE_TEST_NAME_Villain"), "Blofeld");
        Assert.assertEquals(exampleFile1Feature.get("XSV_LOCATABLE_TEST_NAME_chr"), "chr19");
        Assert.assertEquals(exampleFile1Feature.get("XSV_LOCATABLE_TEST_NAME_test_val"), "test_val_chr19");
        Assert.assertEquals(exampleFile1Feature.get("XSV_LOCATABLE_TEST_NAME_start"), "8959519");
        Assert.assertEquals(exampleFile1Feature.get("XSV_LOCATABLE_TEST_NAME_end"), "9092018");
        Assert.assertEquals(exampleFile1Feature.get("XSV_LOCATABLE_TEST_NAME_Bond"), "Connery");

    }

    @Test
    public void testGet_index_based() {
        Assert.assertEquals(exampleFile1Feature.get(0), "Blofeld");
        Assert.assertEquals(exampleFile1Feature.get(1), "chr19");
        Assert.assertEquals(exampleFile1Feature.get(2), "test_val_chr19");
        Assert.assertEquals(exampleFile1Feature.get(3), "8959519");
        Assert.assertEquals(exampleFile1Feature.get(4), "9092018");
        Assert.assertEquals(exampleFile1Feature.get(5), "Connery");
    }

    @Test
    public void testGetHeader() {
        Assert.assertEquals(exampleFile1Feature.getHeader(), file1Headers);
    }

    @Test
    public void testGetHeaderWithoutLocationColumns () {
        Assert.assertEquals(
            exampleFile1Feature.getHeaderWithoutLocationColumns(),
            Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain",
                    "XSV_LOCATABLE_TEST_NAME_test_val",
                    "XSV_LOCATABLE_TEST_NAME_Bond")
        );
    }

    @Test
    public void testGetValues() {
        Assert.assertEquals(exampleFile1Feature.getValues(), file1Line1);
    }

    @Test
    public void testGetValuesWithoutLocationColumns() {
        Assert.assertEquals(
            exampleFile1Feature.getValuesWithoutLocationColumns(),
            Arrays.asList("Blofeld", "test_val_chr19", "Connery")
        );
    }

    @Test
    public void testSize() {
        Assert.assertEquals(exampleFile1Feature.size(), 6);
    }


}
