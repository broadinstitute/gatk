package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class StripMateNumberTransformerTest extends GATKBaseTest {


    @DataProvider(name = "testData")
    public Object[][] getTestData() {
        return new Object[][]{
                {"TEST_READ/1", "TEST_READ"},
                {"TEST_READ/2", "TEST_READ"},
                {"TEST_READ/1 ", "TEST_READ"},
                {"TEST_READ/2  ", "TEST_READ"},
                {"TEST_READ /1", "TEST_READ"},
                {"TEST_READ /2", "TEST_READ"},
                {"TEST_READ/1/2", "TEST_READ/1"},
                {"TEST_READ/2/1", "TEST_READ/2"},
                {"TEST/1READ", "TEST/1READ"},
                {"TEST/2READ", "TEST/2READ"},
                {"TEST/1READ/1", "TEST/1READ"},
                {"TEST/2READ/1", "TEST/2READ"},
        };
    }

    @Test(dataProvider = "testData")
    public void test(final String nameIn, final String nameOut) {
        final StripMateNumberTransformer trans = new StripMateNumberTransformer();
        final GATKRead read = new SAMRecordToGATKReadAdapter(new SAMRecord(null));
        read.setName(nameIn);
        Assert.assertEquals(trans.apply(read).getName(), nameOut);
    }

}