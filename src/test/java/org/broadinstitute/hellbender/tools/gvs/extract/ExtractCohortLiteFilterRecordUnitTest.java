package org.broadinstitute.hellbender.tools.gvs.extract;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.gvs.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortLiteFilterRecordUnitTest extends GATKBaseTest {

    @Test
    public void testExtractCohortFilterRecord() {
        GenericRecord inputGenericRecord = new AvroFileReader(new GATKPath(getToolTestDataDir() + "test_input.avro")).next();
        ExtractCohortLiteFilterRecord allDefinedRecord = new ExtractCohortLiteFilterRecord(inputGenericRecord);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr1");
        Assert.assertEquals(allDefinedRecord.getStart(), 183706);
        Assert.assertEquals(allDefinedRecord.getEnd(), 183706);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("1000000183706"));
        Assert.assertEquals(allDefinedRecord.getRefAllele(), "G");
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "GT");
        Assert.assertEquals(allDefinedRecord.getSensitivity(), Double.parseDouble("0.9389"));
        Assert.assertEquals(allDefinedRecord.getYng(), "G");
    }
}
