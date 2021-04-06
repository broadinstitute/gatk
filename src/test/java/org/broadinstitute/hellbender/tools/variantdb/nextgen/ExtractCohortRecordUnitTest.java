package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortRecordUnitTest extends GATKBaseTest {

    @Test
    public void testDefinedExtractCohortRecord() {
        GenericRecord definedInputGenericRecord = new AvroFileReader(getToolTestDataDir() + "test_input_defined.avro").next();
        ExtractCohortRecord allDefinedRecord = new ExtractCohortRecord(definedInputGenericRecord);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr15");
        Assert.assertEquals(allDefinedRecord.getStart(), 28865305);
        Assert.assertEquals(allDefinedRecord.getEnd(), 28865305);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("15000028865305"));
        Assert.assertEquals(allDefinedRecord.getSampleName(), "SM-GXZVY");
        Assert.assertEquals(allDefinedRecord.getState(), "v");
        Assert.assertEquals(allDefinedRecord.getRefAllele(), "CTTT");
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "C");
        Assert.assertEquals(allDefinedRecord.getCallGT(), "1/1");
        Assert.assertEquals(allDefinedRecord.getCallGQ(), "3");
        Assert.assertEquals(allDefinedRecord.getCallRGQ(), "38");
        Assert.assertEquals(allDefinedRecord.getQUALApprox(), "38");
        Assert.assertEquals(allDefinedRecord.getAsQUALApprox(), "38");
        Assert.assertEquals(allDefinedRecord.getCallPL(), "38,3,0,31,3,28");
    }

    @Test
    public void testNulledExtractCohortRecord() {
        GenericRecord nulledInputGenericRecord = new AvroFileReader(getToolTestDataDir() + "test_input_nulls.avro").next();
        ExtractCohortRecord someNullsRecord = new ExtractCohortRecord(nulledInputGenericRecord);

        Assert.assertEquals(someNullsRecord.getContig(), "chr15");
        Assert.assertEquals(someNullsRecord.getStart(), 28865344);
        Assert.assertEquals(someNullsRecord.getEnd(), 28865344);
        Assert.assertEquals(someNullsRecord.getLocation(), Long.parseLong("15000028865344"));
        Assert.assertEquals(someNullsRecord.getSampleName(), "SM-GXZVY");
        Assert.assertEquals(someNullsRecord.getState(), "3");

        // the rest are null
        Assert.assertNull(someNullsRecord.getRefAllele());
        Assert.assertNull(someNullsRecord.getAltAllele());
        Assert.assertNull(someNullsRecord.getCallGT());
        Assert.assertNull(someNullsRecord.getCallGQ());
        Assert.assertNull(someNullsRecord.getCallRGQ());
        Assert.assertNull(someNullsRecord.getQUALApprox());
        Assert.assertNull(someNullsRecord.getAsQUALApprox());
        Assert.assertNull(someNullsRecord.getCallPL());
    }
}
