package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortRecordUnitTest extends GATKBaseTest {

    private static GenericRecord getSingleRecordFromAvroFile (String inputFileName) {
        final AvroFileReader avroFileReader = new AvroFileReader(inputFileName);
        return avroFileReader.next();
    }

    // load in test data
    private static final String testFileDir = "src/test/resources/org/broadinstitute/hellbender/tools/variantdb/nextgen/ExtractCohortRecord/";
    private static final GenericRecord definedInputGenericRecord = getSingleRecordFromAvroFile(testFileDir + "test_input_defined.avro");
    private static final GenericRecord nulledInputGenericRecord = getSingleRecordFromAvroFile(testFileDir + "test_input_nulls.avro");

    // generate records to test
    private ExtractCohortRecord allDefinedRecord = new ExtractCohortRecord(definedInputGenericRecord);
    private ExtractCohortRecord someNullsRecord = new ExtractCohortRecord(nulledInputGenericRecord);

    @Test
    public void testGetContig() {
        String expectedContig = "chr15";
        Assert.assertEquals(allDefinedRecord.getContig(), expectedContig);
        Assert.assertEquals(someNullsRecord.getContig(), expectedContig);
    }

    @Test
    public void testGetStart() {
        int expectedStartDefined = 28865305;
        Assert.assertEquals(allDefinedRecord.getStart(), expectedStartDefined);

        int expectedStartNulled = 28865344;
        Assert.assertEquals(someNullsRecord.getStart(), expectedStartNulled);
    }

    @Test
    public void testGetEnd() {
        int expectedEndDefined = 28865305;
        Assert.assertEquals(allDefinedRecord.getEnd(), expectedEndDefined);


        int expectedEndNulled = 28865344;
        Assert.assertEquals(someNullsRecord.getEnd(), expectedEndNulled);
    }

    @Test
    public void testGetLocation() {
        Long expectedLocationDefined = Long.parseLong("15000028865305");
        Assert.assertEquals(allDefinedRecord.getLocation(), expectedLocationDefined);

        Long expectedLocationNulled = Long.parseLong("15000028865344");
        Assert.assertEquals(someNullsRecord.getLocation(), expectedLocationNulled);
    }

    @Test
    public void testGetSampleName() {
        String expectedSampleName = "SM-GXZVY";
        Assert.assertEquals(allDefinedRecord.getSampleName(), expectedSampleName);
        Assert.assertEquals(someNullsRecord.getSampleName(), expectedSampleName);
    }

    @Test
    public void testGetState() {
        Assert.assertEquals(allDefinedRecord.getState(), "v");

        Assert.assertEquals(someNullsRecord.getState(), "3");
    }

    @Test
    public void testGetRefAllele() {
        String expectedRef = "CTTT";
        Assert.assertEquals(allDefinedRecord.getRefAllele(), expectedRef);

        Assert.assertNull(someNullsRecord.getRefAllele());
    }

    @Test
    public void testGetAltAllele() {
        String expectedAlt = "C";
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "C");

        Assert.assertNull(someNullsRecord.getAltAllele());
    }

    @Test
    public void testGetCallGT() {
        String expectedCallGT = "1/1";
        Assert.assertEquals(allDefinedRecord.getCallGT(), expectedCallGT);

        Assert.assertNull(someNullsRecord.getCallGT());
    }

    @Test
    public void testGetCallGQ() {
        String expectedCallGQ = "3";
        Assert.assertEquals(allDefinedRecord.getCallGQ(), expectedCallGQ);

        Assert.assertNull(someNullsRecord.getCallGQ());
    }

    @Test
    public void testGetCallRGQ() {
        String expectedCallRGQ = "38";
        Assert.assertEquals(allDefinedRecord.getCallRGQ(), expectedCallRGQ);

        Assert.assertNull(someNullsRecord.getCallRGQ());
    }

    @Test
    public void testGetQUALApprox() {
        String expectedQualapprox = "38";
        Assert.assertEquals(allDefinedRecord.getQUALApprox(), expectedQualapprox);

        Assert.assertNull(someNullsRecord.getQUALApprox());
    }

    @Test
    public void testGetAsQUALApprox() {
        String expectedAsQualapprox = "38";
        Assert.assertEquals(allDefinedRecord.getAsQUALApprox(), expectedAsQualapprox);

        Assert.assertNull(someNullsRecord.getAsQUALApprox());
    }

    @Test
    public void testGetCallPL() {
        String expectedCallPL = "38,3,0,31,3,28";
        Assert.assertEquals(allDefinedRecord.getCallPL(), expectedCallPL);

        Assert.assertNull(someNullsRecord.getCallPL());
    }
}
