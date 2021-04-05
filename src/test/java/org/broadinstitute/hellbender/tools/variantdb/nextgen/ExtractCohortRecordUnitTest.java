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
        String sampleName = "SM-GXZVY";
        Assert.assertEquals(allDefinedRecord.getSampleName(), sampleName);
        Assert.assertEquals(someNullsRecord.getSampleName(), sampleName);
    }

    @Test
    public void testGetState() {
        Assert.assertEquals(allDefinedRecord.getState(), "v");

        Assert.assertEquals(someNullsRecord.getState(), "3");
    }

    @Test
    public void testGetRefAllele() {
        String ref = "CTTT";
        Assert.assertEquals(allDefinedRecord.getRefAllele(), ref);

        Assert.assertNull(someNullsRecord.getRefAllele());
    }

    @Test
    public void testGetAltAllele() {
        String alt = "C";
        Assert.assertEquals(allDefinedRecord.getAltAllele(), alt);

        Assert.assertNull(someNullsRecord.getAltAllele());
    }

    @Test
    public void testGetCallGT() {
        String callGT = "1/1";
        Assert.assertEquals(allDefinedRecord.getCallGT(), callGT);

        Assert.assertNull(someNullsRecord.getCallGT());
    }

    @Test
    public void testGetCallGQ() {
        String callGQ = "3";
        Assert.assertEquals(allDefinedRecord.getCallGQ(), callGQ);

        Assert.assertNull(someNullsRecord.getCallGQ());
    }

    @Test
    public void testGetCallRGQ() {
        String callRGQ = "38";
        Assert.assertEquals(allDefinedRecord.getCallRGQ(), callRGQ);

        Assert.assertNull(someNullsRecord.getCallRGQ());
    }

    @Test
    public void testGetQUALApprox() {
        String qualapprox = "38";
        Assert.assertEquals(allDefinedRecord.getQUALApprox(), qualapprox);

        Assert.assertNull(someNullsRecord.getQUALApprox());
    }

    @Test
    public void testGetAsQUALApprox() {
        String asQualapprox = "38";
        Assert.assertEquals(allDefinedRecord.getAsQUALApprox(), asQualapprox);

        Assert.assertNull(someNullsRecord.getAsQUALApprox());
    }

    @Test
    public void testGetCallPL() {
        String callPL = "38,3,0,31,3,28";
        Assert.assertEquals(allDefinedRecord.getCallPL(), callPL);

        Assert.assertNull(someNullsRecord.getCallPL());
    }
}
