package org.broadinstitute.hellbender.tools.variantdb.nextgen;

//import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

//import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortRecordUnitTest extends GATKBaseTest {

//    private GenericRecord loadGenericRecord() {
//         // TODO: read in a real avro file, create a GenericRecord input
//        final String testAvroFileName =
//        final AvroFileReader avroFileReader = new AvroFileReader(testAvroFileName);
//    }

    private final Long location = Long.parseLong("21000008416376");
    private final String sampleName = "28";
    private final String state = "v";
    private final String ref = "G";
    private final String alt = "C";
    private final String callGT = "1|1";
    private final String callGQ = "3";
    private final String callRGQ = "???RGQ";
    private final String qualapprox = "19";
    private final String asQualapprox = "19";
    private final String callPL = "???PL";

    private ExtractCohortRecord allDefinedRecord = new ExtractCohortRecord(location, sampleName,
            state, ref, alt, callGT, callGQ, callRGQ, qualapprox, asQualapprox, callPL);

    private ExtractCohortRecord someNullsRecord = new ExtractCohortRecord(location, sampleName,
            state, ref, alt, null, null, null, null, null, null);

    @Test
    public void testGetContig() {
        String expectedContig = "chr21";
        Assert.assertEquals(allDefinedRecord.getContig(), expectedContig);
        Assert.assertEquals(someNullsRecord.getContig(), expectedContig);
    }

    @Test
    public void testGetStart() {
        int expectedStart = 8416376;
        Assert.assertEquals(allDefinedRecord.getStart(), expectedStart);
        Assert.assertEquals(someNullsRecord.getStart(), expectedStart);
    }

    @Test
    public void testGetEnd() {
        int expectedEnd = 8416376;
        Assert.assertEquals(allDefinedRecord.getEnd(), expectedEnd);
        Assert.assertEquals(someNullsRecord.getEnd(), expectedEnd);
    }

    @Test
    public void testGetLocation() {
        Assert.assertEquals(allDefinedRecord.getLocation(), location);
        Assert.assertEquals(someNullsRecord.getLocation(), location);
    }

    @Test
    public void testGetSampleName() {
        Assert.assertEquals(allDefinedRecord.getSampleName(), sampleName);
        Assert.assertEquals(someNullsRecord.getSampleName(), sampleName);
    }

    @Test
    public void testGetState() {
        Assert.assertEquals(allDefinedRecord.getState(), state);
        Assert.assertEquals(someNullsRecord.getState(), state);
    }

    @Test
    public void testGetRefAllele() {
        Assert.assertEquals(allDefinedRecord.getRefAllele(), ref);
        Assert.assertEquals(someNullsRecord.getRefAllele(), ref);
    }

    @Test
    public void testGetAltAllele() {
        Assert.assertEquals(allDefinedRecord.getAltAllele(), alt);
        Assert.assertEquals(someNullsRecord.getAltAllele(), alt);
    }

    @Test
    public void testGetCallGT() {
        Assert.assertEquals(allDefinedRecord.getCallGT(), callGT);
        Assert.assertNull(someNullsRecord.getCallGT());
    }

    @Test
    public void testGetCallGQ() {
        Assert.assertEquals(allDefinedRecord.getCallGQ(), callGQ);
        Assert.assertNull(someNullsRecord.getCallGQ());
    }

    @Test
    public void testGetCallRGQ() {
        Assert.assertEquals(allDefinedRecord.getCallRGQ(), callRGQ);
        Assert.assertNull(someNullsRecord.getCallRGQ());
    }

    @Test
    public void testGetQUALApprox() {
        Assert.assertEquals(allDefinedRecord.getQUALApprox(), qualapprox);
        Assert.assertNull(someNullsRecord.getQUALApprox());
    }

    @Test
    public void testGetAsQUALApprox() {
        Assert.assertEquals(allDefinedRecord.getAsQUALApprox(), asQualapprox);
        Assert.assertNull(someNullsRecord.getAsQUALApprox());
    }

    @Test
    public void testGetCallPL() {
        Assert.assertEquals(allDefinedRecord.getCallPL(), callPL);
        Assert.assertNull(someNullsRecord.getCallPL());
    }
}
