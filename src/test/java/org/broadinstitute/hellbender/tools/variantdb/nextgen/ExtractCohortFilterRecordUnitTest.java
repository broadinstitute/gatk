package org.broadinstitute.hellbender.tools.variantdb.nextgen;

//import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

//import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortFilterRecordUnitTest extends GATKBaseTest {

//    private GenericRecord loadGenericRecord() {
//         // TODO: read in a real avro file, create a GenericRecord input
//        final String testAvroFileName =
//        final AvroFileReader avroFileReader = new AvroFileReader(testAvroFileName);
//    }

    private final Long location = Long.parseLong("21000008416376");
    private final String sampleName = "28";
    private final String ref = "G";
    private final String alt = "C";


    private final String vqslod = "100000";
    private final String yng_status = "G";


    private ExtractCohortFilterRecord allDefinedRecord = new ExtractCohortFilterRecord(location,
            sampleName, ref, alt, vqslod, yng_status);

    private ExtractCohortFilterRecord someNullsRecord = new ExtractCohortFilterRecord(location,
            sampleName, null, null, vqslod, yng_status);

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
    public void testGetRefAllele() {
        Assert.assertEquals(allDefinedRecord.getRefAllele(), ref);
        Assert.assertNull(someNullsRecord.getRefAllele());
    }

    @Test
    public void testGetAltAllele() {
        Assert.assertEquals(allDefinedRecord.getAltAllele(), alt);
        Assert.assertNull(someNullsRecord.getAltAllele());
    }

    @Test
    public void testGetVqslod() {
        Double expectedVqslod = Double.parseDouble(vqslod);
        Assert.assertEquals(allDefinedRecord.getVqslod(), expectedVqslod);
        Assert.assertEquals(someNullsRecord.getVqslod(), expectedVqslod);
    }

    @Test
    public void testGetYng() {
        Assert.assertEquals(allDefinedRecord.getYng(), yng_status);
        Assert.assertEquals(someNullsRecord.getYng(), yng_status);
    }
}
