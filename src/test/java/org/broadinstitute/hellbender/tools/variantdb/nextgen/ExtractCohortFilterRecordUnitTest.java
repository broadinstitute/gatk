package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortFilterRecordUnitTest extends GATKBaseTest {

    // TODO refactor - this is copied from ExtractCohortRecordUnitTest
    private static GenericRecord getSingleRecordFromAvroFile (String inputFileName) {
        final AvroFileReader avroFileReader = new AvroFileReader(inputFileName);
        return avroFileReader.next();
    }

    // load in test data
    private static final String testFileDir = "src/test/resources/org/broadinstitute/hellbender/tools/variantdb/nextgen/ExtractCohortFilterRecord/";
    private static final GenericRecord inputGenericRecord = getSingleRecordFromAvroFile(testFileDir + "test_input.avro");


    private ExtractCohortFilterRecord allDefinedRecord = new ExtractCohortFilterRecord(inputGenericRecord);

    @Test
    public void testGetContig() {
        String expectedContig = "chr1";
        Assert.assertEquals(allDefinedRecord.getContig(), expectedContig);
    }

    @Test
    public void testGetStart() {
        int expectedStart = 183706;
        Assert.assertEquals(allDefinedRecord.getStart(), expectedStart);
    }

    @Test
    public void testGetEnd() {
        int expectedEnd = 183706;
        Assert.assertEquals(allDefinedRecord.getEnd(), expectedEnd);
    }

    @Test
    public void testGetLocation() {
        Long expectedLocation = Long.parseLong("1000000183706");
        Assert.assertEquals(allDefinedRecord.getLocation(), expectedLocation);
    }

    @Test
    public void testGetRefAllele() {
        String expectedRef = "G";
        Assert.assertEquals(allDefinedRecord.getRefAllele(), expectedRef);
    }

    @Test
    public void testGetAltAllele() {
        String expectedAlt = "GT";
        Assert.assertEquals(allDefinedRecord.getAltAllele(), expectedAlt);
    }

    @Test
    public void testGetVqslod() {
        Double expectedVqslod = Double.parseDouble("20.2295");
        Assert.assertEquals(allDefinedRecord.getVqslod(), expectedVqslod);
    }

    @Test
    public void testGetYng() {
        String expectedYng = "G";
        Assert.assertEquals(allDefinedRecord.getYng(), expectedYng);
    }
}
