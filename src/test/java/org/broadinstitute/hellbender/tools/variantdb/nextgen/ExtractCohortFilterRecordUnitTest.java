package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortFilterRecordUnitTest extends GATKBaseTest {

    @Test
    public void testDefinedRecords() {
        GenericRecord inputGenericRecord = new AvroFileReader(getToolTestDataDir() + "test_input.avro").next();
        ExtractCohortFilterRecord allDefinedRecord = new ExtractCohortFilterRecord(inputGenericRecord);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr1");
        Assert.assertEquals(allDefinedRecord.getStart(), 183706);
        Assert.assertEquals(allDefinedRecord.getEnd(), 183706);

        // TODO figure out how to inline this
        Long expectedLocation = Long.parseLong("1000000183706");
        Assert.assertEquals(allDefinedRecord.getLocation(), expectedLocation);

        Assert.assertEquals(allDefinedRecord.getRefAllele(), "G");
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "GT");

        // TODO figure out how to inline this
        Double expectedVqslod = Double.parseDouble("20.2295");
        Assert.assertEquals(allDefinedRecord.getVqslod(), expectedVqslod);
        
        Assert.assertEquals(allDefinedRecord.getYng(), "G");
    }
}
