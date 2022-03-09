package org.broadinstitute.hellbender.tools.gvs.extract;

import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortRecordUnitTest extends GATKBaseTest {
    // KCIBUL - 2021-06-04 For the historical record, I had to recreate these AVROs for this PR and wanted to document my reverse
    // engineering.  This query was run, saving results to a temp table and then using the BQ GUI to export that to Avro
    //
    //     SELECT location, 34, state, ref, alt, call_GT, call_GQ, call_RGQ, QUALapprox, AS_QUALapprox, call_PL
    //     FROM `spec-ops-aou.gvs_tieout_acmg_v2.exported_cohort_all_samples`
    //     WHERE location = XXX and sample_name = 'SM-GXZVY'
    // using 15000028865305 as the location for the "defined" test and 15000028865344 for the "null" test

    @Test
    public void testDefinedExtractCohortRecord() {
        GenericRecord definedInputGenericRecord = new AvroFileReader(new GATKPath(getToolTestDataDir() + "test_input_defined.avro")).next();
        ExtractCohortRecord allDefinedRecord = new ExtractCohortRecord(definedInputGenericRecord);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr15");
        Assert.assertEquals(allDefinedRecord.getStart(), 28865305);
        Assert.assertEquals(allDefinedRecord.getEnd(), 28865305);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("15000028865305"));
        Assert.assertEquals(allDefinedRecord.getSampleId().longValue(), 34);
        Assert.assertEquals(allDefinedRecord.getState(), "v");
        Assert.assertEquals(allDefinedRecord.getRefAllele(), "CTTT");
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "C");
        Assert.assertEquals(allDefinedRecord.getCallGT(), "1/1");
        // Assert.assertEquals(allDefinedRecord.getCallAD(), "0,1"); // TODO need to add this to the extract Avro!
        Assert.assertEquals(allDefinedRecord.getCallGQ(), "3");
        Assert.assertEquals(allDefinedRecord.getCallRGQ(), "38");
        Assert.assertEquals(allDefinedRecord.getQUALApprox(), "38");
        Assert.assertEquals(allDefinedRecord.getAsQUALApprox(), "38");
        Assert.assertEquals(allDefinedRecord.getCallPL(), "38,3,0,31,3,28");
    }

    @Test
    public void testNulledExtractCohortRecord() {
        GenericRecord nulledInputGenericRecord = new AvroFileReader(new GATKPath(getToolTestDataDir() + "test_input_nulls.avro")).next();
        ExtractCohortRecord someNullsRecord = new ExtractCohortRecord(nulledInputGenericRecord);

        Assert.assertEquals(someNullsRecord.getContig(), "chr15");
        Assert.assertEquals(someNullsRecord.getStart(), 28865344);
        Assert.assertEquals(someNullsRecord.getEnd(), 28865344);
        Assert.assertEquals(someNullsRecord.getLocation(), Long.parseLong("15000028865344"));
        Assert.assertEquals(someNullsRecord.getSampleId().longValue(), 34);
        Assert.assertEquals(someNullsRecord.getState(), "3");

        // the rest are null
        Assert.assertNull(someNullsRecord.getRefAllele());
        Assert.assertNull(someNullsRecord.getAltAllele());
        Assert.assertNull(someNullsRecord.getCallGT());
        Assert.assertNull(someNullsRecord.getCallAD());
        Assert.assertNull(someNullsRecord.getCallGQ());
        Assert.assertNull(someNullsRecord.getCallRGQ());
        Assert.assertNull(someNullsRecord.getQUALApprox());
        Assert.assertNull(someNullsRecord.getAsQUALApprox());
        Assert.assertNull(someNullsRecord.getCallPL());
    }
}
