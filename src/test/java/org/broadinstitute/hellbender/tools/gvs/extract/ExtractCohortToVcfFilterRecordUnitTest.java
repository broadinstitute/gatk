package org.broadinstitute.hellbender.tools.gvs.extract;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.gvs.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractCohortToVcfFilterRecordUnitTest extends GATKBaseTest {

    @Test
    public void testExtractCohortFilterVQSRClassicRecord() {
        GenericRecord inputGenericRecord = new AvroFileReader(new GATKPath(getToolTestDataDir() + "test_input_vqsr_classic.avro")).next();
        ExtractCohortFilterRecord allDefinedRecord = new ExtractCohortFilterRecord(inputGenericRecord, SchemaUtils.VQSLOD, null);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr1");
        Assert.assertEquals(allDefinedRecord.getStart(), 183706);
        Assert.assertEquals(allDefinedRecord.getEnd(), 183706);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("1000000183706"));
        Assert.assertEquals(allDefinedRecord.getRefAllele(), "G");
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "GT");
        Assert.assertEquals(allDefinedRecord.getVqScore(), Double.parseDouble("20.2295"));
        Assert.assertEquals(allDefinedRecord.getYng(), "G");
    }

    @Test
    public void testExtractCohortFilterVQSRLiteRecord() {
        GenericRecord inputGenericRecord = new AvroFileReader(new GATKPath(getToolTestDataDir() + "test_input_vqsr_lite.avro")).next();
        ExtractCohortFilterRecord allDefinedRecord = new ExtractCohortFilterRecord(inputGenericRecord, SchemaUtils.CALIBRATION_SENSITIVITY, SchemaUtils.SCORE);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr1");
        Assert.assertEquals(allDefinedRecord.getStart(), 183706);
        Assert.assertEquals(allDefinedRecord.getEnd(), 183706);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("1000000183706"));
        Assert.assertEquals(allDefinedRecord.getRefAllele(), "G");
        Assert.assertEquals(allDefinedRecord.getAltAllele(), "GT");
        Assert.assertEquals(allDefinedRecord.getVqScore(), Double.parseDouble("0.9389"));
        Assert.assertEquals(allDefinedRecord.getYng(), "G");
    }
}
