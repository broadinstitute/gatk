package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractFeaturesRecordUnitTest extends GATKBaseTest {

    @Test
    public void testExtractFeaturesRecordDefined() {
        GenericRecord inputGenericRecord = new AvroFileReader(getToolTestDataDir() + "feature_extract_test_defined.avro").next();
        ExtractFeaturesRecord allDefinedRecord = new ExtractFeaturesRecord(inputGenericRecord);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr1");
        Assert.assertEquals(allDefinedRecord.getStart(), 112719503);
        Assert.assertEquals(allDefinedRecord.getEnd(), 112719503);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("10000112719503"));
        Assert.assertEquals(allDefinedRecord.getRef(), "G");
        Assert.assertEquals(allDefinedRecord.getAllele(), "T");
        Assert.assertEquals(allDefinedRecord.getRawQual(), Double.valueOf("15272")); // Double
        Assert.assertEquals(allDefinedRecord.getRefAD(), "361");
        Assert.assertEquals(allDefinedRecord.getAsMQRankSum(), Float.valueOf("0.0")); // Float
        Assert.assertEquals(allDefinedRecord.getAsReadPosRankSum(), Float.valueOf("-0.2")); // Float
        Assert.assertEquals(allDefinedRecord.getRawMQ(), Double.valueOf("1431225")); // Double
        Assert.assertEquals(allDefinedRecord.getRawAD(), Double.valueOf("398")); // Double
        Assert.assertEquals(allDefinedRecord.getRawADGT1(), Double.valueOf("361")); // Double
        Assert.assertEquals(allDefinedRecord.getSbRefPlus(), 205); // int
        Assert.assertEquals(allDefinedRecord.getSbRefMinus(), 156); // int
        Assert.assertEquals(allDefinedRecord.getSbAltPlus(), 221); // int
        Assert.assertEquals(allDefinedRecord.getSbAltMinus(), 177); // int
//        Assert.assertEquals(allDefinedRecord.getNumHetSamples(), ""); // int
//        Assert.assertEquals(allDefinedRecord.getNumHomvarSamples(), ""); // int
    }

    @Test
    public void testExtractFeaturesRecordNulled() {
        GenericRecord inputGenericRecord = new AvroFileReader(getToolTestDataDir() + "feature_extract_test_nulls.avro").next();
        ExtractFeaturesRecord someNullsRecord = new ExtractFeaturesRecord(inputGenericRecord);

        Assert.assertEquals(someNullsRecord.getContig(), "chr1");
        Assert.assertEquals(someNullsRecord.getStart(), 183706);
        Assert.assertEquals(someNullsRecord.getEnd(), 183706);
        Assert.assertEquals(someNullsRecord.getLocation(), Long.parseLong("1000000183706"));
        Assert.assertEquals(someNullsRecord.getRef(), "G");
        Assert.assertEquals(someNullsRecord.getAllele(), "GT");
        Assert.assertEquals(someNullsRecord.getRawQual(), ""); // Double
        Assert.assertEquals(someNullsRecord.getRefAD(), "");
        Assert.assertEquals(someNullsRecord.getAsMQRankSum(), ""); // Float
        Assert.assertEquals(someNullsRecord.getAsReadPosRankSum(), ""); // Float
        Assert.assertEquals(someNullsRecord.getRawMQ(), ""); // Double
        Assert.assertEquals(someNullsRecord.getRawAD(), ""); // Double
        Assert.assertEquals(someNullsRecord.getRawADGT1(), ""); // Double
        Assert.assertEquals(someNullsRecord.getSbRefPlus(), ""); // int
        Assert.assertEquals(someNullsRecord.getSbRefMinus(), ""); // int
        Assert.assertEquals(someNullsRecord.getSbAltPlus(), ""); // int
        Assert.assertEquals(someNullsRecord.getSbAltMinus(), ""); // int
        Assert.assertEquals(someNullsRecord.getNumHetSamples(), ""); // int
        Assert.assertEquals(someNullsRecord.getNumHomvarSamples(), ""); // int
    }
}
