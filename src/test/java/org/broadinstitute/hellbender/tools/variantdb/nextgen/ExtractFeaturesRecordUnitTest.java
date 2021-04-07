package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.bigquery.AvroFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractFeaturesRecordUnitTest extends GATKBaseTest {

    @Test
    public void testExtractFeaturesRecordDefined() {
        GenericRecord inputGenericRecord = new AvroFileReader(getToolTestDataDir() + "extract_features_input_defined.avro").next();
        ExtractFeaturesRecord allDefinedRecord = new ExtractFeaturesRecord(inputGenericRecord);

        Assert.assertEquals(allDefinedRecord.getContig(), "chr10");
        Assert.assertEquals(allDefinedRecord.getStart(), 112719503);
        Assert.assertEquals(allDefinedRecord.getEnd(), 112719503);
        Assert.assertEquals(allDefinedRecord.getLocation(), Long.parseLong("10000112719503"));
        Assert.assertEquals(allDefinedRecord.getRef(), "G");
        Assert.assertEquals(allDefinedRecord.getAllele(), "T");
        Assert.assertEquals(allDefinedRecord.getRawQual(), Double.valueOf("15272"));
        Assert.assertEquals(allDefinedRecord.getRefAD(), Double.valueOf("361"));
        Assert.assertEquals(allDefinedRecord.getAsMQRankSum(), Float.valueOf("0.0"));
        Assert.assertEquals(allDefinedRecord.getAsMQRankSumFreqTable(), "[{\"value\": -9, \"freq\": 1}, {\"value\": 0, \"freq\": 16}]");
        Assert.assertEquals(allDefinedRecord.getAsReadPosRankSum(), Float.valueOf("-0.2"));
        Assert.assertEquals(allDefinedRecord.getAsReadPosRankSumFreqTable(), "[{\"value\": -16, \"freq\": 1}, {\"value\": -9, \"freq\": 1}, {\"value\": -7, \"freq\": 1}, {\"value\": -5, \"freq\": 2}, {\"value\": -3, \"freq\": 2}, {\"value\": -2, \"freq\": 2}, {\"value\": -1, \"freq\": 1}, {\"value\": 0, \"freq\": 1}, {\"value\": 1, \"freq\": 1}, {\"value\": 3, \"freq\": 1}, {\"value\": 7, \"freq\": 1}, {\"value\": 9, \"freq\": 1}, {\"value\": 13, \"freq\": 1}, {\"value\": 14, \"freq\": 1}]");
        Assert.assertEquals(allDefinedRecord.getRawMQ(), Double.valueOf("1431225"));
        Assert.assertEquals(allDefinedRecord.getRawAD(), Double.valueOf("398"));
        Assert.assertEquals(allDefinedRecord.getRawADGT1(), Double.valueOf("398"));
        Assert.assertEquals(allDefinedRecord.getSbRefPlus(), 205);
        Assert.assertEquals(allDefinedRecord.getSbRefMinus(), 156);
        Assert.assertEquals(allDefinedRecord.getSbAltPlus(), 221);
        Assert.assertEquals(allDefinedRecord.getSbAltMinus(), 177);
        Assert.assertEquals(allDefinedRecord.getNumHetSamples(), 17); // int
        Assert.assertEquals(allDefinedRecord.getNumHomvarSamples(), 0); // int
    }

    @Test
    public void testExtractFeaturesRecordNulled() {
        GenericRecord inputGenericRecord = new AvroFileReader(getToolTestDataDir() + "extract_features_input_nulls.avro").next();
        ExtractFeaturesRecord someNullsRecord = new ExtractFeaturesRecord(inputGenericRecord);

        System.out.println(someNullsRecord.getAsMQRankSumFreqTable());

        Assert.assertEquals(someNullsRecord.getContig(), "chr14");
        Assert.assertEquals(someNullsRecord.getStart(), 64822351);
        Assert.assertEquals(someNullsRecord.getEnd(), 64822351);
        Assert.assertEquals(someNullsRecord.getLocation(), Long.parseLong("14000064822351"));
        Assert.assertEquals(someNullsRecord.getRef(), "TTCTC");
        Assert.assertEquals(someNullsRecord.getAllele(), "T");
        Assert.assertEquals(someNullsRecord.getRawQual(), Double.valueOf("1"));
        Assert.assertEquals(someNullsRecord.getRefAD(), Double.valueOf("0"));
        Assert.assertNull(someNullsRecord.getAsMQRankSum());
        Assert.assertEquals(someNullsRecord.getAsMQRankSumFreqTable(), "[]");
        Assert.assertNull(someNullsRecord.getAsReadPosRankSum());
        Assert.assertEquals(someNullsRecord.getAsReadPosRankSumFreqTable(), "[]");
        Assert.assertEquals(someNullsRecord.getRawMQ(), Double.valueOf("0"));
        Assert.assertEquals(someNullsRecord.getRawAD(), Double.valueOf("0"));
        Assert.assertEquals(someNullsRecord.getRawADGT1(), Double.valueOf("0"));
        Assert.assertEquals(someNullsRecord.getSbRefPlus(), 1);
        Assert.assertEquals(someNullsRecord.getSbRefMinus(), 2);
        Assert.assertEquals(someNullsRecord.getSbAltPlus(), 0);
        Assert.assertEquals(someNullsRecord.getSbAltMinus(), 0);
        Assert.assertEquals(someNullsRecord.getNumHetSamples(), 2);
        Assert.assertEquals(someNullsRecord.getNumHomvarSamples(), 2);
    }
}
