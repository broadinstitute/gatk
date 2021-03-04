package org.broadinstitute.hellbender.utils.bigquery;


import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.apache.avro.Schema;

import java.util.*;

/**
 * A class to test the functionality of {@link AvroFileReader}.
 */
public class AvroFileReaderUnitTest extends GATKBaseTest {

    private static final String avroFileURI = "gs://broad-dsp-spec-ops/kcibul/acmg.35.chr20.1mb.cohort.extract.avro";
    private static final AvroFileReader AvroFile = new AvroFileReader(avroFileURI);

    @Test()
    public void testGetSchema() {
        Schema AvroFileSchema = AvroFile.getSchema();
        String  TestAvroFileSchema  = "{\"type\":\"record\",\"name\":\"Root\",\"fields\":[{\"name\":\"location\",\"type\":[\"null\",\"long\"]},{\"name\":\"sample_name\",\"type\":[\"null\",\"string\"]},{\"name\":\"state\",\"type\":[\"null\",\"string\"]},{\"name\":\"ref\",\"type\":[\"null\",\"string\"]},{\"name\":\"alt\",\"type\":[\"null\",\"string\"]},{\"name\":\"call_GT\",\"type\":[\"null\",\"string\"]},{\"name\":\"call_GQ\",\"type\":[\"null\",\"long\"]},{\"name\":\"call_RGQ\",\"type\":[\"null\",\"long\"]},{\"name\":\"QUALapprox\",\"type\":[\"null\",\"string\"]},{\"name\":\"AS_QUALapprox\",\"type\":[\"null\",\"string\"]},{\"name\":\"call_PL\",\"type\":[\"null\",\"string\"]}]}";
        Assert.assertEquals(AvroFileSchema.toString(), TestAvroFileSchema, "AvroFileSchema did not match.");
    }

    @Test()
    public void testAvroFileHasNext() {
        boolean AvroFileHasNext = AvroFile.hasNext();
        Assert.assertTrue(AvroFileHasNext, "Arvo File Reader didn't detect");
    }

    @Test()
    public void testAvroFileNext() {
        GenericRecord AvroFileNext = AvroFile.next();
        String  TestAvroFileNext  =   "{\"location\": 20000000060998, \"sample_name\": \"SM-GXZUY\", \"state\": \"0\", \"ref\": null, \"alt\": null, \"call_GT\": null, \"call_GQ\": null, \"call_RGQ\": null, \"QUALapprox\": null, \"AS_QUALapprox\": null, \"call_PL\": null}";
        Assert.assertEquals(AvroFileNext.toString(), TestAvroFileNext, "AvroFile Next did not match.");
    }

    @Test()
    public void testiterator() {
        Iterator<GenericRecord> AvroFileTest = AvroFile.iterator();
        Assert.assertEquals(AvroFile, AvroFileTest, "AvroFile did not match.");
    }

}
