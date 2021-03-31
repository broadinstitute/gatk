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

    private static final String avroFileName = "src/test/java/org/broadinstitute/hellbender/utils/bigquery/avro_test_avro_test_file.avro";
    private static final AvroFileReader avroFile = new AvroFileReader(avroFileName);

    @Test()
    public void testGetSchema() {
        Schema avroFileSchema = avroFile.getSchema();
        String  testAvroFileSchema  = "{\"type\":\"record\",\"name\":\"Root\",\"fields\":[{\"name\":\"test_string\",\"type\":[\"null\",\"string\"]},{\"name\":\"test_float\",\"type\":[\"null\",\"double\"]},{\"name\":\"test_integer\",\"type\":[\"null\",\"long\"]}]}";
        Assert.assertEquals(avroFileSchema.toString(), testAvroFileSchema, "AvroFileSchema did not match.");
    }

    @Test()
    public void testAvroFileHasNext() {
        boolean avroFileHasNext = avroFile.hasNext();
        Assert.assertTrue(avroFileHasNext, "Arvo File Reader didn't detect");
    }

    @Test()
    public void testAvroFileNext() {
        GenericRecord avroFileNext = avroFile.next();
        String  testAvroFileNext  = "{\"test_string\": \"one\", \"test_float\": 1111111.0, \"test_integer\": 1}";
        Assert.assertEquals(avroFileNext.toString(), testAvroFileNext, "AvroFile Next did not match.");
    }

}
