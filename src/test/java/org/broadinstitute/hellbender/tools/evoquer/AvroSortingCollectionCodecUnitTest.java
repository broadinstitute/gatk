package org.broadinstitute.hellbender.tools.evoquer;

import org.apache.avro.Schema;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

public class AvroSortingCollectionCodecUnitTest extends GATKBaseTest {

    @Test
    public void testAvroSortingCollectionCodecRoundtrip() {
        Schema.Parser schemaParser = new Schema.Parser();
        final Schema schema = schemaParser.parse("{" +
                "\"type\": \"record\"," +
                "\"name\": \"Pair\"," +
                "\"doc\": \"A pair of strings.\"," +
                "\"fields\": [" +
                "{\"name\": \"first\", \"type\": \"string\"}," +
                "{\"name\": \"second\", \"type\": \"string\"}" +
                "]" +
                "}");

        final GenericRecord record1 = new GenericData.Record(schema);
        record1.put("first", "hello");
        record1.put("second", "world");

        final GenericRecord record2 = new GenericData.Record(schema);
        record2.put("first", "number");
        record2.put("second", "two");

        Assert.assertEquals(record1.get("first").toString(), "hello");
        Assert.assertEquals(record1.get("second").toString(), "world");

        final AvroSortingCollectionCodec codec = new AvroSortingCollectionCodec(schema);

        final File tempFile = createTempFile("testAvroSortingCollectionCodecRoundtrip", ".avro");
        try ( final FileOutputStream outStream = new FileOutputStream(tempFile)) {
            codec.setOutputStream(outStream);
            codec.encode(record1);
            codec.encode(record2);
            codec.flushOutput();
        } catch ( IOException e ) {
            Assert.fail("Error writing to avro output stream", e);
        }

        try ( final FileInputStream inStream = new FileInputStream(tempFile) ) {
            codec.setInputStream(inStream);
            final GenericRecord roundtrippedRecord1 = codec.decode();

            Assert.assertNotNull(roundtrippedRecord1);
            Assert.assertEquals(roundtrippedRecord1.get("first").toString(), record1.get("first").toString());
            Assert.assertEquals(roundtrippedRecord1.get("second").toString(), record1.get("second").toString());

            final GenericRecord roundtrippedRecord2 = codec.decode();

            Assert.assertNotNull(roundtrippedRecord2);
            Assert.assertEquals(roundtrippedRecord2.get("first").toString(), record2.get("first").toString());
            Assert.assertEquals(roundtrippedRecord2.get("second").toString(), record2.get("second").toString());

            Assert.assertNull(codec.decode());
        } catch ( IOException e ) {
            Assert.fail("Error reading from avro input stream", e);
        }
    }
}
