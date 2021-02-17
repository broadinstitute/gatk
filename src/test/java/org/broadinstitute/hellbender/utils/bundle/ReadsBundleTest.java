package org.broadinstitute.hellbender.utils.bundle;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;

public class ReadsBundleTest extends BaseTest {

    public static final String EXPECTED_JSON_STRING = "{\n" +
            "  \"reads\": {\n" +
            "    \"path\": \"a file\",\n" +
            "    \"fileType\": \"bam\"\n" +
            "  },\n" +
            "  \"index\": {\n" +
            "    \"path\": \"an index\",\n" +
            "    \"fileType\": \"bai\"\n" +
            "  },\n" +
            "  \"schemaVersion\": \"0.1.0\",\n" +
            "  \"schemaType\": \"ReadsBundle\"\n" +
            "}";

    @Test
    public void testWrite() {
        final ReadsBundle readsBundle = new ReadsBundle(new BundleFile("a file", "bam"), new BundleFile("an index", "bai"));
        final String s = readsBundle.toJson();
        Assert.assertEquals(s, EXPECTED_JSON_STRING);
    }

    @Test
    public void testRead(){
        final String json = EXPECTED_JSON_STRING;
        final ReadsBundle readsBundle1 = ReadsBundle.fromJson(json);
        final GATKPath path = new GATKPath("a file");
        path.setTagAttributes(Collections.singletonMap(ReadsBundle.FILE_TYPE_KEY, "bam"));
        Assert.assertEquals(readsBundle1.getReads(), path);
    }

    @Test
    public void testReadFromFile(){
        final GATKPath json = getTestFileGATKPath("reads1.json");
        final ReadsBundle readsBundle1 = ReadsBundle.fromPath(json);
        final GATKPath path = new GATKPath("a file");
        path.setTagAttributes(Collections.singletonMap(ReadsBundle.FILE_TYPE_KEY, "bam"));
        Assert.assertEquals(readsBundle1.getReads(), path);
    }

}