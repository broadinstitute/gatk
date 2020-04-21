package org.broadinstitute.hellbender.utils.bundle;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ReadsBundleTest extends BaseTest {

    @Test
    public void testWrite() {
        final ReadsBundle readsBundle = new ReadsBundle(new BundleFile("a file", "bam"), new BundleFile("an index", "bai"));
        final String s = readsBundle.toJson();
        Assert.assertEquals(s, "{\n" +
                "  \"reads\": {\n" +
                "    \"path\": \"a file\",\n" +
                "    \"fileType\": \"bam\"\n" +
                "  },\n" +
                "  \"index\": {\n" +
                "    \"path\": \"an index\",\n" +
                "    \"fileType\": \"bai\"\n" +
                "  },\n" +
                "  \"schemaVersion\": \"0.1.0\"\n" +
                "}");
    }

    @Test
    public void testRead(){
        final String json = "{\"reads\":{\"path\":\"a file\",\"fileType\":\"bam\"},\"index\":{\"path\":\"an index\",\"fileType\":\"bai\"}}";
        final ReadsBundle readsBundle1 = ReadsBundle.fromJson(json);
        Assert.assertEquals(readsBundle1.getReads(), "a file");
    }

    @Test
    public void testReadFromFile(){
        final GATKPath json = getTestPath("reads1.json");
        final ReadsBundle readsBundle1 = ReadsBundle.fromPath(json);
        Assert.assertEquals(readsBundle1.getReads(), "a file");
    }

}