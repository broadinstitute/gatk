package org.broadinstitute.hellbender.utils.bundle;

import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadsBundle;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ReadsBundleTest extends BaseTest {

    private final static String BAM_FILE = "reads.bam";
    private final static String INDEX_FILE = "reads.bai";

    public static final String EXPECTED_JSON_STRING =
            "{\"schemaVersion\":\"0.1.0\",\"INDEX\":{\"path\":\"reads.bai\"},\"schemaName\":\"htsbundle\",\"READS\":{\"path\":\"reads.bam\"}}";

    @Test
    public void testWrite() {
        final ReadsBundle readsBundle = new ReadsBundle(
                new GATKPath(BAM_FILE),
                new GATKPath(INDEX_FILE));
        final String s = readsBundle.toJSON();
        Assert.assertEquals(s, EXPECTED_JSON_STRING);
    }

    @Test
    public void testRead(){
        final String json = EXPECTED_JSON_STRING;
        final ReadsBundle readsBundle1 = new ReadsBundle(json);
        //TODO: verify that this space gets correctly encoded ?
        final GATKPath path = new GATKPath(BAM_FILE);
        //TODO: propagate the tags and attributes
        //path.setTagAttributes(Collections.singletonMap(BundleResourceType.READS_BAM, "bam"));
        Assert.assertEquals(readsBundle1.getReads(), path);
    }

    @Test
    public void testReadFromFile(){
        final GATKPath jsonFilePath = getTestFileGATKPath("reads1.json");
        final ReadsBundle readsBundle1 = ReadsBundle.getReadsBundleFromJSON(jsonFilePath);
        final GATKPath path = new GATKPath(BAM_FILE);
        //path.setTagAttributes(Collections.singletonMap(ReadsBundle.FILE_TYPE_KEY, "bam"));
        Assert.assertEquals(readsBundle1.getReads(), path);
    }

}