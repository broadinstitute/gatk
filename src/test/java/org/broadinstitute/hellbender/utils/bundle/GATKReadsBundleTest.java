package org.broadinstitute.hellbender.utils.bundle;

import htsjdk.beta.plugin.IOUtils;
import htsjdk.beta.plugin.bundle.BundleJSON;
import org.broadinstitute.hellbender.cmdline.argumentcollections.GATKReadsBundle;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class GATKReadsBundleTest extends BaseTest {

    private final static String BAM_FILE = "reads.bam";
    private final static String INDEX_FILE = "reads.bai";

    public static final String EXPECTED_JSON_STRING =
            "{\"schemaVersion\":\"0.1.0\",\"INDEX\":{\"path\":\"reads.bai\"},\"schemaName\":\"htsbundle\",\"READS\":{\"path\":\"reads.bam\",\"subtype\":\"BAM\"},\"primary\":\"READS\"}";

    @Test
    public void testReadFromJSONString(){
        final String json = EXPECTED_JSON_STRING;
        final GATKReadsBundle gatkReadsBundleFromString = GATKReadsBundle.getGATKReadsBundleFromString(json);
        final GATKPath path = new GATKPath(BAM_FILE);

        Assert.assertTrue(gatkReadsBundleFromString.getReads().getIOPath().get() instanceof GATKPath);
        Assert.assertEquals(gatkReadsBundleFromString.getReads().getIOPath().get(), path);
    }

    @Test
    public void testWriteToJSONString() {
        final GATKReadsBundle readsBundle = new GATKReadsBundle(
                new GATKPath(BAM_FILE),
                new GATKPath(INDEX_FILE));
        final String s = BundleJSON.toJSON(readsBundle);
        Assert.assertEquals(s, EXPECTED_JSON_STRING);
    }

    @Test
    public void testReadFromJSONFile(){
        final GATKPath jsonFilePath = getTestFileGATKPath("reads1.json");
        IOUtils.writeStringToPath(jsonFilePath, EXPECTED_JSON_STRING);
        final GATKReadsBundle gatkReadsBundleFromPath = GATKReadsBundle.getGATKReadsBundleFromPath(jsonFilePath);
        Assert.assertTrue(gatkReadsBundleFromPath.getReads().getIOPath().isPresent());
        Assert.assertTrue(gatkReadsBundleFromPath.getReads().getIOPath().get() instanceof GATKPath);
        final GATKPath path = new GATKPath(BAM_FILE);
        Assert.assertEquals(gatkReadsBundleFromPath.getReads().getIOPath().get(), path);
    }

}