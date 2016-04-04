package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Unit tests for {@link AllelicCountCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    private static final File SNPS_FILE = new File(TEST_SUB_DIR, "snps.tsv");
    private static final File SNPS_WITH_MISSING_COLUMN_FILE = new File(TEST_SUB_DIR, "snps-with-missing-column.tsv");

    @Test(expectedExceptions = UserException.class)
    public void testReadWithMissingColumns() throws Exception {
        new AllelicCountCollection(SNPS_WITH_MISSING_COLUMN_FILE);
    }

    @Test
    public void testAddAndGetCounts() throws Exception {
        final AllelicCountCollection counts = new AllelicCountCollection();
        final AllelicCount count = new AllelicCount(new SimpleInterval("1", 1, 1), 1, 1);
        counts.add(new SimpleInterval("1", 1, 1), 1, 1);

        Assert.assertEquals(count, counts.getCounts().get(0));
    }

    @Test
    public void testReadAndWrite() throws Exception {
        final File tempFile = createTempFile("allelic-count-collection-test", "tsv");
        final AllelicCountCollection snps = new AllelicCountCollection(SNPS_FILE);
        snps.write(tempFile);

        final AllelicCountCollection snpsWritten = new AllelicCountCollection(tempFile);

        final AllelicCountCollection snpsExpected = new AllelicCountCollection();
        snpsExpected.add(new SimpleInterval("1", 881918, 881918), 14, 21);
        snpsExpected.add(new SimpleInterval("1", 909238, 909238), 13, 11);
        snpsExpected.add(new SimpleInterval("1", 934940, 934940), 14, 0);
        snpsExpected.add(new SimpleInterval("1", 949608, 949608), 20, 14);

        Assert.assertEquals(snpsWritten, snpsExpected);
    }
}