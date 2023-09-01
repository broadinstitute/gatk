package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class TrimmedReadsReaderUnitTest extends GATKBaseTest {

    private static String testDir = publicTestDir + FlowTestConstants.VARIANT_CALLING_DATA_DIR;

    private static class BamSource {
        final String      filename;
        final int         readCount;
        final String      firstReadName;
        final String      lastReadName;

        BamSource(final String filename, final int readCount, final String firstReadName, final String lastReadName) {
            this.filename = filename;
            this.readCount = readCount;
            this.firstReadName = firstReadName;
            this.lastReadName = lastReadName;
        }
    }

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        final Object[][]        testData = {
                {
                    new BamSource[] {
                            new BamSource("chr5.bam1.rename.bam", 324, "2265833312", "0608076045"),
                            new BamSource("chr5.bam2.rename.bam", 299, "2507511136", "1255139282")
                    },
                    new SimpleInterval("chr5:70036483-70036764"),
                    new SimpleInterval("chr5:70036625-70036625")
                }
        };

        return testData;
    }

    @Test(dataProvider = "testData")
    public void testBasic(final BamSource bamSources[], final Locatable span, final Locatable vcLoc) throws Exception {

        // establish paths
        List<Path>  paths = new LinkedList<>();
        for ( int i = 0 ; i < bamSources.length ; i++ ) {
            paths.add(FileSystems.getDefault().getPath(testDir, bamSources[i].filename));
        }

        // create reader
        TrimmedReadsReader reader = new TrimmedReadsReader(paths, null, 40);

        // global
        Assert.assertNotNull(reader.getSamSequenceDictionary(null));
        Assert.assertNotNull(reader.getHeader(null));

        // reads
        Map<SamReader, Collection<FlowBasedRead>> reads = reader.getReads(span, vcLoc);
        Assert.assertEquals(reads.size(), bamSources.length);
        int             bamSourceIndex = 0;
        for ( Map.Entry<SamReader, Collection<FlowBasedRead>> entry : reads.entrySet() ) {

            final BamSource       bamSource = bamSources[bamSourceIndex++];

            // sanity
            Assert.assertNotNull(reader.getSamSequenceDictionary(entry.getKey()));
            Assert.assertNotNull(reader.getHeader(entry.getKey()));

            // verify size
            Assert.assertEquals(entry.getValue().size(), bamSource.readCount);

            // verify first and last
            FlowBasedRead       firstRead = entry.getValue().iterator().next();
            FlowBasedRead       lastRead = entry.getValue().stream().reduce((prev, next) -> next).orElse(null);
            Assert.assertEquals(firstRead.getName(), bamSource.firstReadName);
            Assert.assertEquals(lastRead.getName(), bamSource.lastReadName);
        }
    }
}
