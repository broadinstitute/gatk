package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.Random;

public class ReadsKeyUnitTest extends GATKBaseTest {

   private GATKRead createTestRead(String name, String contig, int startPosition, String cigar, String readGroupRecord, boolean reverse) {
       GATKRead read = ArtificialReadUtils.createSamBackedRead(name, contig, startPosition, 100);
       read.setCigar(cigar);
       read.setReadGroup(readGroupRecord);
       read.setIsReverseStrand(reverse);
       return read;
   }

    @Test (enabled = false)
    // This test makes the argument that using a StringBuilder in ReadsKey is faster than an equivalent call to String.format()
    /* an example of the results:
        With Format:	13909867289
        With Builder:	1130077051
        With String +:	1155343164
        With Joiner:	2013227121
     */
    public void testKeyForFragmentPerformanceEvidence() throws Exception {
        Random rand = new Random();
        int numTrials = 10000000;
        for (int i = 0; i < numTrials; i++) {
            rand.nextInt();
        }

        long startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = String.format(
                    "%s|%d|%d|%s",
                    rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-",
                    rand.nextInt(30),
                    rand.nextInt(100),
                    rand.nextBoolean() ? "r" : "f");
        }
        long endTime = System.nanoTime();
        System.out.println("With Format:\t" + (endTime - startTime));

        startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = new StringBuilder().append(rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-")
                    .append("|")
                    .append(rand.nextInt(30))
                    .append("|")
                    .append(rand.nextInt(100))
                    .append("|")
                    .append( rand.nextBoolean() ? "r" : "f")
                    .toString();

        }
        endTime = System.nanoTime();
        System.out.println("With Builder:\t" + (endTime - startTime));

        startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = (rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-") +
                    "|" +
                    rand.nextInt(30) +
                    "|" +
                    rand.nextInt(100) +
                    "|" +
                    (rand.nextBoolean() ? "r" : "f");

        }
        endTime = System.nanoTime();
        System.out.println("With String +:\t" + (endTime - startTime));

        startTime = System.nanoTime();
        for (int i = 0; i < numTrials; i++) {
            String s = String.join("|", rand.nextBoolean() ? "SomeLongTestLibraryName-----x" : "-",
                    Integer.toString(rand.nextInt(30)),
                    Integer.toString(rand.nextInt(100)),
                    rand.nextBoolean() ? "r" : "f");

        }
        endTime = System.nanoTime();
        System.out.println("With Joiner:\t" + (endTime - startTime));
    }

    @DataProvider
    public Object[][] artificalReadsForKeys() {
        SAMReadGroupRecord library1 = new SAMReadGroupRecord("rg1");
        library1.setLibrary("library1");
        SAMReadGroupRecord library2 = new SAMReadGroupRecord("rg2");
        library2.setLibrary("library2");

        SAMFileHeader header = hg19Header.clone();
        header.addReadGroup(library1);
        header.addReadGroup(library2);

        return new Object[][]{
                // Test of two groups with different start positions
                {header, createTestRead("name1", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        false,
                        createTestRead("name2", "1", 1010, "100M", library1.getReadGroupId(), false),
                        createTestRead("name2", "1", 1200, "100M", library1.getReadGroupId(), true),},

                // Test of two equivalent groups
                {header, createTestRead("name1", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        true,
                        createTestRead("name2", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name2", "1", 1200, "100M", library1.getReadGroupId(), true),},

                // Test of two equivalent group, different contig
                {header, createTestRead("name1", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        false,
                        createTestRead("name2", "2", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name2", "2", 1200, "100M", library1.getReadGroupId(), true),},

                // Test of two equivalent groups, different orientation
                {header, createTestRead("name1", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        false,
                        createTestRead("name2", "1", 1000, "100M", library1.getReadGroupId(), true),
                        createTestRead("name2", "1", 1200, "100M", library1.getReadGroupId(), true),},

                // Test of two equivalent grops, but different libraries
                {header, createTestRead("name1", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        false,
                        createTestRead("name2", "1", 1000, "100M", library2.getReadGroupId(), false),
                        createTestRead("name2", "1", 1200, "100M", library2.getReadGroupId(), true),},

                // Test of two equivalent groups, but one was soft-clipped to a different start
                {header, createTestRead("name1", "1", 1010, "10S90M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        true,
                        createTestRead("name2", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name2", "1", 1200, "100M", library1.getReadGroupId(), true),},

                // Test of two equivalent groups, but one read is on a different contig
                {header, createTestRead("name1", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name1", "1", 1200, "100M", library1.getReadGroupId(), true),
                        false,
                        createTestRead("name2", "1", 1000, "100M", library1.getReadGroupId(), false),
                        createTestRead("name2", "2", 1200, "100M", library1.getReadGroupId(), true),},

        };
    }

    @Test(dataProvider = "artificalReadsForKeys")
    public void testReadsHaveDifferentKeys(SAMFileHeader header, GATKRead pair1r1, GATKRead pair1r2,
                                           boolean shouldEqual, GATKRead pair2r1, GATKRead pair2r2) {
        Map<String, Byte> libraryIndex = MarkDuplicatesSparkUtils.constructLibraryIndex(header);
        ReadsKey key1 = ReadsKey.getKeyForPair(header, pair1r1, pair1r2, libraryIndex);
        ReadsKey key2 = ReadsKey.getKeyForPair(header, pair2r1, pair2r2, libraryIndex);

        Assert.assertEquals(key1.equals(key2), shouldEqual);
    }

}