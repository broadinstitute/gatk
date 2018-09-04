package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMTag;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;


public class ReadQueryNameComparatorUnitTest extends GATKBaseTest {

    public static final SAMFileHeader HEADER =ArtificialReadUtils.createArtificialSamHeader();
    public static final String NAME = "NAME";

    /**
     * Tests that the ordering produced by {@link ReadQueryNameComparator} matches queryname ordering
     * as produced by htsjdk's {@link SAMRecordQueryNameComparator} for a representative selection of reads. Ignores
     * differences in tie-breaking done for reads with the same position -- just asserts that the reads are
     * queryname-sorted according to htsjdk, including unmapped reads with and without an assigned position.
     */
    @Test
    public void testComparatorOrderingMatchesHtsjdkFileOrdering() throws IOException {
        final String inputBam = publicTestDir + "org/broadinstitute/hellbender/utils/read/comparator_test_with_unmapped.bam";
        final List<GATKRead> reads = new ArrayList<>();
        SAMFileHeader header;

        try ( final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(inputBam)) ) {
            header = readsSource.getHeader();

            for ( GATKRead read : readsSource ) {
                reads.add(read);
            }
        }

        // Randomize ordering and then re-sort
        Collections.shuffle(reads);
        reads.sort(new ReadQueryNameComparator());

        final SAMRecordQueryNameComparator samComparator = new SAMRecordQueryNameComparator();
        GATKRead previousRead = null;
        for ( final GATKRead currentRead : reads ) {
            if ( previousRead != null ) {
                Assert.assertTrue(samComparator.compare(previousRead.convertToSAMRecord(header), currentRead.convertToSAMRecord(header)) <= 0,
                                  "Reads are out of order: " + previousRead + " and " + currentRead);
            }
            previousRead = currentRead;
        }
    }

    @DataProvider
    public Object[][] getNames(){
        return new Object[][]{
                {"A", "B", -1},
                {"A","A", 0},
                {"AA", "A", 1},
                {"1","10", -1},
                {"2", "10", 1}
        };
    }



    @Test(dataProvider = "getNames")
    public void testCompareNames(String firstName, String secondName, int expected) throws Exception {
        ReadQueryNameComparator comp = new ReadQueryNameComparator();
        GATKRead first = getRead(firstName);
        GATKRead second = getRead(secondName);
        Assert.assertEquals(comp.compareReadNames(first, second ), expected);
        Assert.assertEquals(comp.compareReadNames(second, first), -expected);
        Assert.assertEquals(comp.compareReadNames(first, first), 0);
        Assert.assertEquals(comp.compareReadNames(second, second), 0);
    }

    private static GATKRead getRead(String firstName) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(HEADER, firstName, 1, 100, 10);
        return read;
    }

    @DataProvider
    public Iterator<Object[]> getReads(){
        final GATKRead differentName = getRead(NAME+NAME);

        final GATKRead unpaired = getRead(NAME);
        unpaired.setIsPaired(false);

        final GATKRead paired = getRead(NAME);
        paired.setIsPaired(true);

        final GATKRead firstOfPair = getRead(NAME);
        firstOfPair.setIsFirstOfPair();

        final GATKRead secondOfPair = getRead(NAME);
        secondOfPair.setIsSecondOfPair();

        final GATKRead reverseStrand = getRead(NAME);
        reverseStrand.setIsReverseStrand(true);

        final GATKRead supplementary = getRead(NAME);
        supplementary.setIsSupplementaryAlignment(true);

        final GATKRead secondary = getRead(NAME);
        secondary.setIsSecondaryAlignment(true);

        final GATKRead tagHI1 = getRead(NAME);
        tagHI1.setAttribute(SAMTag.HI.name(), 1);

        final GATKRead tagHI2 = getRead(NAME);
        tagHI2.setAttribute(SAMTag.HI.name(), 2);

        List<GATKRead> reads = Arrays.asList(differentName, unpaired, paired, firstOfPair, secondOfPair, reverseStrand, supplementary, secondary, tagHI1, tagHI2);
        List<Object[]> tests = new ArrayList<>();
        for(GATKRead left: reads){
            for(GATKRead right: reads){
                tests.add(new Object[]{left, right});
            }
        };
        return tests.iterator();
    }

    @Test(dataProvider = "getReads")
    public void testTieBreakers(GATKRead left, GATKRead right){
        ReadQueryNameComparator readComparator = new ReadQueryNameComparator();
        SAMRecordQueryNameComparator samComparator = new SAMRecordQueryNameComparator();
        Assert.assertEquals(readComparator.compare(left, right), samComparator.compare(left.convertToSAMRecord(HEADER), right.convertToSAMRecord(HEADER)));
    }

}