package org.broadinstitute.hellbender.utils;


// the imports for unit testing.


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class NGSPlatformUnitTest extends GATKBaseTest {

    // example fasta index file, can be deleted if you don't use the reference
    private IndexedFastaSequenceFile seq;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(IOUtils.getPath(exampleReference));
    }

    @DataProvider(name = "TestPrimary")
    public Object[][] makeTestPrimary() {
        List<Object[]> tests = new ArrayList<>();

        for ( final NGSPlatform pl : NGSPlatform.values() ) {
            tests.add(new Object[]{pl, pl.BAM_PL_NAMES[0]});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TestPrimary")
    public void testPrimary(final NGSPlatform pl, final String expectedPrimaryName) {
        Assert.assertEquals(pl.getDefaultPlatform(), expectedPrimaryName, "Failed primary test for " + pl);
    }

    // make sure common names in BAMs are found
    @DataProvider(name = "TestMappings")
    public Object[][] makeTestMappings() {
        List<Object[]> tests = new ArrayList<>();

        final Map<String, NGSPlatform> expected = new LinkedHashMap<>();
        // VALID VALUES ACCORDING TO SAM SPEC: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CC8QFjAA&url=http%3A%2F%2Fsamtools.sourceforge.net%2FSAM1.pdf&ei=Dm8WUbXAEsi10QHYqoDwDQ&usg=AFQjCNFkMtvEi6LeiKgpxQGtHTlqWKw2yw&bvm=bv.42080656,d.dmQ
        expected.put("CAPILLARY", NGSPlatform.CAPILLARY);
        expected.put("LS454", NGSPlatform.LS454);
        expected.put("ILLUMINA", NGSPlatform.ILLUMINA);
        expected.put("SOLID", NGSPlatform.SOLID);
        expected.put("HELICOS", NGSPlatform.HELICOS);
        expected.put("IONTORRENT", NGSPlatform.ION_TORRENT);
        expected.put("PACBIO", NGSPlatform.PACBIO);
        // other commonly seen values out in the wild
        expected.put("SLX", NGSPlatform.ILLUMINA);
        expected.put("SOLEXA", NGSPlatform.ILLUMINA);
        expected.put("454", NGSPlatform.LS454);
        expected.put("COMPLETE", NGSPlatform.COMPLETE_GENOMICS);
        // unknown platforms should map to unknown
        expected.put("MARKS_GENOMICS_TECH", NGSPlatform.UNKNOWN);
        expected.put("RANDOM_PL_VALUE", NGSPlatform.UNKNOWN);
        // critical -- a null platform maps to unknown
        expected.put(null, NGSPlatform.UNKNOWN);

        for ( final Map.Entry<String,NGSPlatform> one : expected.entrySet() ) {
            tests.add(new Object[]{one.getKey(), one.getValue()});

            if ( one.getKey() != null ) {
                // make sure we're case insensitive
                tests.add(new Object[]{one.getKey().toLowerCase(), one.getValue()});
                tests.add(new Object[]{one.getKey().toUpperCase(), one.getValue()});

                // make sure appending GENOMICS works (required for COMPLETE mapping
                tests.add(new Object[]{one.getKey() + " GENOMICS", one.getValue()});
                // make sure that random junk works correctly
                tests.add(new Object[]{one.getKey() + " asdfa", one.getValue()});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TestMappings")
    public void testMappings(final String plField, final NGSPlatform expected) {
        Assert.assertEquals(NGSPlatform.fromReadGroupPL(plField), expected, "Failed primary test for " + plField + " mapping to " + expected);
    }

    @Test(dataProvider = "TestMappings")
    public void testKnown(final String plField, final NGSPlatform expected) {
        Assert.assertEquals(NGSPlatform.isKnown(plField), expected != NGSPlatform.UNKNOWN, "Failed isKnown test for " + plField + " mapping to " + expected);
    }

    /**
     * A unit test that creates an artificial read for testing some code that uses reads
     */
    @Test(dataProvider = "TestMappings")
    public void testPLFromReadWithRG(final String plField, final NGSPlatform expected) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final String rgID = "ID";
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(rgID);
        if ( plField != null )
            rg.setPlatform(plField);
        header.addReadGroup(rg);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        read.setAttribute("RG", rgID);
        Assert.assertEquals(NGSPlatform.fromRead(read, header), expected);
    }

    @Test()
    public void testPLFromReadWithRGButNoPL() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final String rgID = "ID";
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(rgID);
        header.addReadGroup(rg);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        read.setAttribute("RG", rgID);
        Assert.assertEquals(NGSPlatform.fromRead(read, header), NGSPlatform.UNKNOWN);
    }

    @Test
    public void testReadWithoutRG() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        Assert.assertEquals(NGSPlatform.fromRead(read, header), NGSPlatform.UNKNOWN);
    }
}