package org.broadinstitute.hellbender.utils;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.SimpleFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;
import static org.testng.Assert.assertFalse;

/**
 * Test out the functionality of the new genome loc parser
 */
public final class GenomeLocParserUnitTest extends GATKBaseTest {
    private GenomeLocParser genomeLocParser;
    private SAMFileHeader header;

    @BeforeClass
    public void init() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testGetContigIndex() {
        assertEquals(genomeLocParser.getContigIndex("blah"), -1); // should not be in the reference
    }                

    @Test
    public void testGetContigIndexValid() {
        assertEquals(genomeLocParser.getContigIndex("1"), 0); // should be in the reference
    }

    @Test(expectedExceptions=UserException.class)
    public void testGetContigInfoUnknownContig1() {
        assertEquals(null, genomeLocParser.getContigInfo("blah")); // should *not* be in the reference
    }

    @Test(expectedExceptions=UserException.class)
    public void testGetContigInfoUnknownContig2() {
        assertEquals(null, genomeLocParser.getContigInfo(null)); // should *not* be in the reference
    }

    @Test()
    public void testHasContigInfoUnknownContig1() {
        assertEquals(false, genomeLocParser.contigIsInDictionary("blah")); // should *not* be in the reference
    }

    @Test()
    public void testHasContigInfoUnknownContig2() {
        assertEquals(false, genomeLocParser.contigIsInDictionary(null)); // should *not* be in the reference
    }

    @Test()
    public void testHasContigInfoKnownContig() {
        assertEquals(true, genomeLocParser.contigIsInDictionary("1")); // should be in the reference
    }

    @Test
    public void testGetContigInfoKnownContig() {
        assertEquals(0, "1".compareTo(genomeLocParser.getContigInfo("1").getSequenceName())); // should be in the reference
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testParseBadString() {
        genomeLocParser.parseGenomeLoc("Bad:0-1");
    }

    @Test
    public void testParseUnknownSequenceLength() {
        SAMSequenceDictionary seqDict = new SAMSequenceDictionary();
        seqDict.addSequence(new SAMSequenceRecord("1", SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH));
        Assert.assertEquals(seqDict.getSequence("1").getSequenceLength(), SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH);
        GenomeLocParser myLocParser = new GenomeLocParser(seqDict);
        GenomeLoc genomeLoc = myLocParser.parseGenomeLoc("1:1-99");
        Assert.assertEquals(genomeLoc.getEnd(), 99);
    }

    @Test
    public void testContigHasColon() {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        SAMSequenceRecord rec = new SAMSequenceRecord("c:h:r1", 10);
        rec.setSequenceLength(10);
        dict.addSequence(rec);
        header.setSequenceDictionary(dict);

        final GenomeLocParser myGenomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        GenomeLoc loc = myGenomeLocParser.parseGenomeLoc("c:h:r1:4-5");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStart(), 4);
        assertEquals(loc.getStop(), 5);
    }

    @Test
    public void testParseGoodString() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1:1-10");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 10);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc1() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("1", 1, 100);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc1point5() { // in honor of VAAL!
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1:1");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc2() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("1", 1, 100);
        assertEquals("1", loc.getContig());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc3() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("1", 1);
        assertEquals("1", loc.getContig());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc4() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("1", 1);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc5() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("1", 1, 100);
        GenomeLoc copy = genomeLocParser.createGenomeLoc(loc.getContig(),loc.getStart(),loc.getStop());
        assertEquals(0, copy.getContigIndex());
        assertEquals(copy.getStop(), 100);
        assertEquals(copy.getStart(), 1);
    }

    @Test
    public void testGenomeLocPlusSign() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1:1+");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testGenomeLocParseOnlyChrome() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testGenomeLocParseOnlyBadChrome() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("12");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testGenomeLocBad() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1:1-");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testGenomeLocBad2() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1:1-500-0");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testGenomeLocBad3() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("1:1--0");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    // test out the validating methods
    @Test
    public void testValidationOfGenomeLocs() {
        assertTrue(genomeLocParser.isValidGenomeLoc("1",1,1));
        assertFalse(genomeLocParser.isValidGenomeLoc("2",1,1)); // shouldn't have an entry
        assertFalse(genomeLocParser.isValidGenomeLoc("1",1,11)); // past the end of the contig
        assertFalse(genomeLocParser.isValidGenomeLoc("1",-1,10)); // bad start
        assertFalse(genomeLocParser.isValidGenomeLoc("1",1,-2)); // bad stop
        assertTrue( genomeLocParser.isValidGenomeLoc("1",-1,2, false)); // bad stop
        assertFalse(genomeLocParser.isValidGenomeLoc("1",10,11)); // bad start, past end
        assertTrue( genomeLocParser.isValidGenomeLoc("1",10,11, false)); // bad start, past end
        assertFalse(genomeLocParser.isValidGenomeLoc("1",2,1)); // stop < start
    }

    @Test(expectedExceptions = UserException.MalformedGenomeLoc.class)
    public void testValidateGenomeLoc() {
        // bad contig index
        genomeLocParser.validateGenomeLoc("1", 1, 1, 2, false);
    }

    private static class FlankingGenomeLocTestData extends TestDataProvider {
        final GenomeLocParser parser;
        final int basePairs;
        final GenomeLoc original, flankStart, flankStop;

        private FlankingGenomeLocTestData(String name, GenomeLocParser parser, int basePairs, String original, String flankStart, String flankStop) {
            super(FlankingGenomeLocTestData.class, name);
            this.parser = parser;
            this.basePairs = basePairs;
            this.original = parse(parser, original);
            this.flankStart = flankStart == null ? null : parse(parser, flankStart);
            this.flankStop = flankStop == null ? null : parse(parser, flankStop);
        }

        private static GenomeLoc parse(GenomeLocParser parser, String str) {
            return "unmapped".equals(str) ? GenomeLoc.UNMAPPED : parser.parseGenomeLoc(str);
        }
    }

    @DataProvider(name = "flankingGenomeLocs")
    public Object[][] getFlankingGenomeLocs() {
        int contigLength = 10000;
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, contigLength);
        GenomeLocParser parser = new GenomeLocParser(header.getSequenceDictionary());

        new FlankingGenomeLocTestData("atStartBase1", parser, 1,
                "1:1", null, "1:2");

        new FlankingGenomeLocTestData("atStartBase50", parser, 50,
                "1:1", null, "1:2-51");

        new FlankingGenomeLocTestData("atStartRange50", parser, 50,
                "1:1-10", null, "1:11-60");

        new FlankingGenomeLocTestData("atEndBase1", parser, 1,
                "1:" + contigLength, "1:" + (contigLength - 1), null);

        new FlankingGenomeLocTestData("atEndBase50", parser, 50,
                "1:" + contigLength, String.format("1:%d-%d", contigLength - 50, contigLength - 1), null);

        new FlankingGenomeLocTestData("atEndRange50", parser, 50,
                String.format("1:%d-%d", contigLength - 10, contigLength),
                String.format("1:%d-%d", contigLength - 60, contigLength - 11),
                null);

        new FlankingGenomeLocTestData("nearStartBase1", parser, 1,
                "1:2", "1:1", "1:3");

        new FlankingGenomeLocTestData("nearStartRange50", parser, 50,
                "1:21-30", "1:1-20", "1:31-80");

        new FlankingGenomeLocTestData("nearEndBase1", parser, 1,
                "1:" + (contigLength - 1), "1:" + (contigLength - 2), "1:" + contigLength);

        new FlankingGenomeLocTestData("nearEndRange50", parser, 50,
                String.format("1:%d-%d", contigLength - 30, contigLength - 21),
                String.format("1:%d-%d", contigLength - 80, contigLength - 31),
                String.format("1:%d-%d", contigLength - 20, contigLength));

        new FlankingGenomeLocTestData("beyondStartBase1", parser, 1,
                "1:3", "1:2", "1:4");

        new FlankingGenomeLocTestData("beyondStartRange50", parser, 50,
                "1:101-200", "1:51-100", "1:201-250");

        new FlankingGenomeLocTestData("beyondEndBase1", parser, 1,
                "1:" + (contigLength - 3),
                "1:" + (contigLength - 4),
                "1:" + (contigLength - 2));

        new FlankingGenomeLocTestData("beyondEndRange50", parser, 50,
                String.format("1:%d-%d", contigLength - 200, contigLength - 101),
                String.format("1:%d-%d", contigLength - 250, contigLength - 201),
                String.format("1:%d-%d", contigLength - 100, contigLength - 51));

        new FlankingGenomeLocTestData("unmapped", parser, 50,
                "unmapped", null, null);

        new FlankingGenomeLocTestData("fullContig", parser, 50,
                "1", null, null);

        return FlankingGenomeLocTestData.getTests(FlankingGenomeLocTestData.class);
    }

    @Test(dataProvider = "flankingGenomeLocs")
    public void testCreateGenomeLocAtStart(FlankingGenomeLocTestData data) {
        GenomeLoc actual = data.parser.createGenomeLocAtStart(data.original, data.basePairs);
        String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                data.toString(), data.original, actual, data.flankStart);
        assertEquals(actual, data.flankStart, description);
    }

    @Test(dataProvider = "flankingGenomeLocs")
    public void testCreateGenomeLocAtStop(FlankingGenomeLocTestData data) {
        GenomeLoc actual = data.parser.createGenomeLocAtStop(data.original, data.basePairs);
        String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                data.toString(), data.original, actual, data.flankStop);
        assertEquals(actual, data.flankStop, description);
    }

    @DataProvider(name = "parseGenomeLoc")
    public Object[][] makeParsingTest() {
        final List<Object[]> tests = new LinkedList<>();

        tests.add(new Object[]{ "1:10", "1", 10 });
        tests.add(new Object[]{ "1:100", "1", 100 });
        tests.add(new Object[]{ "1:1000", "1", 1000 });
        tests.add(new Object[]{ "1:1,000", "1", 1000 });
        tests.add(new Object[]{ "1:10000", "1", 10000 });
        tests.add(new Object[]{ "1:10,000", "1", 10000 });
        tests.add(new Object[]{ "1:100000", "1", 100000 });
        tests.add(new Object[]{ "1:100,000", "1", 100000 });
        tests.add(new Object[]{ "1:1000000", "1", 1000000 });
        tests.add(new Object[]{ "1:1,000,000", "1", 1000000 });
        tests.add(new Object[]{ "1:1000,000", "1", 1000000 });
        tests.add(new Object[]{ "1:1,000000", "1", 1000000 });

        return tests.toArray(new Object[][]{});
    }

    @Test( dataProvider = "parseGenomeLoc")
    public void testParsingPositions(final String string, final String contig, final int start) {
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 10000000);
        GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        final GenomeLoc loc = genomeLocParser.parseGenomeLoc(string);
        Assert.assertEquals(loc.getContig(), contig);
        Assert.assertEquals(loc.getStart(), start);
        Assert.assertEquals(loc.getStop(), start);
    }

    @Test( )
    public void testCreationFromSAMRecord() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 1, 5);
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(read);
        Assert.assertEquals(loc.getContig(), read.getContig());
        Assert.assertEquals(loc.getContigIndex(), ReadUtils.getReferenceIndex(read, header));
        Assert.assertEquals(loc.getStart(), read.getStart());
        Assert.assertEquals(loc.getStop(), read.getEnd());
    }

    @Test( )
    public void testCreationFromSAMRecordUnmapped() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 1, 5);
        read.setIsUnmapped();
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(read);
        Assert.assertTrue(loc.isUnmapped());
    }

    @Test( )
    public void testCreationFromSAMRecordUnmappedButOnGenome() {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 1, 5);
        read.setIsUnmapped();
        read.setCigar("*");
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(read);
        Assert.assertEquals(loc.getContig(), read.getContig());
        Assert.assertEquals(loc.getContigIndex(), ReadUtils.getReferenceIndex(read, header));
        Assert.assertEquals(loc.getStart(), read.getStart());
        Assert.assertEquals(loc.getStop(), read.getStart());
    }

    @Test
    public void testCreationFromFeature() {
        final Feature feature = new SimpleFeature("1", 1, 5);
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(feature);
        Assert.assertEquals(loc.getContig(), feature.getContig());
        Assert.assertEquals(loc.getStart(), feature.getStart());
        Assert.assertEquals(loc.getStop(), feature.getEnd());
    }

    @Test
    public void testCreationFromLocatable() {
        final Locatable locatable = new SimpleInterval("1", 1, 5);
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(locatable);
        Assert.assertEquals(loc.getContig(), locatable.getContig());
        Assert.assertEquals(loc.getStart(), locatable.getStart());
        Assert.assertEquals(loc.getStop(), locatable.getEnd());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testCreationFromNullLocatable() {
        genomeLocParser.createGenomeLoc((Locatable)null);
    }

    @Test
    public void testCreationFromVariantContext() {
        final VariantContext feature = new VariantContextBuilder("x", "1", 1, 5, Arrays.asList(Allele.create("AAAAA", true))).make();
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(feature);
        Assert.assertEquals(loc.getContig(), feature.getContig());
        Assert.assertEquals(loc.getStart(), feature.getStart());
        Assert.assertEquals(loc.getStop(), feature.getEnd());
    }

    @Test
    public void testcreateGenomeLocOnContig() throws FileNotFoundException {
        final CachingIndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(IOUtils.getPath(exampleReference));
        final SAMSequenceDictionary dict = seq.getSequenceDictionary();
        final GenomeLocParser genomeLocParser = new GenomeLocParser(dict);

        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final GenomeLoc loc = genomeLocParser.createOverEntireContig(rec.getSequenceName());
            Assert.assertEquals(loc.getContig(), rec.getSequenceName());
            Assert.assertEquals(loc.getStart(), 1);
            Assert.assertEquals(loc.getStop(), rec.getSequenceLength());
        }
    }

    @DataProvider(name = "GenomeLocOnContig")
    public Object[][] makeGenomeLocOnContig() {
        final List<Object[]> tests = new LinkedList<>();

        final int contigLength = header.getSequence(0).getSequenceLength();
        for ( int start = -10; start < contigLength + 10; start++ ) {
            for ( final int len : Arrays.asList(1, 10, 20) ) {
                tests.add(new Object[]{ "1", start, start + len });
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test( dataProvider = "GenomeLocOnContig")
    public void testGenomeLocOnContig(final String contig, final int start, final int stop) {
        final int contigLength = header.getSequence(0).getSequenceLength();
        final GenomeLoc loc = genomeLocParser.createGenomeLocOnContig(contig, start, stop);

        if ( stop < 1 || start > contigLength )
            Assert.assertNull(loc, "GenomeLoc should be null if the start/stops are not meaningful");
        else {
            Assert.assertNotNull(loc);
            Assert.assertEquals(loc.getContig(), contig);
            Assert.assertEquals(loc.getStart(), Math.max(start, 1));
            Assert.assertEquals(loc.getStop(), Math.min(stop, contigLength));
        }
    }

    @DataProvider(name = "GenomeLocPadding")
    public Object[][] makeGenomeLocPadding() {
        final List<Object[]> tests = new LinkedList<>();

        final int contigLength = header.getSequence(0).getSequenceLength();
        for ( int pad = 0; pad < contigLength + 1; pad++) {
            for ( int start = 1; start < contigLength; start++ ) {
                for ( int stop = start; stop < contigLength; stop++ ) {
                    tests.add(new Object[]{ genomeLocParser.createGenomeLoc("1", start, stop), pad});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test( dataProvider = "GenomeLocPadding")
    public void testGenomeLocPadding(final GenomeLoc input, final int pad) {
        final int contigLength = header.getSequence(0).getSequenceLength();
        final GenomeLoc padded = genomeLocParser.createPaddedGenomeLoc(input, pad);

        Assert.assertNotNull(padded);
        Assert.assertEquals(padded.getContig(), input.getContig());
        Assert.assertEquals(padded.getStart(), Math.max(input.getStart() - pad, 1));
        Assert.assertEquals(padded.getStop(), Math.min(input.getStop() + pad, contigLength));
    }

    @Test
    public void testQueryAllHG38Intervals() {
        SAMSequenceDictionary sd;
        final File testFile = new File (publicTestDir, "org/broadinstitute/hellbender/engine/Homo_sapiens_assembly38.headerOnly.vcf.gz");

        try (VCFFileReader vcfReader = new VCFFileReader(testFile, false)) {
            sd = vcfReader.getFileHeader().getSequenceDictionary();
        }

        // Test that we can use any contig from hg38 as a query against a VCF with an hg38 sequence dictionary, in any
        // query format, without ambiguity.
        final GenomeLocParser localGenomeLocParser = new GenomeLocParser(sd);
        sd.getSequences().stream().forEach(
                hg38Contig -> {
                    assertValidUniqueInterval(
                            localGenomeLocParser,
                            hg38Contig.getSequenceName(),
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength()));
                    assertValidUniqueInterval(
                            localGenomeLocParser,
                            hg38Contig.getSequenceName() + ":1",
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, 1));
                    assertValidUniqueInterval(
                            localGenomeLocParser,
                            hg38Contig.getSequenceName() + ":1+",
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, hg38Contig.getSequenceLength()));
                    assertValidUniqueInterval(
                            localGenomeLocParser,
                            hg38Contig.getSequenceName() + ":1-1",
                            new SimpleInterval(hg38Contig.getSequenceName(), 1, 1));
                }
        );
    }

    private void assertValidUniqueInterval(
            final GenomeLocParser localGenomeLocParser,
            final String queryString,
            final SimpleInterval expectedInterval) {
        final SimpleInterval actualInterval = new SimpleInterval(localGenomeLocParser.parseGenomeLoc(queryString));
        Assert.assertEquals(actualInterval, expectedInterval);
    }

}
