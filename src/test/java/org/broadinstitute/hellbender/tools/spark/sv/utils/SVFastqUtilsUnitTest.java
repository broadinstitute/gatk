package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils.singletonIterator;

/**
 * Unit tests for components in {@link SVFastqUtils}.
 */
public class SVFastqUtilsUnitTest extends GATKBaseTest {


    private static final String[] TEST_VALID_CHROMOSOMES = {
            "chr1", "chr2", "X", "Y", "MT"
    };

    private static final String[] TEST_MAPPED_VALID_STRINGS = {
            "chr1,100,+,100M,60,2,35", // single interval
            "chr2,101,-,50M,60,.,10", // another single interval missing NM.
            "10,51,+,10S45M,30;20,1213,-,10M45H,.,10,30", // two intervals
            "X,1231,-,10S45M31S,20,.;Y,3112,+,10M76H,255;MT,1,-,55H31M,60,0,112" // three intervals.
    };

    private static final String[] TEST_INVALID_STRINGS = {
            null,
            "chr1", // partials
            "chr1,100",
            "chr1,100,+",
            "chr1,100,+,10M",
            "chr1,-1,-,10M,60", // bad values.
            "chr1,100,blah,10M,60",
            "chr1,100,,10M,60",
            "X,100,-,blah,60",
            "X,100,-,10M,-1",
            "X,100,-,10M,60,.,-1.0",
            "Y,100,-,10H10M,60,30,30;Y,200,10M10H,60,20,.", // bad cigar as first must not contain hard-clips.
    };


    @Test(dataProvider = "mappingFromStringData", dependsOnMethods = "testMappingFromString", groups = "sv")
    public void testMappingToString(final String str, final String[] contig, final int[] start,
                                      final boolean[] forward, final Cigar[] cigarAlongRef,
                                      final int[] mappingQual, final int[] mismatches,
                                      final int[] alignmentScores) {
        final SVFastqUtils.Mapping mapping = new SVFastqUtils.Mapping(str);
        final String backToString = mapping.toString();
        Assert.assertNotNull(backToString);
        if (!backToString.equals(str)) { // this can be the case because of the '.' ambiguitities.
            final SVFastqUtils.Mapping backToMapping = new SVFastqUtils.Mapping(backToString);
            assertMappingIsAsExpected(backToMapping, contig, start, forward, cigarAlongRef, mappingQual, mismatches, alignmentScores);
        }
    }

    @Test(dataProvider = "mappingFromStringData", groups = "sv")
    public void testMappingFromString(final String str, final String[] contig, final int[] start,
                                      final boolean[] forward, final Cigar[] cigarAlongRef,
                                      final int[] mappingQual, final int[] mismatches,
                                      final int[] alignmentScores) {
        final SVFastqUtils.Mapping mapping = new SVFastqUtils.Mapping(str);
        assertMappingIsAsExpected(mapping, contig, start, forward, cigarAlongRef, mappingQual, mismatches, alignmentScores);
    }

    @Test(dataProvider = "mappingFromBadStringData", expectedExceptions = IllegalArgumentException.class, groups = "sv")
    public void testMappingFromBadString(final String str) {
        new SVFastqUtils.Mapping(str);
    }

    @Test(groups = "sv")
    public void testUnmappedFromString() {
        assertUnmapped(new SVFastqUtils.Mapping("*"));
    }

    @Test(expectedExceptions = GATKException.class, groups = "sv")
    public void testNonPrimaryAlignmentReadWithoutSA() {
        final SAMFileHeader header;
        final SAMRecord samRecord;
        final GATKRead read;
        try {
            header = new SAMFileHeader(
                    new SAMSequenceDictionary(Collections.singletonList(
                            new SAMSequenceRecord("seq1", 10_000))));
            samRecord = new SAMRecord(header);
            samRecord.setReadName("test");
            samRecord.setReadUnmappedFlag(false);
            samRecord.setReferenceName("seq1");
            samRecord.setAlignmentStart(100);
            samRecord.setCigar(new Cigar(Collections.singletonList(new CigarElement(100, CigarOperator.M))));
            samRecord.setSupplementaryAlignmentFlag(true);
            read = new SAMRecordToGATKReadAdapter(samRecord);
        } catch (final RuntimeException ex) {
            Assert.fail("exception thrown during preparation");
            throw new GATKException.ShouldNeverReachHereException(ex);
        }

        new SVFastqUtils.Mapping(read);
    }

    @Test(expectedExceptions = GATKException.class, groups = "sv")
    public void testNonPrimaryAlignmentWithSAWhereFirstElementHasHardClips() {
        final SAMFileHeader header;
        final SAMRecord samRecord;
        final GATKRead read;
        try {
            header = new SAMFileHeader(
                    new SAMSequenceDictionary(Collections.singletonList(
                            new SAMSequenceRecord("seq1", 10_000))));
            samRecord = new SAMRecord(header);
            samRecord.setReadName("test");
            samRecord.setReadUnmappedFlag(false);
            samRecord.setReferenceName("seq1");
            samRecord.setAlignmentStart(100);
            samRecord.setCigar(new Cigar(Collections.singletonList(new CigarElement(100, CigarOperator.M))));
            samRecord.setSupplementaryAlignmentFlag(true);
            samRecord.setAttribute(SAMTag.SA.name(), "seq1,200,-,50H100M50H,60;seq1,300,+,200M,30;");
            read = new SAMRecordToGATKReadAdapter(samRecord);
        } catch (final RuntimeException ex) {
            Assert.fail("exception thrown during preparation");
            throw new GATKException.ShouldNeverReachHereException(ex);
        }

        new SVFastqUtils.Mapping(read);
    }

    @Test(expectedExceptions = GATKException.class, groups = "sv")
    public void testPrimaryAlignmentWithHardClips() {
        final SAMFileHeader header;
        final SAMRecord samRecord;
        final GATKRead read;
        try {
            header = new SAMFileHeader(
                    new SAMSequenceDictionary(Collections.singletonList(
                            new SAMSequenceRecord("seq1", 10_000))));
            samRecord = new SAMRecord(header);
            samRecord.setReadName("test");
            samRecord.setReadUnmappedFlag(false);
            samRecord.setReferenceName("seq1");
            samRecord.setAlignmentStart(100);
            samRecord.setCigar(new Cigar(Arrays.asList(new CigarElement(100, CigarOperator.H), new CigarElement(100, CigarOperator.M))));
            samRecord.setSupplementaryAlignmentFlag(false);
            read = new SAMRecordToGATKReadAdapter(samRecord);

        } catch (final RuntimeException ex) {
            Assert.fail("exception thrown during preparation");
            throw new GATKException.ShouldNeverReachHereException(ex);
        }

        new SVFastqUtils.Mapping(read);
    }

    @Test(groups = "sv")
    public void testNonPrimaryAlignmentReadWithSA() {
        final SAMFileHeader header;
        final SAMRecord samRecord;
        header = new SAMFileHeader(
                new SAMSequenceDictionary(Collections.singletonList(
                        new SAMSequenceRecord("seq1", 10_000))));
        samRecord = new SAMRecord(header);
        samRecord.setReadName("test");
        samRecord.setReadUnmappedFlag(false);
        samRecord.setReferenceName("seq1");
        samRecord.setAlignmentStart(100);
        samRecord.setCigar(new Cigar(Arrays.asList(new CigarElement(100, CigarOperator.H), new CigarElement(100, CigarOperator.M))));
        samRecord.setSupplementaryAlignmentFlag(true);
        samRecord.setAttribute(SAMTag.SA.name(), "seq1,200,+,100M100S,60,.,20;");

        final GATKRead read = new SAMRecordToGATKReadAdapter(samRecord);

        final SVFastqUtils.Mapping mapping = new SVFastqUtils.Mapping(read);
        assertMappingIsAsExpected(mapping,
                new String[] { "seq1", "seq1" }, new int[] { 200, 100},
                new boolean[] {true, true},
                new Cigar[] { TextCigarCodec.decode("100M100S"), TextCigarCodec.decode("100H100M")},
                new int[] { 60, 0}, new int[] {AlignmentInterval.NO_NM, AlignmentInterval.NO_NM},
                new int[] {20, AlignmentInterval.NO_AS});
    }

    @Test(groups = "sv")
    public void testUnmappedFromGATKRead() {
        final SAMFileHeader header = new SAMFileHeader(
                new SAMSequenceDictionary(Collections.singletonList(
                        new SAMSequenceRecord("seq1", 10_000))));
        final SAMRecord samRecord = new SAMRecord(header);
        samRecord.setReadName("test");
        samRecord.setReadUnmappedFlag(true);
        final GATKRead read = new SAMRecordToGATKReadAdapter(samRecord);
        assertUnmapped(new SVFastqUtils.Mapping(read));
        Assert.assertEquals(SVFastqUtils.Mapping.toString(read), "*");
    }

    @Test(groups = "sv")
    public void testUnmappedFromFastqRead() {
        final byte[] bases = new byte[100];
        final byte[] quals = new byte[100];
        Arrays.fill(bases, (byte)'A');
        Arrays.fill(quals, (byte)'K');
        final SVFastqUtils.FastqRead fastqRead  = new SVFastqUtils.FastqRead("@test mapping=*", bases, quals);
        final SVFastqUtils.Mapping mapping = fastqRead.getMapping();
        assertUnmapped(mapping);
    }


    private void assertUnmapped(SVFastqUtils.Mapping mapping) {
        Assert.assertFalse(mapping.isMapped());
        Assert.assertNull(mapping.getPrimaryInterval());
        Assert.assertNull(mapping.getContig());
        Assert.assertEquals(mapping.getCigar(), new Cigar());
        Assert.assertEquals(mapping.getStart(), SAMRecord.NO_ALIGNMENT_START);
        Assert.assertEquals(mapping.getEnd(), SAMRecord.NO_ALIGNMENT_START);
        Assert.assertTrue(mapping.getAllIntervals().isEmpty());
        Assert.assertTrue(mapping.getSupplementaryIntervals().isEmpty());
        Assert.assertEquals(mapping.toString(), "*");
    }

    @Test(dataProvider = "mappingFromGATKReadData", groups = "sv")
    public void testMappingFromGATKRead(final GATKRead read, final String[] contig, final int[] start,
                                        final boolean[] forward, final Cigar[] cigarAlongRef,
                                        final int[] mappingQual, final int[] mismatches,
                                        final int[] alignmentScores) {
        final SVFastqUtils.Mapping mapping = new SVFastqUtils.Mapping(read);
        assertMappingIsAsExpected(mapping, contig, start, forward, cigarAlongRef, mappingQual, mismatches, alignmentScores);
    }

    @Test(dataProvider = "mappingFromFastqReadData", groups = "sv")
    public void testMappingFromFastqRead(final SVFastqUtils.FastqRead read, final String[] contig, final int[] start,
                                        final boolean[] forward, final Cigar[] cigarAlongRef,
                                        final int[] mappingQual, final int[] mismatches,
                                        final int[] alignmentScores) {
        final SVFastqUtils.Mapping mapping = read.getMapping();
        assertMappingIsAsExpected(mapping, contig, start, forward, cigarAlongRef, mappingQual, mismatches, alignmentScores);
    }

    @Test(dataProvider = "mappingFromGATKReadData", dependsOnMethods = "testMappingFromGATKRead", groups = "sv")
    public void testMappingToStringRead(final GATKRead read, final String[] contig, final int[] start,
                                        final boolean[] forward, final Cigar[] cigarAlongRef,
                                        final int[] mappingQual, final int[] mismatches,
                                        final int[] alignmentScores) {
        final String asString = SVFastqUtils.Mapping.toString(read);
        final SVFastqUtils.Mapping mapping = new SVFastqUtils.Mapping(asString);
        assertMappingIsAsExpected(mapping, contig, start, forward, cigarAlongRef, mappingQual, mismatches, alignmentScores);
    }

    private void assertMappingIsAsExpected(final SVFastqUtils.Mapping mapping, final String[] contig, final int[] start,
                                      final boolean[] forward, final Cigar[] cigarAlongRef,
                                      final int[] mappingQual, final int[] mismatches,
                                      final int[] alignmentScores) {
        Assert.assertEquals(mapping.getAllIntervals().size(), contig.length);
        Assert.assertEquals(mapping.isMapped(), contig.length > 0);
        if (!mapping.isMapped()) {
            Assert.assertNull(mapping.getPrimaryInterval());
            Assert.assertTrue(mapping.getSupplementaryIntervals().isEmpty());
            Assert.assertEquals(mapping.getCigar(), new Cigar());
            Assert.assertNull(mapping.getContig());
            Assert.assertEquals(mapping.getStart(), SAMRecord.NO_ALIGNMENT_START);
            Assert.assertEquals(mapping.getEnd(), SAMRecord.NO_ALIGNMENT_START);
        } else {
            Assert.assertEquals(mapping.getSupplementaryIntervals().size(), contig.length - 1);
            for (int i = 0; i < contig.length; i++) {
                final AlignmentInterval ai = mapping.getAllIntervals().get(i);
                Assert.assertNotNull(ai);
                Assert.assertEquals(ai.referenceSpan.getContig(), contig[i]);
                Assert.assertEquals(ai.referenceSpan.getStart(), start[i]);
                Assert.assertEquals(ai.forwardStrand, forward[i]);
                Assert.assertEquals(ai.cigarAlongReference(), cigarAlongRef[i]);
                Assert.assertEquals(ai.mapQual, mappingQual[i]);
                Assert.assertEquals(ai.mismatches, mismatches[i]);
                Assert.assertEquals(ai.alnScore, alignmentScores[i]);
                if (i == 0) {
                    Assert.assertSame(ai, mapping.getPrimaryInterval());
                    Assert.assertEquals(mapping.getCigar(), ai.cigarAlongReference());
                    Assert.assertEquals(mapping.getContig(), ai.referenceSpan.getContig());
                    Assert.assertEquals(mapping.getStart(), ai.referenceSpan.getStart());
                    Assert.assertEquals(mapping.getEnd(), ai.referenceSpan.getEnd());
                    Assert.assertEquals(mapping.isForwardStrand(), ai.forwardStrand);
                } else {
                    Assert.assertSame(ai, mapping.getSupplementaryIntervals().get(i - 1));
                }
            }
        }

    }

    @Test
    public void testSimpleRoundTrip() {
        final byte[] calls = "ACGTTGCAACGTTGCAACGT".getBytes();
        final byte[] quals = new byte[calls.length];
        Arrays.fill(quals, (byte) 20);
        final SVFastqUtils.FastqRead read = new SVFastqUtils.FastqRead("@readName", calls, quals);
        final ByteArrayOutputStream os = new ByteArrayOutputStream();
        try {
            SVFastqUtils.writeFastqStream(os, singletonIterator(read));
        }
        catch ( final IOException ioe ) {
            throw new GATKException("failed to write FASTQ file to memory", ioe);
        }
        try {
            final BufferedReader reader =
                    new BufferedReader(new InputStreamReader(new ByteArrayInputStream(os.toByteArray())));
            final List<SVFastqUtils.FastqRead> reads = SVFastqUtils.readFastqStream(reader, "from memory");
            Assert.assertEquals(reads.size(), 1);
            final SVFastqUtils.FastqRead read2 = reads.get(0);
            Assert.assertEquals(read.getHeader(), read2.getHeader());
            Assert.assertTrue(Arrays.equals(read.getBases(), read2.getBases()));
            Assert.assertTrue(Arrays.equals(read.getQuals(), read2.getQuals()));
        }
        catch ( final IOException ioe ) {
            throw new GATKException("failed to read FASTQ from memory", ioe);
        }
    }


    @DataProvider(name="mappingFromBadStringData")
    public static Object[][] mappingFromBadStringData() {
        return Arrays.stream(TEST_INVALID_STRINGS)
                .map(s -> new Object[] { s })
                .toArray(Object[][]::new);
    }

    @DataProvider(name = "mappingFromFastqReadData")
    public Object[][] mappingFromFastqReadData() {
        final RandomDNA rdnDNA = new RandomDNA(113);
        final byte[] quals = new byte[100];
        Arrays.fill(quals, (byte) 'K');
        return mappingFromData(str -> new SVFastqUtils.FastqRead("@test   mapping=" + str, rdnDNA.nextBases(100), quals));
    }

    @DataProvider(name = "mappingFromStringData")
    public Object[][] mappingFromStringData() {
        return mappingFromData(Function.identity());
    }

    @DataProvider(name = "mappingFromGATKReadData")
    public Object[][] mappingFromGATKReadData() {
        final SAMFileHeader header = new SAMFileHeader(
                new SAMSequenceDictionary(Arrays.stream(TEST_VALID_CHROMOSOMES)
                    .map(s -> new SAMSequenceRecord(s, 10_000))
                    .collect(Collectors.toList()))
        );
        final RandomDNA rdnDNA = new RandomDNA(101);
        return mappingFromData(str -> {
            final SAMRecord samRecord = new SAMRecord(header);
            samRecord.setReadName("test");
            final List<AlignmentInterval> alignmentIntervals = Arrays.stream(str.split(";"))
                    .filter(s -> !s.isEmpty())
                    .map(AlignmentInterval::new)
                    .collect(Collectors.toList());
            if (alignmentIntervals.isEmpty()) {
                samRecord.setReadUnmappedFlag(true);
                samRecord.setReadBases(rdnDNA.nextBases(100));
            } else {
                final AlignmentInterval primaryInterval = alignmentIntervals.get(0);
                samRecord.setReadNegativeStrandFlag(!primaryInterval.forwardStrand);
                samRecord.setReadUnmappedFlag(false);
                samRecord.setReferenceName(primaryInterval.referenceSpan.getContig());
                samRecord.setAlignmentStart(primaryInterval.referenceSpan.getStart());
                samRecord.setCigar(primaryInterval.cigarAlongReference());
                samRecord.setReadBases(rdnDNA.nextBases(primaryInterval.cigarAlong5to3DirectionOfContig.getReadLength()));
                samRecord.setMappingQuality(primaryInterval.mapQual);
                if (primaryInterval.mismatches != AlignmentInterval.NO_NM) {
                    samRecord.setAttribute(SAMTag.NM.name(), primaryInterval.mismatches);
                }
                if (primaryInterval.alnScore != AlignmentInterval.NO_AS) {
                    samRecord.setAttribute(SAMTag.AS.name(), primaryInterval.alnScore);
                }
                if (alignmentIntervals.size() > 1) {
                    samRecord.setAttribute(SAMTag.SA.name(),
                            alignmentIntervals.subList(1, alignmentIntervals.size()).stream()
                                .map(AlignmentInterval::toSATagString).collect(Collectors.joining(";","", ";")));
                }
            }
            return new SAMRecordToGATKReadAdapter(samRecord);
        });
    }

    private Object[][] mappingFromData(final Function<String, ?> testStringTrasform) {

        // code below takes care of breaking the tests above into its parts to feed into the
        // the test methods...
        // if you want to add additional test just add the corresponding SA styled string above.
        final String[][] contigs = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(s -> s.split(",")[0]).toArray(String[]::new))
                .toArray(String[][]::new);

        final int[][] starts = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(x -> x.split(",")[1]).mapToInt(Integer::parseInt).toArray())
                .toArray(int[][]::new);

        final boolean[][] forwardStrand = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(x -> x.split(",")[2]).map(x -> {
                    if (!x.equals("+") && !x.equals("-")) {
                        throw new IllegalArgumentException();
                    }
                    return x.equals("+");}).toArray(Boolean[]::new))
                .map(B -> {
                    final boolean[] result = new boolean[B.length];
                    for (int i = 0; i < result.length; i++) {
                        result[i] = B[i];
                    }
                    return result;
                })
                .toArray(boolean[][]::new);

        final Cigar[][] cigars = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(s -> TextCigarCodec.decode(s.split(",")[3])).toArray(Cigar[]::new))
                .toArray(Cigar[][]::new);

        final int[][] mappingQuals = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(x -> x.split(",")[4]).mapToInt(s -> s.equals(".") ? SAMRecord.UNKNOWN_MAPPING_QUALITY : Integer.parseInt(s)).toArray())
                .toArray(int[][]::new);

        final int[][] mismatches = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(x -> {
                    final String[] split = x.split(",");
                    return (split.length < 6) ? "." : split[5];
                }).mapToInt(s -> s.equals(".") ? AlignmentInterval.NO_NM : Integer.parseInt(s)).toArray())
                .toArray(int[][]::new);

        final int[][] alignmentScores = Arrays.stream(TEST_MAPPED_VALID_STRINGS)
                .map(s -> Arrays.stream(s.split(";")))
                .map(parts -> parts.map(x ->  {
                    final String[] split = x.split(",");
                    return (split.length < 7) ? "." : split[6];
                }).mapToInt(s -> s.equals(".") ? AlignmentInterval.NO_NM : Integer.parseInt(s)).toArray())
                .toArray(int[][]::new);

        return IntStream.range(0, TEST_MAPPED_VALID_STRINGS.length)
                .mapToObj(i -> new Object[] { testStringTrasform.apply(TEST_MAPPED_VALID_STRINGS[i]), contigs[i], starts[i], forwardStrand[i],
                  cigars[i], mappingQuals[i], mismatches[i], alignmentScores[i]})
                .toArray(Object[][]::new);

    }
}
