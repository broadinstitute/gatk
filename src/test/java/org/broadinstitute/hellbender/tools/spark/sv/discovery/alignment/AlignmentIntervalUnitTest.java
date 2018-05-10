package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.utils.Strand;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignmentUtils;
import org.broadinstitute.hellbender.utils.read.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

public class AlignmentIntervalUnitTest extends GATKBaseTest {
    private static final String dummyRefName = "1";
    private static final List<String> refNames = Collections.singletonList(dummyRefName);

    @DataProvider(name = "TestDataForAIOverlaps")
    Object[][] testDataForAIOverlaps() {
        final List<Object[]> data = new ArrayList<>(20);

        AlignmentInterval ar1 = new AlignmentInterval(new SimpleInterval("1",1,5), 1,5, TextCigarCodec.decode("5M5H"),true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval ar2 = new AlignmentInterval(new SimpleInterval("1",11,16), 5,10, TextCigarCodec.decode("4S6M"),true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        data.add(new Object[]{ar1, ar2, 1, 0});

        ar1 = new AlignmentInterval(new SimpleInterval("1",1,5), 1,5, TextCigarCodec.decode("5M5H"),true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("1",11,15), 6,10, TextCigarCodec.decode("5S5M"),true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 0, 0});

        // overlaps on ref only
        ar1 = new AlignmentInterval(new SimpleInterval("chr1",4938770,4939497), 1,728, TextCigarCodec.decode("728M61S"),true, 60, 0, 728, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr1",4939439,4939498), 730,789, TextCigarCodec.decode("729H60M"),true, 60, 2, 50, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 0, 59});

        ar1 = new AlignmentInterval(new SimpleInterval("chr1",9170350,9171390), 1,1041, TextCigarCodec.decode("1041M1298H"),false, 60, 4, 1021, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr1",9169370,9170505), 1204,2339, TextCigarCodec.decode("1203S1136M"),false, 60, 22, 1026, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 0, 505-350+1});

        // overlaps on read only
        ar1 = new AlignmentInterval(new SimpleInterval("chr1",933803,934119), 1,317, TextCigarCodec.decode("317M302H"),true, 60, 7, 282, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr1",934806,935261), 164,619, TextCigarCodec.decode("163S456M"),true, 60, 8, 416, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 317-164+1, 0});

        ar1 = new AlignmentInterval(new SimpleInterval("chr1",964783,965113), 1,331, TextCigarCodec.decode("331M1028H"),false, 60, 2, 321, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr1",963604,964692), 270,1359, TextCigarCodec.decode("269S69M1I1020M"),false, 60, 9, 1032, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 331-270+1, 0});

        // overlaps on read & ref
        ar1 = new AlignmentInterval(new SimpleInterval("chr1",66659809,66660176), 1,354, TextCigarCodec.decode("124M10D106M3D16M2I75M3D31M241H"),true, 60, 35, 185, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr1",66659958,66660262), 301,595, TextCigarCodec.decode("300S179M10D116M"),true, 60, 24, 199, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 54, 66660176-66659958+1});

        ar1 = new AlignmentInterval(new SimpleInterval("chr1",156328046,156328757), 1,712, TextCigarCodec.decode("712M444S"),false, 60, 2, 702, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr1",156327744,156328331), 588,1156, TextCigarCodec.decode("587H127M15I131M34D296M"),false, 60, 68, 378, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 712-588+1, 331-46+1});

        // overlap with strand switch
        ar1 = new AlignmentInterval(new SimpleInterval("chr6",148696358,148697176), 1,815, TextCigarCodec.decode("725M4D90M472S"),true, 60, 10, 765, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr6",4090017,4090739), 567,1287, TextCigarCodec.decode("566H80M2D641M"),false, 60, 7, 678, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 815-567+1, 0});

        ar1 = new AlignmentInterval(new SimpleInterval("chr5",180678871,180679093), 1,223, TextCigarCodec.decode("223M44S"),false, 60, 0, 223, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr5",180678907,180678954), 220,267, TextCigarCodec.decode("219H48M"),true, 60, 0, 48, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 4, 48});

        // different chr
        ar1 = new AlignmentInterval(new SimpleInterval("chr1",9170350,9171390), 1,1041, TextCigarCodec.decode("1041M1298H"),false, 60, 4, 1021, ContigAlignmentsModifier.AlnModType.NONE);
        ar2 = new AlignmentInterval(new SimpleInterval("chr2",9169370,9170505), 1204,2339, TextCigarCodec.decode("1203S1136M"),false, 60, 22, 1026, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{ar1, ar2, 0, 0});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "TestDataForAIOverlaps", groups = "sv")
    public void testAlignmentIntervalOverlap(final AlignmentInterval ar1, final AlignmentInterval ar2,
                                             final int expectedOverlapOnRead, final int expectedOverlapOnRef) {

        Assert.assertEquals(AlignmentInterval.overlapOnContig(ar1, ar2), expectedOverlapOnRead);
        Assert.assertEquals(AlignmentInterval.overlapOnRefSpan(ar1, ar2), expectedOverlapOnRef);
    }

    /**
     * These alignment records are supposed to be associated with the 4 possible types of evidence we could see for an inversion,
     *   where the chr1:101-200 bases are inverted, namely
     * INV55, where lower  contig coordinate is associated with a forward  strand lower reference coordinate, and
     *              higher contig coordinate is associated with a negative strand higher reference/contig coordinate
     * INV55, where lower  contig coordinate is associated with a forward  strand higher reference coordinate, and
     *              higher contig coordinate is associated with a negative strand lower reference/contig coordinate
     * INV33, where lower  contig coordinate is associated with a negative strand lower reference coordinate, and
     *              higher contig coordinate is associated with a forward  strand higher reference/contig coordinate
     * INV33, where lower  contig coordinate is associated with a forward  strand higher reference coordinate, and
     *              higher contig coordinate is associated with a negative strand lower reference/contig coordinate
     * Finally, one must be aware of the fact that BWA always outputs CIGAR with a '+'-strand representation,
     *   therefore we must use such in constructing the BwaMemAlignment's* @return objects stored in each array
     * @return an array of arrays, each composed of
     * [0] {@link BwaMemAlignment} object,
     * [1] expected reference interval,
     * [2] expected cigar,
     * [3] expected strandedness,
     * [4] expected start in assembled contig, 1-based, inclusive
     * [5] expected end in assembled contig, 1-based, inclusive
     * [6] expected contig length,
     * [7] expected {@link AlignmentInterval} object (generated manually with all fields explicitly spell out and given to
     *                                      {@link AlignmentInterval#AlignmentInterval(SimpleInterval, int, int, Cigar, boolean, int, int, int, ContigAlignmentsModifier.AlnModType)}
     *                                      intended to be used for testing concordance between the two constructors)
     */
    @DataProvider(name = "AlignmentIntervalCtorTestForSimpleInversion")
    Object[][] createInputsAndExpectedResults_BwaMemAlignmentConstruction() {

        final int[] alignmentStartsOnRef_0Based = {96, 196, 195, 95, 101, 201, 101, 201};
        final int[] alignmentStartsOnTig_0BasedInclusive = {0, 4, 0, 5, 0, 6, 0, 7};
        final int[] alignmentEndsOnTig_0BasedExclusive = {4, 8, 5, 10, 6, 12, 7, 14};
        final int[] seqLen = {8, 8, 10, 10, 12, 12, 14, 14};
        final int[] mapQualForBwaMemAlgn = {-1, 0, 10, 20, 30, 40, 50, 60};
        final boolean[] strandedness = {true, false, true, false, false, true, false, true};
        final String[] cigarStrings = {"4M4S", "4M4H", "5M5S", "5M5H", "6S6M", "6H6M", "7S7M", "7H7M"}; // each different number represent a different contig's pair of chimeric alignments
        final Cigar[] cigars = Arrays.stream(cigarStrings).map(TextCigarCodec::decode).toArray(Cigar[]::new);


        final Object[][] data = new Object[cigars.length][];
        for(int i=0; i<cigars.length; ++i) {
            int samFlag = 0;
            if ( !strandedness[i] ) samFlag = SAMFlag.READ_REVERSE_STRAND.intValue();
            if ( cigarStrings[i].indexOf('H') != -1 ) samFlag |= SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue();
            final BwaMemAlignment bwaMemAlignment = new BwaMemAlignment(samFlag,
                    0, alignmentStartsOnRef_0Based[i], alignmentStartsOnRef_0Based[i]+cigars[i].getReferenceLength(),
                    strandedness[i] ? alignmentStartsOnTig_0BasedInclusive[i] : seqLen[i]-alignmentEndsOnTig_0BasedExclusive[i],
                    strandedness[i] ? alignmentEndsOnTig_0BasedExclusive[i] : seqLen[i]-alignmentStartsOnTig_0BasedInclusive[i],
                    mapQualForBwaMemAlgn[i], 0, 1, 1, cigarStrings[i],
                    null, null, 0, Integer.MIN_VALUE, Integer.MAX_VALUE);
            final SimpleInterval referenceInterval = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[i]+1, bwaMemAlignment.getRefEnd());
            final AlignmentInterval alignmentInterval = new AlignmentInterval(referenceInterval, alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i],
                    strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]),
                    strandedness[i], Math.max(SAMRecord.NO_MAPPING_QUALITY, bwaMemAlignment.getMapQual()), bwaMemAlignment.getNMismatches(), bwaMemAlignment.getAlignerScore(), ContigAlignmentsModifier.AlnModType.NONE);
            data[i] = new Object[]{bwaMemAlignment, referenceInterval, strandedness[i] ? cigars[i] : CigarUtils.invertCigar(cigars[i]),
                    strandedness[i], alignmentStartsOnTig_0BasedInclusive[i]+1, alignmentEndsOnTig_0BasedExclusive[i], seqLen[i], mapQualForBwaMemAlgn[i], alignmentInterval};
        }
        return data;
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testCigarAlongTheReference(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                           final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                           final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final AlignmentInterval alignmentInterval = new AlignmentInterval(bwaMemAlignment, refNames, expectedContigLength);
        final Cigar fiveToThree = alignmentInterval.cigarAlong5to3DirectionOfContig;
        final Cigar alongRef = alignmentInterval.cigarAlongReference();
        if (alignmentInterval.forwardStrand) {
            Assert.assertEquals(new ArrayList<>(alongRef.getCigarElements()), new ArrayList<>(fiveToThree.getCigarElements()));
        } else {
            Assert.assertEquals(new ArrayList<>(CigarUtils.invertCigar(alongRef).getCigarElements()), new ArrayList<>(fiveToThree.getCigarElements()));
        }
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testToSATagString(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                  final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                  final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final AlignmentInterval alignmentInterval = new AlignmentInterval(bwaMemAlignment, refNames, expectedContigLength);
        final String saTagString = alignmentInterval.toSATagString();
        final AlignmentInterval fromString = new AlignmentInterval(saTagString);
        Assert.assertEquals(fromString, alignmentInterval);
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testAppendToBuilder(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                    final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                    final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final StringBuilder builder = new StringBuilder();
        builder.append("BEFORE->");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(bwaMemAlignment, refNames, expectedContigLength);
        alignmentInterval.appendSATagString(builder);
        builder.append("<-AFTER");
        final String fullString = builder.toString();
        final String aiString = fullString.replace("BEFORE->","").replace("<-AFTER","");
        Assert.assertEquals(aiString, alignmentInterval.toSATagString());
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromBwaMemAlignment(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                                    final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                                    final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final AlignmentInterval alignmentInterval = new AlignmentInterval(bwaMemAlignment, refNames, expectedContigLength);
        Assert.assertEquals(alignmentInterval.referenceSpan, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromSAMRecord(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                              final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                              final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final SAMRecord samRecord = BwaMemAlignmentUtils.applyAlignment("whatever", TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(expectedContigLength, (byte)'A'), null, null, bwaMemAlignment, refNames, hg19Header, false, false);
        final AlignmentInterval alignmentInterval = new AlignmentInterval(samRecord);
        Assert.assertEquals(alignmentInterval.referenceSpan, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromStr(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                        final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                        final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final SAMRecord samRecord = BwaMemAlignmentUtils.applyAlignment("whatever", TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(expectedContigLength, (byte)'A'), null, null, bwaMemAlignment, refNames, hg19Header, false, false);
        final StringBuilder strBuilder = new StringBuilder(String.join(",", samRecord.getContig(),
                "" + samRecord.getStart(), samRecord.getReadNegativeStrandFlag() ? "-" : "+", samRecord.getCigarString(), "" + samRecord.getMappingQuality()));
        if (samRecord.getAttribute(SAMTag.NM.name()) != null || samRecord.getAttribute(SAMTag.AS.name()) != null) {
            strBuilder.append("," + samRecord.getIntegerAttribute(SAMTag.NM.name()));
            if (samRecord.getAttribute(SAMTag.AS.name()) != null) {
                strBuilder.append("," + samRecord.getIntegerAttribute(SAMTag.AS.name()));
            }
        }
        AlignmentInterval alignmentInterval = new AlignmentInterval(strBuilder.toString());
        Assert.assertEquals(alignmentInterval.referenceSpan, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive, bwaMemAlignment.getCigar());
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testConstructionFromGATKRead(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                             final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                             final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final SAMRecord samRecord = BwaMemAlignmentUtils.applyAlignment("whatever", TestUtilsForAssemblyBasedSVDiscovery.makeDummySequence(expectedContigLength, (byte)'A'), null, null, bwaMemAlignment, refNames, hg19Header, false, false);
        final GATKRead read = new SAMRecordToGATKReadAdapter(samRecord);
        final AlignmentInterval alignmentInterval = new AlignmentInterval(read);
        Assert.assertEquals(alignmentInterval.referenceSpan, expectedReferenceInterval);
        Assert.assertEquals(alignmentInterval.cigarAlong5to3DirectionOfContig, expectedCigar);
        Assert.assertEquals(alignmentInterval.forwardStrand, expectedIsPositiveStrand);
        Assert.assertEquals(alignmentInterval.startInAssembledContig, expectedStartOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.endInAssembledContig, expectedEndOnContig_1BasedInclusive);
        Assert.assertEquals(alignmentInterval.mapQual, Math.max(SAMRecord.NO_MAPPING_QUALITY,expectedMapQualInBwaMemAlignment));
        Assert.assertEquals(alignmentInterval, expectedAlignmentInterval);
    }

    @Test(dataProvider = "AlignmentIntervalCtorTestForSimpleInversion", groups = "sv")
    public void testToSAMRecord(final BwaMemAlignment bwaMemAlignment, final SimpleInterval expectedReferenceInterval, final Cigar expectedCigar,
                                final boolean expectedIsPositiveStrand, final int expectedStartOnContig_1BasedInclusive, final int expectedEndOnContig_1BasedInclusive,
                                final int expectedContigLength, final int expectedMapQualInBwaMemAlignment, final AlignmentInterval expectedAlignmentInterval) {

        final byte[] randomContigBases = new RandomDNA(13).nextBases(expectedContigLength);
        final SAMRecord samRecord = BwaMemAlignmentUtils.applyAlignment("whatever", randomContigBases, null, null, bwaMemAlignment, refNames, hg19Header, false, false);
        final AlignmentInterval alignmentInterval = new AlignmentInterval(samRecord);
        final SAMRecord backSamRecord = alignmentInterval.toSAMRecord(samRecord.getHeader(), samRecord.getReadName(), randomContigBases, samRecord.getCigar().containsOperator(CigarOperator.H), samRecord.getFlags() , samRecord.getAttributes());
        Assert.assertEquals(backSamRecord.getReadName(), samRecord.getReadName());
        Assert.assertEquals(backSamRecord.getFlags(), samRecord.getFlags());
        Assert.assertEquals(backSamRecord.getAttributes().stream().collect(Collectors.toMap(x -> x.tag, x -> x.value)), samRecord.getAttributes().stream().collect(Collectors.toMap(x -> x.tag, x -> x.value)));
        // currently toSAMRecord does not have an option to keep a cigar with a combo of H and S operators.
        // so in case this is added.
        if (samRecord.getCigar().containsOperator(CigarOperator.H) &&
                samRecord.getCigar().containsOperator(CigarOperator.S)) {
            Assert.fail("test case with a combo of H and S operators, either you should silence this error or enhance toSAMRecord to handle such scenario");
        } else {
            Assert.assertEquals(backSamRecord.getCigar(), samRecord.getCigar());
            Assert.assertEquals(backSamRecord.getReadBases(), samRecord.getReadBases());
            Assert.assertEquals(backSamRecord.getBaseQualities(), samRecord.getBaseQualities());
        }
    }

    @DataProvider(name = "alignmentIntervalStrings")
    public Object[][] alignmentIntervalStrings() {
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[]{ "chr1", 10, Strand.NEGATIVE, "10M1I30M100H", 10, 3, 2 });
        result.add(new Object[]{ "chrX", 10_000_000, Strand.POSITIVE, "31H10S10M1I30M230N4M100H", 34, 31, 0 });
        result.add(new Object[]{ "chr20", 3456, Strand.POSITIVE, "31M", 3, 310, 5 });
        return result.toArray(new Object[result.size()][]);
    }

    @Test(dataProvider = "alignmentIntervalStrings", groups = "sv")
    public void testAlignmentIntervalStrings(final String contig, final int start, final Strand strand, final String cigarString, final int mq, final int nm, final int as) {
        final String fullStr = String.join(",", contig, "" + start, strand == Strand.NEGATIVE ? "-" : "+", cigarString, "" + mq, "" + nm, "" + as);
        final AlignmentInterval fullInterval = new AlignmentInterval(fullStr);
        Assert.assertEquals(fullInterval.referenceSpan.getContig(), contig);
        Assert.assertEquals(fullInterval.referenceSpan.getStart(), start);
        Assert.assertEquals(fullInterval.forwardStrand, strand == Strand.POSITIVE);
        Assert.assertEquals(fullInterval.cigarAlong5to3DirectionOfContig,
                fullInterval.forwardStrand ? TextCigarCodec.decode(cigarString) : CigarUtils.invertCigar(TextCigarCodec.decode(cigarString)));
        Assert.assertEquals(fullInterval.mapQual, mq);
        Assert.assertEquals(fullInterval.mismatches, nm);
        Assert.assertEquals(fullInterval.alnScore, as);

        final String basicStr = String.join(",", contig, "" + start, strand == Strand.NEGATIVE ? "-" : "+", cigarString, "" + mq);
        final AlignmentInterval basicInterval = new AlignmentInterval(basicStr);
        Assert.assertEquals(basicInterval.referenceSpan.getContig(), contig);
        Assert.assertEquals(basicInterval.referenceSpan.getStart(), start);
        Assert.assertEquals(basicInterval.forwardStrand, strand == Strand.POSITIVE);
        Assert.assertEquals(basicInterval.cigarAlong5to3DirectionOfContig,
                basicInterval.forwardStrand ? TextCigarCodec.decode(cigarString) : CigarUtils.invertCigar(TextCigarCodec.decode(cigarString)));
        Assert.assertEquals(basicInterval.mapQual, mq);
        Assert.assertEquals(basicInterval.mismatches, AlignmentInterval.NO_NM);
        Assert.assertEquals(basicInterval.alnScore, AlignmentInterval.NO_AS);

    }

    @Test(dataProvider = "randomValidCigars", groups = "sv")
    public void testSoftClip(final Cigar cigar) {
        final Cigar actual = AlignmentInterval.softOrHardReclip(cigar, CigarOperator.S);
        final Cigar expected = CigarUtils.combineAdjacentCigarElements(new Cigar(
                cigar.getCigarElements().stream()
                        .map(ce -> ce.getOperator().isClipping() ? new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP) : ce)
                        .collect(Collectors.toList())
        ));

        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "randomValidCigars", groups = "sv")
    public void testHardClip(final Cigar cigar) {
        final Cigar actual = AlignmentInterval.softOrHardReclip(cigar, CigarOperator.H);
        final Cigar expected = CigarUtils.combineAdjacentCigarElements(new Cigar(
                cigar.getCigarElements().stream()
                        .map(ce -> ce.getOperator().isClipping() ? new CigarElement(ce.getLength(), CigarOperator.HARD_CLIP) : ce)
                        .collect(Collectors.toList())
        ));
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name = "randomValidCigars")
    public static Object[][] randomValidCigars() {
        final List<Cigar> cigars = CigarTestUtils.randomValidCigars(new Random(13), 1000, 10, 100, new Cigar());
        return cigars.stream().map(x -> new Object[] { x }).toArray(Object[][]::new);
    }

    @Test(groups = "sv", expectedExceptions = IllegalArgumentException.class)
    public void testConstructFromUnmappedRead() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead artificialUnmappedRead = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{}, new byte[]{});
        new AlignmentInterval(artificialUnmappedRead);
    }

    @Test(groups = "sv", expectedExceptions = IllegalArgumentException.class)
    public void testConstructFromUnmappedSAM() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final SAMRecord unmappedSam = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{}, new byte[]{}).convertToSAMRecord(header);
        new AlignmentInterval(unmappedSam);
    }

    @DataProvider(name = "TestDataForReadIntervalAlignedToRefSpan")
    Object[][] testDataForReadIntervalAlignedToRefSpan() {
        final List<Object[]> data = new ArrayList<>(20);

        final Tuple2<Integer, Integer> noOverlap = new Tuple2<>(-1,-1);

        SimpleInterval alignmentRefSpan = new SimpleInterval("chr1",4938770,4939497);
        AlignmentInterval alignment = new AlignmentInterval(alignmentRefSpan, 1,728, TextCigarCodec.decode("728M61S"),true, 60, 0, 728, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleInterval otherRefSpan = new SimpleInterval("chr2", 4938770, 4939497);
        data.add(new Object[]{alignment, otherRefSpan, noOverlap});

        otherRefSpan = new SimpleInterval("chr1", 1000, 2000);
        data.add(new Object[]{alignment, otherRefSpan, noOverlap});

        Tuple2<Integer, Integer> readInterval = new Tuple2<>(1, 728);
        otherRefSpan = new SimpleInterval("chr1",4938770,4939497);
        data.add(new Object[]{alignment, otherRefSpan, readInterval});

        alignmentRefSpan = new SimpleInterval("chr1",9170350,9171390);
        alignment = new AlignmentInterval(alignmentRefSpan, 1,1041, TextCigarCodec.decode("1041M1298H"),false, 60, 4, 1021, ContigAlignmentsModifier.AlnModType.NONE);
        otherRefSpan = new SimpleInterval("chr1",934806,935261);
        data.add(new Object[]{alignment, otherRefSpan, noOverlap});



        alignmentRefSpan = new SimpleInterval("chr1",66659809,66660176);
        alignment = new AlignmentInterval(alignmentRefSpan, 1,354, TextCigarCodec.decode("124M10D106M3D16M2I75M3D31M241H"),true, 60, 35, 185, ContigAlignmentsModifier.AlnModType.NONE);
        otherRefSpan = new SimpleInterval("chr1",66659958,66660262);
        readInterval = new Tuple2<>(140,354);
        data.add(new Object[]{alignment, otherRefSpan, readInterval});

        alignmentRefSpan = new SimpleInterval("chr1",156328046,156328757);
        alignment = new AlignmentInterval(alignmentRefSpan, 1,712, TextCigarCodec.decode("712M444S"),false, 60, 2, 702, ContigAlignmentsModifier.AlnModType.NONE);
        otherRefSpan = new SimpleInterval("chr1",156327744,156328331);
        readInterval = new Tuple2<>(427,712);
        data.add(new Object[]{alignment, otherRefSpan, readInterval});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "TestDataForReadIntervalAlignedToRefSpan", groups = "sv")
    public void testComputeReadIntervalAlignedToRefSpan(final AlignmentInterval alignment, final SimpleInterval otherRefSpan,
                                                        final Tuple2<Integer, Integer> expectedReadInterval) {
        Assert.assertEquals(alignment.readIntervalAlignedToRefSpan(otherRefSpan), expectedReadInterval);
    }

    @DataProvider(name = "testDataForContainsGapOfEqualOrLargerSize")
    private Object[][] testDataForContainsGapOfEqualOrLargerSize(){
        final List<Object[]> data = new ArrayList<>(20);

        final AlignmentInterval one = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00002\t0\tchr1\t180969\t60\t400S237M58D320M4D76M\t*\t0\t0\tCCCTGCAATGTCCCTAGCTGCCAGCAGGCGGCGTGCCACCACTATACAGTAAGCAAGAGGGCCCTGCAGTGCCCCGGCGCCAGCAGGGGGCGCTGGCCACCACTCTAAGCAAGAGAGCCCTGCAGTTGCCCTAGTCGCCAGCAGGGGGCGCCCTGGCACAGCACCGTGAGCAAGCGGGTCCTGTAGTGCCCGGCTGCAAGCAAGGGGCTGTCGATCCCGGCGTTTCGGATTACTGAGGTTCCACCCGTCTCTGCGCCGCGCCGCCGTGACGTGAGTTTCTGCGCGTGCACGGCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGCGCCCTCCCCTCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCAACCCTAACCCCAACCCTAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGAGGAGGCGTGGCACAGGCGCAGAGACACATGCTAGCGCGCCCAGGGGAGGAGGCGTGGCGCAGGCGCAGAGAGGCGCGCCGTGCTGCCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGGTGGAGGCGTGGCGCAGGCGCAGAGACGCACGCCTACGGGCGGGGTTGGGGGGGGCGTGTGTTACAGGAGCAAAGTCGCACGGCGCCGGGCTGGGGGCGGGGGCGGGGGGGCGCCGTGCACGCGCAGAAACTCACGTCACGGCGGCGCGGCGCAGAGACGGGTGGAACCTCAGT\t*\tSA:Z:chr1,181578,-,716S317M,6,1;chr2,32916352,-,640S99M294S,0,5;\tMD:Z:11T132A92^GCCGGCGCAGGCGCAGAGAGGCGCGCCTCGCCGGCGCAGGCGCAGAGAGGCGCGCCGG320^GGGG5G70\tRG:Z:GATKSVContigAlignments\tNM:i:65\tAS:i:524\tXS:i:426",
                true);
        data.add(new Object[]{one, 50, true});

        final AlignmentInterval two = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000001:tig00005\t16\tchr2\t113602055\t60\t129M29D136M\t*\t0\t0\tCAAGAGGACCCTGCAATGTCCCTAGCTGCCAGCAGGCGGCGTGCCACCACTATACAGTAAGCAAGAGGGCCCTGCAGTGCCCCGGCGCCAGCAGGGGGCGCTGGCGACCACTGTAAGCAAGAGAGCCCTGCGCCTCTCTGCGCCGGCGCGGCGCGGGGTGCCTTTGCGACGGCGGAGTTGCGTTCTCCTCAGCACAGACCCGGAGAGCACCGCGAGGGCGGACCTGCGTTGTCCTCTGCACAGATTTCAGTGGTACTGCGAAGGC\t*\tMD:Z:129^GCGCCTCTCTGCGCCGGCGCCGGCGCGGC136\tRG:Z:GATKSVContigAlignments\tNM:i:29\tAS:i:2",
                false);
        data.add(new Object[]{two, 50, false});

        final AlignmentInterval three = TestUtilsForAssemblyBasedSVDiscovery.fromSAMRecordString("asm000004:tig00026\t0\tchr1\t224015348\t40\t104S109M20I190M10D79M15D69M\t*\t0\t0\tAGGAGGGAATGTATTCGAGTGGAATGGAAAGGAATGGAATCAACCATAGTGGAATGGAAGGGAATGGAATGGAATGGTAACTAATAGAATGGAATCAACCCGAGTGGAATGGAATGGAAAGGACTGGAATGGAATGGAATGGAATGGAATGGACTCAGATGGAATGGAATGGAATGGACTCGAAAGGATTGGGATGGAATACAATGGAATGGTCTCGAATGGAATGGAATGGACTCGAATGGAATGGAATGCAATGGAATGGACTCAAATGGAATGGAATGGAACTGACTCGAATGGAATTGAATGGAATGGACCCGAATGGAATGGAATGGAATGGACTGGGCTCAAATTGAATGGAATGGAAAAGAATGGAATGGAATAGAATGGACTGGAATGTAATGAGTTTGGAATGGACTTGAATGCAATGGAATGGAATGGACTCAAATGGAATAGCATGGAATGGAATGGACTCAAATGCATTGGAATGGAATGGACTTGAATGGAATGGAATGGAGTTGAATGGACTCATATGGAATGGAATGGCATTGAATGGACTCGAATGGAATAGAAT\t*\tSA:Z:chr5,49661411,-,418S150M3S,4,9;chrUn_KN707896v1_decoy,5378,-,418S153M,0,3;\tMD:Z:7C7T5C277^AATGGAATGG79^GAATGGAATGGACTC69\tRG:Z:GATKSVContigAlignments\tNM:i:48\tAS:i:339\tXS:i:308",
                true);
        data.add(new Object[]{three, 10, true});
        data.add(new Object[]{three, 20, true});
        data.add(new Object[]{three, 21, false});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "testDataForContainsGapOfEqualOrLargerSize")
    public void testContainsGapOfEqualOrLargerSize(final AlignmentInterval alignment, final int gapSize,
                                                   final boolean expectedResult) {
        Assert.assertEquals(alignment.containsGapOfEqualOrLargerSize(gapSize), expectedResult);
    }

    @DataProvider(name = "forTestCtorArgChecking")
    private Object[][] forTestCtorArgChecking() {
        final List<Object[]> data = new ArrayList<>(20);

        data.add(new Object[]{TextCigarCodec.decode("1155M1154S"), new SimpleInterval("chr22", 47043976, 47045130), 1, 1155, null});
        data.add(new Object[]{TextCigarCodec.decode("1424M1424S"), new SimpleInterval("chr15", 80355809, 80357232), 1, 1424, null});

        data.add(new Object[]{TextCigarCodec.decode("1155M1154S"), new SimpleInterval("chr22", 47043976, 47045131), 1, 1155, IllegalArgumentException.class});
        data.add(new Object[]{TextCigarCodec.decode("1424M1424S"), new SimpleInterval("chr15", 80355809, 80357232), 1, 1429, IllegalArgumentException.class});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forTestCtorArgChecking")
    @SuppressWarnings("rawtypes")
    public void testCtorArgChecking(final Cigar cigar, final SimpleInterval referenceSpan, final int readStart, final int readEnd,
                                    final Class expectedExceptionClass) {
        try {
            AlignmentInterval.checkValidArgument(cigar, referenceSpan, readStart, readEnd);
        } catch (final Exception e) {
            Assert.assertEquals(e.getClass(), expectedExceptionClass);
        }
    }
}
