package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link SingleSequenceReferenceAligner}
 */
public class SingleContigReferenceAlignerUnitTest extends BaseTest {

    private static final int READ_LENGTH = 250;
    private static final int NUM_ALIGNS = 1000;
    private static final String REF_NAME = "ref00";

    @Test(expectedExceptions = IllegalStateException.class, expectedExceptionsMessageRegExp = "^.*closed.*$")
    public void testClosedException() {
        final RandomDNA rdn = new RandomDNA(131);
        final byte[] refBases = rdn.nextBases(1000);
        final SingleSequenceReferenceAligner<byte[], AlignedContig> aligner = SingleSequenceReferenceAligner.contigsAligner("test", refBases,
                a -> "contig", a -> a);
        try {
            aligner.close();
        } catch (final IOException ex) {
            Assert.fail("unexpected exception when closing");
        }
        aligner.align(Collections.singletonList(Arrays.copyOfRange(refBases, 100, 200)));
    }


    @Test(dataProvider = "testAlignmentData")
    public void testAlignment(final boolean pairAlignment, final byte[] reference, final String referenceName)
        throws IOException
    {
        try (final SingleSequenceReferenceAligner<byte[], AlignedContig> aligner =  SingleSequenceReferenceAligner.contigsAligner(referenceName, reference,
                a -> "ctg", a -> a);) {
            Assert.assertNotNull(aligner.getAligner());
            if (pairAlignment) {
                aligner.getAligner().alignPairs();
            }
            final Random rdn = new Random(13111);
            final RandomDNA rdnDNA = new RandomDNA(rdn);
            final List<List<AlignmentInterval>> expected = new ArrayList<>(NUM_ALIGNS << 1);
            final List<byte[]> seqs = new ArrayList<>(NUM_ALIGNS << 1);
            for (int i = 0; i < NUM_ALIGNS; i++) {
                final int start = rdn.nextInt(reference.length - READ_LENGTH) + 1;
                final boolean forward = rdn.nextBoolean();
                final boolean insert = rdn.nextDouble() < 0.1;
                final boolean deletion = !insert && rdn.nextDouble() < 0.1;
                final int indelLength = insert || deletion ? rdn.nextInt(10) + 10 : 0;
                final int indelStart = insert || deletion ? rdn.nextInt((int) (READ_LENGTH  * .50)) + 25 : -1;
                final int end = start + READ_LENGTH - 1;
                final byte[] templateSeq = Arrays.copyOfRange(reference, start - 1, end);
                final byte[] actualSeq;
                if (insert) {
                    actualSeq = Arrays.copyOf(templateSeq, templateSeq.length + indelLength);
                    System.arraycopy(actualSeq, indelStart - 1, actualSeq, indelStart - 1 + indelLength, templateSeq.length - indelStart + 1);
                    rdnDNA.nextBases(actualSeq, indelStart - 1, indelLength);
                } else if (deletion) {
                    actualSeq = Arrays.copyOf(templateSeq, templateSeq.length - indelLength);
                    System.arraycopy(templateSeq, indelStart - 1 + indelLength, actualSeq, indelStart - 1, templateSeq.length - indelStart + 1 - indelLength);
                } else {
                    actualSeq = templateSeq.clone();
                }

                while (insert && actualSeq[indelStart - 1] == actualSeq[indelStart + indelLength - 1]) {
                    actualSeq[indelStart + indelLength - 1] = rdnDNA.nextBase();
                }
                if (!forward) {
                    SequenceUtil.reverseComplement(actualSeq);
                }
                seqs.add(actualSeq);
                final Cigar cigar = (!insert && !deletion) ? TextCigarCodec.decode(READ_LENGTH + "M"):
                        (insert ? TextCigarCodec.decode( "" + (indelStart - 1) + "M" + indelLength + "I" + (READ_LENGTH - indelStart + 1) + "M")
                                : TextCigarCodec.decode( "" + (indelStart - 1) + "M" + indelLength + "D" + (READ_LENGTH - indelStart + 1 - indelLength) + "M"));

                expected.add(Collections.singletonList(new AlignmentInterval(new SimpleInterval(referenceName, start, end), 1, actualSeq.length, !forward ? CigarUtils.invertCigar(cigar) : cigar
                        , forward, 0, 0, 0, null)));
            }
            final List<AlignedContig> results = aligner.align(seqs);
            final Map<byte[], AlignedContig> mapResult = aligner.align(seqs, (b, a) -> b);
            Assert.assertEquals(results, new ArrayList<>(mapResult.values()));
            Assert.assertEquals(new ArrayList<>(mapResult.keySet()), mapResult.values().stream().map(AlignedContig::getContigSequence).collect(Collectors.toList()));
            for (int i = 0; i < NUM_ALIGNS; i++) {
                final List<AlignmentInterval> actualValue = results.get(i).getAlignments();
                final List<AlignmentInterval> expectedValue = expected.get(i);
                Assert.assertEquals(actualValue.size(), 1);
                Assert.assertEquals(actualValue.get(0).forwardStrand, expectedValue.get(0).forwardStrand);
                Assert.assertEquals(actualValue.get(0).referenceSpan, expectedValue.get(0).referenceSpan, expectedValue.get(0).cigarAlong5to3DirectionOfContig.toString());
                Assert.assertEquals(actualValue.get(0).startInAssembledContig, expectedValue.get(0).startInAssembledContig);
                Assert.assertEquals(actualValue.get(0).endInAssembledContig, expectedValue.get(0).endInAssembledContig);
                final Cigar expectedCigar = expectedValue.get(0).cigarAlong5to3DirectionOfContig;
                final Cigar actualCigar = actualValue.get(0).cigarAlong5to3DirectionOfContig;
                if (!expectedCigar.equals(actualCigar)) { // small differences may occur due to ambiguous indel location. So we check that they are small differences indeed:
                    Assert.assertEquals(expectedCigar.numCigarElements(), actualCigar.numCigarElements()); // same number of elements
                    Assert.assertEquals(expectedCigar.getCigarElements().stream().map(CigarElement::getOperator).collect(Collectors.toList()),
                            actualCigar.getCigarElements().stream().map(CigarElement::getOperator).collect(Collectors.toList())); // same operators sequence.
                    // then we check the total lengths per operator (must be the same):
                    final Map<CigarOperator, Integer> expectedLengthByOperator = expectedCigar.getCigarElements().stream()
                            .collect(Collectors.groupingBy(CigarElement::getOperator,
                                    Collectors.reducing(0, CigarElement::getLength, (a, b) -> a + b)));
                    final Map<CigarOperator, Integer> actualLengthByOperator = actualCigar.getCigarElements().stream()
                            .collect(Collectors.groupingBy(CigarElement::getOperator,
                                    Collectors.reducing(0, CigarElement::getLength, (a, b) -> a + b)));
                    Assert.assertEquals(actualLengthByOperator, expectedLengthByOperator);
                    // finally we don't allow more than 5 bases length difference for any given element.
                    for (int j = 0; j < expectedCigar.numCigarElements(); j++) {
                        Assert.assertTrue(Math.abs(expectedCigar.getCigarElement(j).getLength() - actualCigar.getCigarElement(j).getLength()) < 10, "actual: " + actualCigar + " != expected: " + expectedCigar);
                    }
                }
            }
        }
    }

    @DataProvider(name="testAlignmentData")
    public Object[][] testAlignmentData() {
        final List<Object[]> result = new ArrayList<>();
        final RandomDNA randomDNA = new RandomDNA(1301);
        result.add(new Object[]{ false, randomDNA.nextBases(1000), REF_NAME});
        result.add(new Object[]{ true, randomDNA.nextBases(10000), REF_NAME});
        return result.toArray(new Object[result.size()][]);
    }

}
