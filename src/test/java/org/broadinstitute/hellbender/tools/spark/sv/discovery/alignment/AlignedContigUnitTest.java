package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by valentin on 4/23/18.
 */
public class AlignedContigUnitTest extends GATKBaseTest {

    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    private static final ReferenceSparkSource REFERENCE = new ReferenceMultiSparkSource(twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);

    public static void toFasta(final File file, final ReferenceSparkSource reference, final AlignedContig contig, final int padding, final boolean matchedBasedAsDots) throws IOException {
        try (final Writer writer = new FileWriter(file)){
            toFasta(writer, reference, contig, padding, matchedBasedAsDots);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, ex.getMessage(), ex);
        }
    }

    private static void toFasta(final Writer writer, final ReferenceSparkSource reference, final AlignedContig contig, final int padding, final boolean matchedBasesAsDots) throws IOException {
        writer.append(toFastaString(reference, contig, padding, matchedBasesAsDots));
    }

    private static String toFastaString(final ReferenceSparkSource reference, final AlignedContig contig, final int padding, final boolean matchedBasesAsDots)
            throws IOException
    {
        Utils.nonNull(reference);
        Utils.nonNull(contig);
        final SAMSequenceDictionary dictionary = Utils.nonNull(reference.getReferenceSequenceDictionary());
        final List<SimpleInterval> refereneSpans = contig.getAlignments().stream()
                .map(ai -> {
                    final SimpleInterval referenceSpan = ai.referenceSpan;
                    final int leftSoftClip = ai.cigarAlongReference().getCigarElement(0).getOperator() == CigarOperator.S ?
                            ai.cigarAlongReference().getCigarElement(0).getLength() : 0;
                    final int rightSoftClip = ai.cigarAlongReference().getLastCigarElement().getOperator() == CigarOperator.S ?
                            ai.cigarAlongReference().getLastCigarElement().getLength() : 0;
                    if (leftSoftClip == 0 && rightSoftClip == 0 && padding == 0) {
                        return referenceSpan;
                    } else {
                        return new SimpleInterval(referenceSpan.getContig(), referenceSpan.getStart() - leftSoftClip - padding, referenceSpan.getEnd() + rightSoftClip + padding);
                    }
                }).collect(Collectors.toList());
        Collections.sort(refereneSpans, IntervalUtils.getDictionaryOrderComparator(reference.getReferenceSequenceDictionary()));
        final Deque<SimpleInterval> mergedSpans = new ArrayDeque<>(refereneSpans.size());
        SimpleInterval last;
        mergedSpans.add(last = refereneSpans.get(0));
        for (int i = 1; i < refereneSpans.size(); i++) {
            final SimpleInterval next = refereneSpans.get(i);
            if (last.overlapsWithMargin(next, 1)) {
                mergedSpans.removeLast();
                mergedSpans.addLast(last = last.mergeWithContiguous(next));
            } else {
                mergedSpans.addLast(last = next);
            }
        }
        final StringBuilder output = new StringBuilder(10000);
        for (final SimpleInterval referenceSpan : mergedSpans) {
            final Stream<ImmutablePair<Integer, Integer>> insertions = contig.getAlignments().stream()
                    .filter(ai -> referenceSpan.contains(ai.referenceSpan))
                    .filter(ai -> ai.cigarAlong5to3DirectionOfContig.containsOperator(CigarOperator.I))
                    .flatMap(ai -> {
                        final List<ImmutablePair<Integer, Integer>> insertionAndLengths = new ArrayList<>(ai.cigarAlong5to3DirectionOfContig.numCigarElements());
                        int referencePostion = ai.referenceSpan.getStart();
                        for (final CigarElement e : ai.cigarAlongReference()) {
                            if (e.getOperator().consumesReferenceBases()) {
                                referencePostion += e.getLength();
                            } else if (e.getOperator() == CigarOperator.I) {
                                insertionAndLengths.add(new ImmutablePair<>(referencePostion, e.getLength()));
                            }
                        }
                        return insertionAndLengths.stream();
                    });
            final Map<Integer, Integer> insertionPositionAndMaxLength = insertions
                    //           .collect(Collectors.groupingBy(p -> p.getLeft()));
                    .collect(Collectors.groupingBy(ImmutablePair::getLeft, Collectors.mapping(ImmutablePair::getRight, Collectors.reducing(0, Math::max))));
            final int width = referenceSpan.size() + insertionPositionAndMaxLength.values().stream().mapToInt(i -> i).sum();
            Utils.nonNull(referenceSpan.getContig());
            final int beyondContigEndPadding = referenceSpan.getEnd() - dictionary.getSequence(referenceSpan.getContig()).getSequenceLength();
            final SimpleInterval actualReferenceSpan;
            output.append('>');
            if (beyondContigEndPadding <= 0) {
                output.append(referenceSpan);
                actualReferenceSpan = referenceSpan;
            } else {
                actualReferenceSpan = new SimpleInterval(referenceSpan.getContig(), referenceSpan.getStart(), dictionary.getSequence(referenceSpan.getContig()).getSequenceLength());
                output.append(actualReferenceSpan);
                output.append("\t(").append(beyondContigEndPadding).append(" gap padded beyond end");
            }
            output.append('\n');
            final int refSeqOffset = output.length();
            final byte[] refBases = reference.getReferenceBases(actualReferenceSpan).getBases();
            for (int i = 0, pos = referenceSpan.getStart(); pos <= actualReferenceSpan.getEnd(); ++i, ++pos) {
                final int insertLength = insertionPositionAndMaxLength.getOrDefault(pos, 0);
                if (insertLength > 0) {
                    for (int j = 0; j < insertLength; ++j) {
                        output.append('-');
                    }
                }
                output.append(Character.toUpperCase((char) refBases[i]));
            }
            for (int i = 0; i < beyondContigEndPadding; i++) {
                output.append('-');
            }
            final byte[] referenceBasesWithGaps = output.subSequence(refSeqOffset, output.length()).toString().getBytes();
            output.append('\n');

            contig.getAlignments().stream()
                    .filter(ai -> actualReferenceSpan.contains(ai.referenceSpan.expandWithinContig(padding, dictionary)))
                    .forEach(ai -> {
                        output.append('>').append(contig.getContigName()).append(':').append(ai.startInAssembledContig).append('-').append(ai.endInAssembledContig);
                        output.append('\t');
                        ai.appendSATagString(output);
                        output.append('\n');
                        final List<CigarElement> cigar = new ArrayList<>(ai.cigarAlong5to3DirectionOfContig.numCigarElements() + 2);

                        // Since we need to skip to the first aligned reference position we artificailly add a deletion operation
                        // that would consume those references bases before the first aligned based to the contig:
                        if (actualReferenceSpan.getStart() < ai.referenceSpan.getStart()) {
                            cigar.add(new CigarElement(ai.referenceSpan.getStart() - actualReferenceSpan.getStart(), CigarOperator.D));
                        }
                        cigar.addAll(ai.cigarAlongReference().getCigarElements());
                        int linePos = 0;
                        final int contigBaseIncrease = ai.forwardStrand ? 1 : -1; // we move forward or backward in the contig sequence?
                        int nextContigBase = (ai.forwardStrand ? ai.startInAssembledContig - CigarUtils.countLeftClippedBases(ai.cigarAlong5to3DirectionOfContig)
                                : ai.endInAssembledContig + CigarUtils.countRightSoftClippedBases(ai.cigarAlong5to3DirectionOfContig)) - 1;
                        for (final CigarElement e : cigar) {
                            final int length = e.getLength();
                            final CigarOperator o = e.getOperator();
                            if (o.consumesReferenceBases() && !o.consumesReadBases()) {
                                for (int remaining = length; remaining > 0; ++linePos) {
                                    if (referenceBasesWithGaps[linePos] != '-') {
                                        remaining--;
                                    } else {
                                        System.err.append('.');
                                    }
                                    output.append('-');
                                }
                            } else if (o == CigarOperator.S) {
                                if (output.length() == 0 || output.charAt(output.length() - 1) == '-') {
                                    for (; referenceBasesWithGaps[linePos] == '-'; ++linePos) {
                                        output.append('-');
                                    } // skip to the next ref base.
                                    output.setLength(output.length() - length);
                                    linePos -= length;
                                }
                                for (int i = 0; i < length; i++, ++linePos) {
                                    final byte contigBase = ai.forwardStrand
                                            ? contig.getContigSequence()[nextContigBase]
                                            : Nucleotide.complement(contig.getContigSequence()[nextContigBase]);
                                    output.append(Character.toLowerCase((char) contigBase));
                                    nextContigBase += contigBaseIncrease;
                                }
                            } else if (o.consumesReadBases() && o.consumesReferenceBases()) {
                                for (int remaining = length; remaining > 0; ++linePos) {
                                    if (referenceBasesWithGaps[linePos] != '-') {
                                        remaining--;
                                        final byte contigBase = ai.forwardStrand
                                                ? contig.getContigSequence()[nextContigBase]
                                                : Nucleotide.complement(contig.getContigSequence()[nextContigBase]);
                                        if (matchedBasesAsDots && Nucleotide.same(referenceBasesWithGaps[linePos], contigBase)) {
                                            output.append('.');
                                        } else {
                                            output.append(Character.toUpperCase((char) contigBase));
                                        }
                                        nextContigBase += contigBaseIncrease;
                                    } else {
                                        output.append('-');
                                    }
                                }
                            } else if (o.consumesReadBases()) {// && !o.consumesReferenceBases())
                                for (int i = 0; i < length; i++, ++linePos) {
                                    final byte contigBase = ai.forwardStrand
                                            ? contig.getContigSequence()[nextContigBase]
                                            : Nucleotide.complement(contig.getContigSequence()[nextContigBase]);
                                    output.append(Character.toUpperCase((char) contigBase));
                                    nextContigBase += contigBaseIncrease;
                                }
                                while (referenceBasesWithGaps[linePos] == '-') {
                                    output.append('-');
                                    linePos++;
                                }
                            }
                        }
                        for (; linePos < width; ++linePos) {
                            output.append('-');
                        }
                        output.append('\n');
                    });
            while (Character.isWhitespace(output.charAt(output.length() - 1))) {
                output.setLength(output.length() - 1);
            }
        }
        return output.toString();
    }


    @Test(dataProvider = "toFastaData")
    public void testToFasta(final AlignedContig contig, final int padding, final boolean matchesAsDot, final String expected)
        throws IOException {
        final String actual = toFastaString(REFERENCE, contig, padding, matchesAsDot);
//        if (!actual.equals(expected)) {
//            Assert.assertEquals(actual.length(), expected.length());
//            for (int i = 0; i < actual.length(); i++) {
//                Assert.assertEquals(actual.charAt(i), expected.charAt(i), "Character at " + (i+1) + "/" + (actual.length()) +  " is different " + actual.charAt(i) + " != " + expected.charAt(i));
 //           }
 //       }
        Assert.assertEquals(actual, expected);
    }

    @DataProvider(name="toFastaData")
    public Object[][] toFastaData() throws IOException {
        int i = 0;
        int j = 0;
        try {
            final Random rdn = new Random(13);
            final RandomDNA randomDNA = new RandomDNA(rdn);

            final List<Object[]> result = new ArrayList<>();
            final SAMSequenceRecord chr1 = REFERENCE.getReferenceSequenceDictionary().getSequence(0);
            final SAMSequenceRecord chr2 = REFERENCE.getReferenceSequenceDictionary().getSequence(1);

            final ReferenceBases chr1_101000_102000 = REFERENCE.getReferenceBases(new SimpleInterval(chr1.getSequenceName(), 101000, 102000));
            final ReferenceBases chr1_101000_102000_10padding = REFERENCE.getReferenceBases(new SimpleInterval(chr1.getSequenceName(), 101000 - 10, 102000 + 10));
            final AlignedContig chr1_101000_102000_pertch_match = new AlignedContig("ctg001", chr1_101000_102000.getBases(), Collections.singletonList(new AlignmentInterval(chr1.getSequenceName() + ",101000,+,1001M,40,0,40")));
            result.add(new Object[]{chr1_101000_102000_pertch_match, 0, false,
                    chr1_101000_102000.appendBasesTo(new StringBuilder(10000).append('>').append(chr1.getSequenceName()).append(":101000-102000\n"), chr1_101000_102000.getInterval().getStart(), chr1_101000_102000.getInterval().getEnd())
                            .append("\n>ctg001:1-1001\t").append(chr1.getSequenceName()).append(",101000,+,1001M,40,0,40\n").append(chr1_101000_102000.toBaseString()).toString()});

            result.add(new Object[]{chr1_101000_102000_pertch_match, 10, true,
                    '>' + chr1.getSequenceName() + ":100990-102010\n" + chr1_101000_102000_10padding.toBaseString() + '\n' +
                            ">ctg001:1-1001\t" + chr1.getSequenceName() + ",101000,+,1001M,40,0,40\n" + StringUtils.repeat("-", 10) + StringUtils.repeat(".", 1001) + StringUtils.repeat("-", 10)});

            final byte[] originalBases = chr1_101000_102000.getBases().clone();
            byte[] referenceBases = chr1_101000_102000_10padding.getBases().clone();
            byte[] mutatedBases = originalBases.clone();
            final List<CigarElement> cigar = new ArrayList<>();
            double mutationRate = 0.01;
            char[] dotsMutatedBases = new char[mutatedBases.length];
            double indelRate = 0.005;
            double extensionRate = 0.75;
            int matchStart = -1;
            int r;
            for (i = 0, j = 0, r = 10; i < originalBases.length; i++, r++) {
                if (rdn.nextDouble() < indelRate) {
                    if (matchStart >= 0) cigar.add(new CigarElement(j - matchStart, CigarOperator.M));
                    matchStart = -1;
                    final boolean deletion = rdn.nextBoolean();
                    int length = 1;
                    while (j + length < mutatedBases.length && rdn.nextDouble() < extensionRate) {
                        length++;
                    }
                    if (deletion) {
                        Arrays.fill(mutatedBases, j, j + length, (byte) '-');
                        Arrays.fill(dotsMutatedBases, j, j + length, '-');
                        cigar.add(new CigarElement(length, CigarOperator.D));
                        i += length - 1;
                        r += length - 1;
                    } else {
                        final byte[] insert = randomDNA.nextBases(length);
                        mutatedBases = Arrays.copyOf(mutatedBases, mutatedBases.length + length);
                        System.arraycopy(mutatedBases, j, mutatedBases, j + length, mutatedBases.length - length - j);
                        dotsMutatedBases = Arrays.copyOf(dotsMutatedBases, dotsMutatedBases.length + length);
                        System.arraycopy(dotsMutatedBases, j, dotsMutatedBases, j + length, dotsMutatedBases.length - length - j);
                        System.arraycopy(insert, 0, mutatedBases, j, length);
                        for (int k = 0; k < length; k++) {
                            dotsMutatedBases[j + k] = (char) insert[k];
                        }
                        referenceBases = Arrays.copyOf(referenceBases, referenceBases.length + length);
                        System.arraycopy(referenceBases, r, referenceBases, r + length, referenceBases.length - r - length);
                        Arrays.fill(referenceBases, r, r + length, (byte) '-');
                        cigar.add(new CigarElement(length, CigarOperator.I));
                        r += length - 1;
                    }
                    j += length;
                    indelRate = 0; // make sure we don't have consecutive indels.
                } else if (rdn.nextDouble() < mutationRate) {
                    indelRate = 0.01; // reset indel-rate.
                    if (matchStart < 0) matchStart = j;
                    mutatedBases[j] = randomDNA.mutate(originalBases[i]);
                    if (Nucleotide.same(originalBases[i],mutatedBases[j])) {
                        throw new IllegalArgumentException("");
                    }
                    dotsMutatedBases[j] = (char) mutatedBases[j];
                    mutationRate = 0.2; // make it more likely to have a follow up mutation.
                    j++;
                } else {
                    indelRate = 0.01; // reset indel-rate.
                    if (matchStart < 0) matchStart = j;
                    dotsMutatedBases[j] = '.';
                    mutationRate = Math.max(0.01, mutationRate * .5);
                    j++;
                }
            }
            if (matchStart >= 0) cigar.add(new CigarElement(mutatedBases.length - matchStart, CigarOperator.M));
            final byte[] mutatedBasesFinal = mutatedBases;
            final byte[] contigSequence = ArrayUtils.removeAll(mutatedBases, IntStream.range(0, mutatedBases.length)
                     .filter(idx -> mutatedBasesFinal[idx] == '-')
                     .toArray());
            final AlignedContig chr1_101000_102000_with_mutations = new AlignedContig("ctg001", contigSequence, Collections.singletonList(new AlignmentInterval(chr1.getSequenceName() + ",101000,+," + new Cigar(cigar) + ",40,0,40")));
            result.add(new Object[]{chr1_101000_102000_with_mutations, 10, true,
                    '>' + chr1.getSequenceName() + ":100990-102010\n" + new String(referenceBases) + '\n' +
                            ">ctg001:1-1001\t" + chr1.getSequenceName() + ",101000,+," + new Cigar(cigar) + ",40,0,40\n" + StringUtils.repeat("-", 10) + new String(dotsMutatedBases) + StringUtils.repeat("-", 10)});

            return result.stream().toArray(Object[][]::new);
        } catch (final RuntimeException ex) {
            throw ex;
        }
    }
}
