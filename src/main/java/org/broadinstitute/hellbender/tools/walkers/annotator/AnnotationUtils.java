package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

public final class AnnotationUtils {
    private AnnotationUtils(){}

    /**
     * Helper function to parse the list into the annotation string
     * @param valueList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    public static String encodeValueList(final List<Double> valueList, final String precisionFormat ) {
        List<String> outputList = new ArrayList<>();
        for (Double d : valueList) {
            outputList.add(String.format(precisionFormat, d));
        }
        return StringUtils.join(outputList, ",");
    }

    /**
     * Helper function to convert a List of Strings to a comma-separated String
     * @param stringList the ArrayList with String data
     * @return a comma-separated String
     */
    public static String encodeStringList( final List<String> stringList) {
        return StringUtils.join(stringList, ",");
    }

    /**
     * Get the position of a variant within a read with respect to the closer end, accounting for hard clipped bases and low quality ends
     * Used by ReadPosRankSum annotations
     *
     * @param read  a read containing the variant
     * @param initialReadPosition   the position based on the modified, post-hard-clipped CIGAR
     * @return read position
     */
    public static int getFinalVariantReadPosition(final GATKRead read, final int initialReadPosition) {
        final int numAlignedBases = getNumAlignedBases(read);

        int readPos = initialReadPosition;
        //TODO: this doesn't work for the middle-right position if we index from zero
        if (initialReadPosition > numAlignedBases / 2) {
            readPos = numAlignedBases - (initialReadPosition + 1);
        }
        return readPos;

    }

    /**
     *
     * @param read  a read containing the variant
     * @return  the number of hard clipped and low qual bases at the read start (where start is the leftmost end w.r.t. the reference)
     */
    public static int getNumClippedBasesAtStart(final GATKRead read) {
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        final CigarElement first = c.getCigarElement(0);

        int numStartClippedBases = 0;
        if (first.getOperator() == CigarOperator.H) {
            numStartClippedBases = first.getLength();
        }
        final byte[] unclippedReadBases = read.getBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        //TODO: this code may not even get used because HaplotypeCaller already hard clips low quality tails
        for (int i = numStartClippedBases; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < 20) { //TODO the 20 hard value here is in lieu of PairHMM.BASE_QUALITY_SCORE_THRESHOLD in order to directly match GATK3 output
                numStartClippedBases++;
            } else {
                break;
            }
        }

        return numStartClippedBases;
    }


    /**
     *
     * @param read  a read containing the variant
     * @return  number of non-hard clipped, aligned bases (excluding low quality bases at either end)
     */
    //TODO: this is bizarre -- this code counts hard clips, but then subtracts them from the read length, which already doesn't count hard clips
    public static int getNumAlignedBases(final GATKRead read) {
        return read.getLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }

    /**
     *
     * @param read  a read containing the variant
     * @return  number of hard clipped and low qual bases at the read end (where end is right end w.r.t. the reference)
     */
    public static int getNumClippedBasesAtEnd(final GATKRead read) {
        // check for hard clips (never consider these bases):
        final Cigar c = read.getCigar();
        CigarElement last = c.getCigarElement(c.numCigarElements() - 1);

        int numEndClippedBases = 0;
        if (last.getOperator() == CigarOperator.H) {
            numEndClippedBases = last.getLength();
        }
        final byte[] unclippedReadBases = read.getBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        //TODO: this code may not even get used because HaplotypeCaller already hard clips low quality tails
        for (int i = unclippedReadBases.length - numEndClippedBases - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < 20) { //TODO the 20 hard value here is in lieu of PairHMM.BASE_QUALITY_SCORE_THRESHOLD in order to directly match GATK3 output

                numEndClippedBases++;
            } else {
                break;
            }
        }

        return numEndClippedBases;
    }
}
