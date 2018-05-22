package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;


import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class RealignmentEngine {
    private final BwaMemAligner aligner;
    private final int maxReasonableFragmentLength;
    private final int minAlignerScoreDifference;
    private final int numberOfRegularContigs;

    public RealignmentEngine(final RealignmentArgumentCollection rfac) {
        final BwaMemIndex index = new BwaMemIndex(rfac.bwaMemIndexImage);
        maxReasonableFragmentLength = rfac.maxReasonableFragmentLength;
        minAlignerScoreDifference = rfac.minAlignerScoreDifference;
        numberOfRegularContigs = rfac.numRegularContigs;
        aligner = new BwaMemAligner(index);
        aligner.setMinSeedLengthOption(rfac.minSeedLength);
        aligner.setDropRatioOption((float) rfac.dropRatio);
        aligner.setSplitFactorOption((float) rfac.splitFactor);
        aligner.setFlagOption(BwaMemAligner.MEM_F_ALL);
    }

    public static boolean supportsVariant(final GATKRead read, final VariantContext vc, int indelStartTolerance) {
        final byte[] readBases = read.getBasesNoCopy();

        final int variantPositionInRead = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(),
                vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( variantPositionInRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), variantPositionInRead)) {
            return false;
        }

        for (final Allele allele : vc.getAlternateAlleles()) {
            final int referenceLength = vc.getReference().length();
            if (allele.length() == referenceLength) {   // SNP or MNP -- check whether read bases match variant allele bases
                if (allele.basesMatch(ArrayUtils.subarray(readBases, variantPositionInRead, Math.min(variantPositionInRead + allele.length(), readBases.length)))) {
                    return true;
                }
            } else {    // indel -- check if the read has the right CIGAR operator near position -- we don't require exact an exact match because indel representation is non-unique
                final boolean isDeletion = allele.length() < referenceLength;
                int readPosition = 0;    // offset by 1 because the variant position is one base before the first indel base
                for (final CigarElement cigarElement : read.getCigarElements()) {
                    if (Math.abs(readPosition - variantPositionInRead) <= indelStartTolerance) {
                        if (isDeletion && mightSupportDeletion(cigarElement) || !isDeletion && mightSupportInsertion(cigarElement)) {
                            return true;
                        }
                    } else {
                        readPosition += cigarElement.getLength();
                    }
                }
            }
        }
        return false;
    }

    private static boolean mightSupportDeletion(final CigarElement cigarElement) {
        final CigarOperator cigarOperator = cigarElement.getOperator();
        return cigarOperator == CigarOperator.D || cigarOperator == CigarOperator.S;
    }

    private final static boolean mightSupportInsertion(final CigarElement cigarElement) {
        final CigarOperator cigarOperator = cigarElement.getOperator();
        return cigarOperator == CigarOperator.I || cigarOperator == CigarOperator.S;
    }

    public RealignmentResult realign(final GATKRead read) {
        final List<BwaMemAlignment> alignments = aligner.alignSeqs(Arrays.asList(read), GATKRead::getBasesNoCopy).get(0);

        final List<BwaMemAlignment> nonAltAlignments = alignments.size() == 1 ? alignments :
                alignments.stream().filter(a -> a.getRefId() < numberOfRegularContigs).collect(Collectors.toList());
        return checkAlignments(nonAltAlignments, minAlignerScoreDifference);
    }

    @VisibleForTesting
    static RealignmentResult checkAlignments(final List<BwaMemAlignment> alignments, int minAlignerScoreDifference) {
        if (alignments.isEmpty()) {
            return new RealignmentResult(false, Collections.emptyList());
        } else if (alignments.get(0).getRefId() < 0) {
            return new RealignmentResult(false, alignments);
        } else if (alignments.size() == 1) {
            return new RealignmentResult(true, alignments);
        } else {
            final int scoreDifference = alignments.get(0).getAlignerScore() - alignments.get(1).getAlignerScore();
            return new RealignmentResult(scoreDifference >= minAlignerScoreDifference, alignments);
        }

        // TODO: we need to check that contig is the same or equivalent up to hg38 alt contig
        // TODO: do this by looking at index.getReferenceContigNames() and bamHeader.getSequenceDictionary().getSequences()
        // TODO: in IDE and seeing what the correspondence could be
        // TODO: put in check that start position is within eg 10 Mb of original mapping
    }

    // note that we realigned the original aligned bases of the read and mate, not the raw sequenced bases
    // thus one of them has been reverse-complemented and we look for pairs that map to the same strand
    public static List<Pair<BwaMemAlignment, BwaMemAlignment>> findPlausiblePairs(List<BwaMemAlignment> readRealignments, List<BwaMemAlignment> mateRealignments, int maxReasonableFragmentLength) {
        return readRealignments.stream()
                .flatMap(r -> mateRealignments.stream()
                        .filter(m -> m.getRefId() == r.getRefId())
                        .filter(m -> Math.abs(m.getRefStart() - r.getRefStart()) < maxReasonableFragmentLength)
                        .filter(m -> SAMFlag.READ_REVERSE_STRAND.isSet(m.getSamFlag()) == SAMFlag.READ_REVERSE_STRAND.isSet(r.getSamFlag()))
                        .map(m -> new Pair<>(r,m)))
                .collect(Collectors.toList());
    }

    public static class RealignmentResult {
        private final boolean mapsToSupposedLocation;
        private final List<BwaMemAlignment> realignments;

        public RealignmentResult(boolean mapsToSupposedLocation, List<BwaMemAlignment> realignments) {
            this.mapsToSupposedLocation = mapsToSupposedLocation;
            this.realignments = realignments;
        }

        public boolean isGood() { return mapsToSupposedLocation;  }

        public List<BwaMemAlignment> getRealignments() { return realignments; }
    }
}