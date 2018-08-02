package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Trims (hard clips) soft-clipped bases due to the following artifact:
 *
 * When a sequence and its reverse complement occur near opposite ends of a fragment DNA damage (especially in the case
 * of FFPE samples and ancient DNA) can disrupt base-pairing causing a single-strand loop of the sequence and its reverse
 * complement, after which end repair copies the true 5' end of the fragment onto the 3' end of the fragment.  That is, the
 * artifact looks like this (A' denotes the reverse complement of A)
 *
 * Biological sequence:
 * Forward strand 3' A  B  . . . B' C 5'
 * Reverse strand 5' A' B' . . . B  C' 3'
 *
 * Loop structure (B/B', C/C' between forward and strands are *not* hybridized due to damage)
 * Reverse strand 3' C' B   . . . . . .
 * Forward strand 5' C  B'  . . . . . .
 *                      |           . .
 * Forward strand 3'    B   . . . . . .
 * Reverse strand 5'    B'  . . . . . .
 *
 * After end repair of 5' overhang of sequence C on self-looped forward strand
 * Reverse strand 3' C' B   . . . . . .
 * Forward strand 5' C  B'  . . . . . .
 *                   |  |           . .
 * Forward strand 3' C' B   . . . . . .
 * Reverse strand 5'    B'  . . . . . .
 *
 * Forward strand with artifact after denaturing: (C' replaces A)
 * Forward strand 3' C' B . . . B' C
 *
 * Since there is no good way early in GATK tools to collect reads and mates, here we filter only for the case where
 * sequence C matches the reference, so that the reference sequence is a proxy for the 5' end of the mate.  This catches
 * most errors and saves a lot of runtime by reducing false active regions and by simplifying the assembly.
 */
public final class PalindromeArtifactClipReadTransformer implements ReadTransformer {

    private static final long serialVersionUID = 1L;

    public static final double MIN_FRACTION_OF_MATCHING_BASES = 0.9;

    private final ReferenceDataSource referenceDataSource;

    // the artifact requires a minimum number of palindromic bases to create the single-strand loop
    private final int minPalindromeSize;

    public PalindromeArtifactClipReadTransformer(final ReferenceDataSource referenceDataSource, final int minPalindromeSize) {
        this.referenceDataSource = referenceDataSource;
        this.minPalindromeSize = minPalindromeSize;
    }

    @Override
    public GATKRead apply(final GATKRead read) {
        final int adaptorBoundary = read.getAdaptorBoundary();
        if (!read.isProperlyPaired() || adaptorBoundary == ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY) {
            return read;
        }

        final Cigar cigar = read.getCigar();
        final CigarOperator firstOperator = cigar.getFirstCigarElement().getOperator();
        final CigarOperator lastOperator = cigar.getLastCigarElement().getOperator();
        final boolean readIsUpstreamOfMate = read.getFragmentLength() > 0;

        if ((readIsUpstreamOfMate && firstOperator != CigarOperator.SOFT_CLIP && firstOperator != CigarOperator.INSERTION) ||
                (!readIsUpstreamOfMate && lastOperator != CigarOperator.SOFT_CLIP && lastOperator != CigarOperator.INSERTION)) {
            return read;
        }

        final int potentialArtifactBaseCount = readIsUpstreamOfMate ? cigar.getFirstCigarElement().getLength() :
                cigar.getLastCigarElement().getLength();


        final int numBasesToCompare = Math.min(potentialArtifactBaseCount + minPalindromeSize, read.getLength());

        // the reference position of bases that are the reverse complement of the suspected artifact
        final int refStart = readIsUpstreamOfMate ? adaptorBoundary - numBasesToCompare : adaptorBoundary + 1;
        final int refEnd = readIsUpstreamOfMate ? adaptorBoundary - 1 : adaptorBoundary + numBasesToCompare;

        // if the reference bases overlap the soft-clipped bases, it's not an artifact.  Note that read.getStart() is the unclipped start
        // this can happen, for a read with a huge soft-clip that overlaps its mate significantly.
        if ( (readIsUpstreamOfMate && refStart < read.getStart()) || (!readIsUpstreamOfMate && read.getEnd() < refEnd)) {
            return read;
        }
        final SimpleInterval refInterval = new SimpleInterval(read.getContig(), refStart, refEnd);
        final MutableInt numMatch = new MutableInt(0);

        // we traverse the reference in the forward direction, hence the read in the reverse direction
        final MutableInt readIndex = new MutableInt(readIsUpstreamOfMate ? numBasesToCompare - 1 : read.getLength() - 1);
        referenceDataSource.query(refInterval).forEachRemaining(refBase -> {
            if (BaseUtils.getComplement(refBase) == read.getBase(readIndex.getValue())) {
                numMatch.increment();
            }
            readIndex.decrement();
        });

        if (numMatch.doubleValue() / numBasesToCompare >= MIN_FRACTION_OF_MATCHING_BASES) {
            final ReadClipper readClipper = new ReadClipper(read);
            final ClippingOp clippingOp = readIsUpstreamOfMate ? new ClippingOp(0, potentialArtifactBaseCount - 1) :
                    new ClippingOp(read.getLength() - potentialArtifactBaseCount, read.getLength());
            readClipper.addOp(clippingOp);
            return readClipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);
        } else {
            return read;
        }
    }
}

