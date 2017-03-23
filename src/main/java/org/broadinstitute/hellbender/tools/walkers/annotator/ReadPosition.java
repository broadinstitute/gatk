package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Created by davidben on 3/20/17.
 */
public class ReadPosition extends InfoFieldAnnotation {
    public static final String REFERENCE_MEDIAN_POSITION_KEY = "REF_MED_POS";
    public static final String ALT_MEDIAN_POSITION_KEY = "ALT_MED_POS";

    @Override
    public List<String> getKeyNames() { return Arrays.asList(REFERENCE_MEDIAN_POSITION_KEY, ALT_MEDIAN_POSITION_KEY); }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {

        final List<Double> refPositions = new ArrayList<>();
        final List<Double> altPositions = new ArrayList<>();

        final int refLoc = vc.getStart();

        for (final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles()) {
            final GATKRead read = bestAllele.read;
            final Allele allele = bestAllele.allele;
            if (bestAllele.isInformative() && isUsableRead(read, refLoc)) {
                final OptionalDouble value = getDistanceFromEnd(read, refLoc);
                // Bypass read if the clipping goal is not reached or the refloc is inside a spanning deletion
                if ( value.isPresent() ) {
                    (allele.isReference() ? refPositions : altPositions).add(value.getAsDouble());
                }
            }
        }

        final double refMedian = refPositions.size() == 0 ? -1 : new Median().evaluate(Doubles.toArray(refPositions));
        final double altMedian = altPositions.size() == 0 ? -1 : new Median().evaluate(Doubles.toArray(altPositions));

        return ImmutableMap.of(REFERENCE_MEDIAN_POSITION_KEY, (int) refMedian, ALT_MEDIAN_POSITION_KEY, (int) altMedian);
    }

    private OptionalDouble getDistanceFromEnd(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), offset)) {
            return OptionalDouble.empty();
        }

        int readPosition = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), offset, false, 0, 0);
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
        final int distanceFromEnd = Math.min(readPosition, numAlignedBases - readPosition - 1);
        return OptionalDouble.of(distanceFromEnd);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(REFERENCE_MEDIAN_POSITION_KEY, 1, VCFHeaderLineType.Integer,
                        "Median distance of reference alleles to end of read."),
                new VCFInfoHeaderLine(ALT_MEDIAN_POSITION_KEY, 1, VCFHeaderLineType.Integer,
                        "Median distance of alt alleles to end of read."));
    }

    private boolean isUsableRead(final GATKRead read, final int refLoc) {
        return read.getMappingQuality() != 0 && !read.isUnmapped() && ReadUtils.getSoftEnd(read) >= refLoc;
    }
}