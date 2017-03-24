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
 * Created by David Benjamin on 3/20/17.
 */
public class BaseQuality extends InfoFieldAnnotation {
    public static final String REFERENCE_MEDIAN_QUALITY_KEY = "REF_MED_QUAL";
    public static final String ALT_MEDIAN_QUALITY_KEY = "ALT_MED_QUAL";

    @Override
    public List<String> getKeyNames() { return Arrays.asList(REFERENCE_MEDIAN_QUALITY_KEY, ALT_MEDIAN_QUALITY_KEY); }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {

        final List<Double> refQuals = new ArrayList<>();
        final List<Double> altQuals = new ArrayList<>();

        final int refLoc = vc.getStart();

        for (final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles()) {
            final GATKRead read = bestAllele.read;
            final Allele allele = bestAllele.allele;
            if (bestAllele.isInformative() && isUsableRead(read, refLoc)) {
                final OptionalDouble value = getQual(read, refLoc);
                // Bypass read if the clipping goal is not reached or the refloc is inside a spanning deletion
                if ( value.isPresent() ) {
                    (allele.isReference() ? refQuals : altQuals).add(value.getAsDouble());
                }
            }
        }

        final double refMedian = refQuals.size() == 0 ? -1 : new Median().evaluate(Doubles.toArray(refQuals));
        final double altMedian = altQuals.size() == 0 ? -1 : new Median().evaluate(Doubles.toArray(altQuals));

        return ImmutableMap.of(REFERENCE_MEDIAN_QUALITY_KEY, (int) refMedian, ALT_MEDIAN_QUALITY_KEY, (int) altMedian);
    }

    private OptionalDouble getQual(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);

        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        return offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), offset) ?
                OptionalDouble.empty() : OptionalDouble.of(read.getBaseQuality(offset));
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(REFERENCE_MEDIAN_QUALITY_KEY, 1, VCFHeaderLineType.Integer,
                        "Median distance of reference alleles to end of read."),
                new VCFInfoHeaderLine(ALT_MEDIAN_QUALITY_KEY, 1, VCFHeaderLineType.Integer,
                        "Median distance of alt alleles to end of read."));
    }

    private boolean isUsableRead(final GATKRead read, final int refLoc) {
        return read.getMappingQuality() != 0 && !read.isUnmapped() && ReadUtils.getSoftEnd(read) >= refLoc;
    }
}