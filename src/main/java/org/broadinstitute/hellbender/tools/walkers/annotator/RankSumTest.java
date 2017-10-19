package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;


/**
 * Abstract root for all RankSum based annotations
 */
public abstract class RankSumTest extends InfoFieldAnnotation implements Annotation {
    protected static double INVALID_ELEMENT_FROM_READ = Double.NEGATIVE_INFINITY;
    private boolean useDithering = true;

    public RankSumTest(final boolean useDithering){
        this.useDithering = useDithering;
    }

    public RankSumTest(){
        this(false);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc, "vc is null");

        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.isEmpty()) {
            return Collections.emptyMap();
        }

        final List<Double> refQuals = new ArrayList<>();
        final List<Double> altQuals = new ArrayList<>();

        final int refLoc = vc.getStart();

        if( likelihoods != null) {
            if (likelihoods.hasFilledLikelihoods()) {
                // Default to using the likelihoods to calculate the rank sum
                fillQualsFromLiklihood(vc, likelihoods, refQuals, altQuals, refLoc);

            // Use the pileup to stratify otherwise
            } else {
                for (final PileupElement p : likelihoods.getStratifiedPileups(vc).values().stream().flatMap(Collection::stream).collect(Collectors.toList())) {
                    final OptionalDouble value = getElementForPileupElement(p, refLoc);
                    if (value.isPresent() && value.getAsDouble() != INVALID_ELEMENT_FROM_READ && isUsableBase(p)) {
                        if (vc.getReference().equals(Allele.create(p.getBase(), true))) {
                            refQuals.add(value.getAsDouble());
                        } else if (vc.hasAllele(Allele.create(p.getBase()))) {
                            altQuals.add(value.getAsDouble());
                        }
                    }
                }
            }
        }


        if ( refQuals.isEmpty() && altQuals.isEmpty() ) {
            return Collections.emptyMap();
        }

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final MannWhitneyU.Result result = mannWhitneyU.test(Doubles.toArray(altQuals), Doubles.toArray(refQuals), MannWhitneyU.TestType.FIRST_DOMINATES);
        final double zScore = result.getZ();

        if (Double.isNaN(zScore)) {
            return Collections.emptyMap();
        } else {
            return Collections.singletonMap(getKeyNames().get(0), String.format("%.3f", zScore));
        }
    }

    protected void fillQualsFromLiklihood(VariantContext vc, ReadLikelihoods<Allele> likelihoods, List<Double> refQuals, List<Double> altQuals, int refLoc) {
        for (final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles()) {
            final GATKRead read = bestAllele.read;
            final Allele allele = bestAllele.allele;
            if (bestAllele.isInformative() && isUsableRead(read, refLoc)) {
                final OptionalDouble value = getElementForRead(read, refLoc, bestAllele);
                // Bypass read if the clipping goal is not reached or the refloc is inside a spanning deletion
                if (value.isPresent() && value.getAsDouble() != INVALID_ELEMENT_FROM_READ) {
                    if (allele.isReference()) {
                        refQuals.add(value.getAsDouble());
                    } else if (vc.hasAllele(allele)) {
                        altQuals.add(value.getAsDouble());
                    }
                }
            }
        }
    }

    /**
     * Get the element for the given read at the given reference position
     *
     * @param read     the read
     * @param refLoc   the reference position
     * @param bestAllele the most likely allele for this read
     * @return a Double representing the element to be used in the rank sum test, or null if it should not be used
     */
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc, ReadLikelihoods<Allele>.BestAllele bestAllele) {
        return getElementForRead(read, refLoc);
    }

    /**
     * Get the element for the given read at the given reference position
     *
     * @param read     the read
     * @param refLoc   the reference position
     * @return an OptionalDouble representing the element to be used in the rank sum test, empty if it should not be used
     */
    protected abstract OptionalDouble getElementForRead(final GATKRead read, final int refLoc);

    /**
     * Can the read be used in comparative tests between ref / alt bases?
     *
     * @param read   the read to consider
     * @param refLoc the reference location
     * @return true if this read is meaningful for comparison, false otherwise
     */
    protected boolean isUsableRead(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    /**
     * Can the base in this pileup element be used in comparative tests between ref / alt bases?
     *
     * Note that this function by default does not allow deletion pileup elements
     *
     * @param p the pileup element to consider
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    protected boolean isUsableBase(final PileupElement p) {
        return !(p.isDeletion() ||
                p.getMappingQual() == 0 ||
                p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                ((int) p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE); // need the unBAQed quality score here
    }

    /**
     * Get the element for the given read at the given reference position
     *
     * By default this function returns null, indicating that the test doesn't support the old style of pileup calculations
     *
     * @param p        the pileup element
     * @return a Double representing the element to be used in the rank sum test, or null if it should not be used
     */
    protected OptionalDouble getElementForPileupElement(final PileupElement p, final int refLoc) {
        // default to returning the same value
        return getElementForRead(p.getRead(), refLoc);
    }

}
