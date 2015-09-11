package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;


/**
 * Abstract root for all RankSum based annotations
 */
public abstract class RankSumTest extends InfoFieldAnnotation {
    private boolean useDithering = true;

    public RankSumTest(final boolean useDithering){
        this.useDithering = useDithering;
    }

    public RankSumTest(){
        this(true);
    }

    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(stratifiedPerReadAlleleLikelihoodMap, "stratifiedPerReadAlleleLikelihoodMap has to be non-null");
        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.isEmpty()) {
            return null;
        }

        final List<Double> refQuals = new ArrayList<>();
        final List<Double> altQuals = new ArrayList<>();

        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
                final PerReadAlleleLikelihoodMap likelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());
                if ( likelihoodMap != null && !likelihoodMap.isEmpty() ) {
                    fillQualsFromLikelihoodMap(vc.getAlleles(), vc.getStart(), likelihoodMap, refQuals, altQuals);
                }
        }

        if ( refQuals.isEmpty() && altQuals.isEmpty() ) {
            return null;
        }

        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final double p = MannWhitneyU.runOneSidedTest(useDithering, altQuals, refQuals).getLeft();
        if (Double.isNaN(p)) {
            return Collections.emptyMap();
        } else {
            return Collections.singletonMap(getKeyNames().get(0), String.format("%.3f", p));
        }
    }

    private void fillQualsFromLikelihoodMap(final List<Allele> alleles,
                                            final int refLoc,
                                            final PerReadAlleleLikelihoodMap likelihoodMap,
                                            final List<Double> refQuals,
                                            final List<Double> altQuals) {
        for ( final Map.Entry<GATKRead, Map<Allele,Double>> el : likelihoodMap.getLikelihoodReadMap().entrySet() ) {
            final MostLikelyAllele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if ( ! a.isInformative() ) {
                continue; // read is non-informative
            }

            final GATKRead read = el.getKey();
            if ( isUsableRead(read, refLoc) ) {
                final OptionalDouble value = getElementForRead(read, refLoc, a);
                if ( !value.isPresent() ) {
                    continue;
                }

                if ( a.getMostLikelyAllele().isReference() ) {
                    refQuals.add(value.getAsDouble());
                } else if ( alleles.contains(a.getMostLikelyAllele()) ) {
                    altQuals.add(value.getAsDouble());
                }
            }
        }
    }

    /**
     * Get the element for the given read at the given reference position
     *
     * @param read     the read
     * @param refLoc   the reference position
     * @param mostLikelyAllele the most likely allele for this read
     * @return a Double representing the element to be used in the rank sum test, or null if it should not be used
     */
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc, final MostLikelyAllele mostLikelyAllele) {
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
        return !( read.getMappingQuality() == 0 ||
                read.getMappingQuality() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE );
    }
}