package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.List;
import java.util.Map;

/**
 * The model representing how we calculate genotype likelihoods
 */
public abstract class GenotypeLikelihoodsCalculationModel {

    public static final String DUMMY_LANE = "Lane1";
    public static final String DUMMY_SAMPLE_NAME = "DummySample1";

    public enum Model {
        SNP,
        INDEL,
        GENERALPLOIDYSNP,
        GENERALPLOIDYINDEL,
        BOTH;
    }

    protected final UnifiedArgumentCollection UAC;
    protected Logger logger;

    /**
     * Create a new object
     * @param logger        logger
     * @param UAC           unified arg collection
     */
    protected GenotypeLikelihoodsCalculationModel(UnifiedArgumentCollection UAC, Logger logger) {
        Utils.nonNull(UAC);
        Utils.nonNull(logger);
        this.UAC = UAC;
        this.logger = logger;
    }

    /**
     * Can be overridden by concrete subclasses
     *
     * @param features              feature data
     * @param ref                   reference context
     * @param contexts              stratified alignment contexts
     * @param contextType           stratified context type
     * @param allAllelesToUse       the alternate allele to use, null if not set
     * @param useBAQedPileup        should we use the BAQed pileup or the raw one?
     * @return variant context where genotypes are no-called but with GLs
     */
    public abstract VariantContext getLikelihoods(final FeatureContext features,
                                                  final ReferenceContext ref,
                                                  final Map<String, AlignmentContext> contexts,
                                                  final AlignmentContext.ReadOrientation contextType,
                                                  final List<Allele> allAllelesToUse,
                                                  final boolean useBAQedPileup,
                                                  final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap);


    protected int getFilteredDepth(ReadPileup pileup) {
        int count = 0;
        for ( PileupElement p : pileup ) {
            if ( BaseUtils.isRegularBase(p.getBase()) ) {
                count++;
            }
        }

        return count;
    }

}
