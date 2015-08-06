package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;

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
    protected GenotypeLikelihoodsCalculationModel(final UnifiedArgumentCollection UAC, final Logger logger) {
        if ( logger == null || UAC == null ) {
            throw new GATKException("Bad arguments");
        }
        this.UAC = UAC;
        this.logger = logger;
    }
}