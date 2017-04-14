package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;

/**
 * This class stores the data required for calculating the emission probability according to the
 * GATK4 coverage model:
 *
 * <p>
 *     Mean log multiplicative bias:
 *     {@link #mu} = m_t + (W.E[z_s])_t
 * </p>
 *
 * <p>
 *     Sum of target- and sample-specific unexplained variances about the mean log multiplicative bias:
 *     {@link #psi} = \Psi_t + \gamma_s
 * </p>
 *
 * <p>
 *     Raw read count:
 *     {@link #readCount} = n_{st}
 * </p>
 *
 * <p>
 *     Mapping error probability:
 *     {@link #mappingErrorProbability} = \varepsilon_{st}
 * </p>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelCopyRatioEmissionData implements Serializable {

    private static final long serialVersionUID = -7363264674200250712L;

    /**
     * Mean log multiplicative bias
     */
    private final double mu;

    /**
     * Total unexplained variance
     */
    private final double psi;

    /**
     * Raw read count
     */
    private final int readCount;

    /**
     * Mapping error probability
     */
    private final double mappingErrorProbability;

    /**
     * Sample metadata (sex genotype, depth of coverage, ...)
     *
     * @implNote The metadata is necessary for calculating the emission probability, though, they are the same
     * for all targets. We have separated such global parameters from target-dependent variables and wrapped them
     * in a class so that they can be "injected" by reference.
     */
    private CopyRatioCallingMetadata copyRatioCallingMetadata = null;

    /**
     * Public constructor.
     *
     * @param mu Mean log multiplicative bias
     * @param psi Total unexplained variance
     * @param readCount Raw read count
     * @param mappingErrorProbability Mapping error probability
     */
    public CoverageModelCopyRatioEmissionData(final double mu,
                                              final double psi,
                                              final int readCount,
                                              final double mappingErrorProbability) {
        this.mu = ParamUtils.isFinite(mu, "Log multiplicative bias must be finite. Bad value: " + mu);
        this.psi = ParamUtils.isPositiveOrZero(psi, "Unexplained variance must be a non-negative real number." +
                " Bad value: " + psi);
        this.readCount = ParamUtils.isPositiveOrZero(readCount, "Read count must be a non-negative real number." +
                " Bad value: " + readCount);
        this.mappingErrorProbability = ParamUtils.isPositiveOrZero(mappingErrorProbability, "Mapping error probability " +
                " must be non-negative. Bad value: " + mappingErrorProbability);
    }

    public double getMu() {
        return mu;
    }

    public double getPsi() {
        return psi;
    }

    public int getReadCount() {
        return readCount;
    }

    public double getMappingErrorProbability() {
        return mappingErrorProbability;
    }

    public CopyRatioCallingMetadata getCopyRatioCallingMetadata() {
        return Utils.nonNull(copyRatioCallingMetadata, "Sample metadata is queried but is not set");
    }

    /**
     * Set metadata.
     *
     * @param copyRatioCallingMetadata the metadata
     */
    public void setCopyRatioCallingMetadata(@Nonnull final CopyRatioCallingMetadata copyRatioCallingMetadata) {
        this.copyRatioCallingMetadata = Utils.nonNull(copyRatioCallingMetadata, "The metadata must be non-null");
        Utils.validateArg(copyRatioCallingMetadata.hasSampleCoverageDepth(), "The metadata must contain depth of coverage field");
        Utils.validateArg(copyRatioCallingMetadata.hasSampleSexGenotypeData(), "The metadata must contain sex genotype data");
    }

    @Override
    public String toString() {
        return "CoverageModelCopyRatioEmissionData{" +
                "mu=" + mu +
                ", psi=" + psi +
                ", readCount=" + readCount +
                ", mappingErrorProbability=" + mappingErrorProbability +
                ", copyRatioCallingMetadata=" + copyRatioCallingMetadata +
                '}';
    }
}
