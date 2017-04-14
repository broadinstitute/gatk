package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hdf5.Utils;
import org.broadinstitute.hellbender.tools.exome.sexgenotyper.SexGenotypeData;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;

/**
 * This class stores a number of global sample-specific metadata such as sample sex genotype, sample name,
 * and sample index, and sample coverage depth.
 *
 * The metadata is constructed via a builder pattern, and is intended to be consumed by
 * {@link org.broadinstitute.hellbender.tools.coveragemodel.interfaces.CopyRatioExpectationsCalculator}
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CopyRatioCallingMetadata implements Serializable {

    private static final long serialVersionUID = -2933334495362553231L;

    /**
     * Various emission probability calculation models
     */
    public enum EmissionCalculationStrategy {
        /**
         * Fully Poisson; this mode is experimental since this model is incompatible with the EM
         * algorithm
         */
        POISSON,

        /**
         * Hybrid Poisson and Laplace+Poisson; this must be used for final posterior calling
         * (i.e. after learning). The Poisson approximation is used on low coverage targets
         * and Laplace approximation on Poisson is used on high coverage targets
         */
        HYBRID_POISSON_GAUSSIAN
    }

    private String sampleName = null;
    private Integer sampleIndex = null;
    private SexGenotypeData sampleSexGenotypeData = null;
    private Double sampleCoverageDepth = null;
    private EmissionCalculationStrategy emissionCalculationStrategy = null;
    private Double sampleAverageMappingErrorProbability = null;

    private CopyRatioCallingMetadata() {
    }

    private void setSampleName(@Nonnull final String sampleName) {
        this.sampleName = Utils.nonNull(sampleName, "Sample name must be non-null");
    }

    private void setSampleIndex(@Nonnull final Integer sampleIndex) {
        this.sampleIndex = Utils.nonNull(sampleIndex, "Sample index must be non-null");
    }

    private void setSampleSexGenotypeData(@Nonnull final SexGenotypeData sexGenotypeData) {
        this.sampleSexGenotypeData = Utils.nonNull(sexGenotypeData, "Sample sex genotype data must be non-null");
    }

    private void setSampleCoverageDepth(@Nonnull final Double sampleCoverageDepth) {
        this.sampleCoverageDepth = ParamUtils.isPositiveOrZero(
                Utils.nonNull(sampleCoverageDepth, "Sample coverage depth must be non-null"),
                "Sample coverage depth must be non-negative");
    }

    private void setEmissionCalculationStrategy(@Nonnull final EmissionCalculationStrategy emissionCalculationStrategy) {
        this.emissionCalculationStrategy = Utils.nonNull(emissionCalculationStrategy,
                "Emission calculation model must be non-null");
    }

    private void setSampleAverageMappingErrorProbability(@Nonnull final Double sampleAverageMappingErrorProbability) {
        this.sampleAverageMappingErrorProbability = ParamUtils.isPositiveOrZero(
                Utils.nonNull(sampleAverageMappingErrorProbability, "Sample average mapping error probability must be non-null"),
                "Sample average mapping error probability must be non-negative");
    }

    public String getSampleName() {
        return Utils.nonNull(sampleName, "Sample name is queried but is not set");
    }

    public Integer getSampleIndex() {
        return Utils.nonNull(sampleIndex, "Sample index is queried but is not set");
    }

    public SexGenotypeData getSampleSexGenotypeData() {
        return Utils.nonNull(sampleSexGenotypeData, "Sample sex genotype data is queried but is not set");
    }

    public Double getSampleCoverageDepth() {
        return Utils.nonNull(sampleCoverageDepth, "Sample coverage depth is queried but is not set");
    }

    public Double getSampleAverageMappingErrorProbability() {
        return Utils.nonNull(sampleAverageMappingErrorProbability, "Sample average mapping error probability is" +
                " queried but is not set");
    }

    public EmissionCalculationStrategy getEmissionCalculationStrategy() {
        return Utils.nonNull(emissionCalculationStrategy, "Emission calculation strategy is queried but is not set");
    }


    public boolean hasSampleName() {
        return sampleName != null;
    }

    public boolean hasSampleIndex() {
        return sampleIndex != null;
    }

    public boolean hasSampleSexGenotypeData() {
        return sampleSexGenotypeData != null;
    }

    public boolean hasSampleCoverageDepth() {
        return sampleCoverageDepth != null;
    }

    public boolean hasEmissionCalculationStrategy() {
        return emissionCalculationStrategy != null;
    }

    public static SampleMetadataBuilder builder() {
        return new SampleMetadataBuilder();
    }


    public static class SampleMetadataBuilder implements Serializable {

        private static final long serialVersionUID = 4431390529312114201L;

        final CopyRatioCallingMetadata metadata;

        SampleMetadataBuilder() {
            this.metadata = new CopyRatioCallingMetadata();
        }

        public SampleMetadataBuilder setSampleName(@Nonnull final String sampleName) {
            metadata.setSampleName(sampleName);
            return this;
        }

        public SampleMetadataBuilder setSampleIndex(@Nonnull final Integer sampleIndex) {
            metadata.setSampleIndex(sampleIndex);
            return this;
        }

        public SampleMetadataBuilder setSampleSexGenotypeData(@Nonnull final SexGenotypeData sexGenotypeData) {
            metadata.setSampleSexGenotypeData(sexGenotypeData);
            return this;
        }

        public SampleMetadataBuilder setSampleCoverageDepth(@Nonnull final Double sampleCoverageDepth) {
            metadata.setSampleCoverageDepth(sampleCoverageDepth);
            return this;
        }

        public SampleMetadataBuilder setSampleAverageMappingErrorProbability(@Nonnull final Double sampleAverageMappingErrorProbability) {
            metadata.setSampleAverageMappingErrorProbability(sampleAverageMappingErrorProbability);
            return this;
        }

        public SampleMetadataBuilder setEmissionCalculationStrategy(@Nonnull final EmissionCalculationStrategy emissionCalculationStrategy) {
            metadata.setEmissionCalculationStrategy(emissionCalculationStrategy);
            return this;
        }

        public CopyRatioCallingMetadata build() {
            return metadata;
        }
    }

    @Override
    public String toString() {
        return "CopyRatioCallingMetadata{" +
                "sampleName='" + sampleName + '\'' +
                ", sampleIndex=" + sampleIndex +
                ", sampleSexGenotypeData=" + sampleSexGenotypeData +
                ", sampleCoverageDepth=" + sampleCoverageDepth +
                ", emissionCalculationStrategy=" + emissionCalculationStrategy +
                ", sampleAverageMappingErrorProbability=" + sampleAverageMappingErrorProbability +
                '}';
    }
}
