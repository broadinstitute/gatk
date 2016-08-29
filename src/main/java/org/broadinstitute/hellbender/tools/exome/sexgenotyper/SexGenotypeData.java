package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

/**
 * A class for storing sex genotypes and their likelihoods.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SexGenotypeData {
    /* mandatory */
    private final String sampleName;
    private final String sexGenotype;
    private final boolean hasExtendedGenotypingInfo;

    /* optional */
    private final Map<String, Double> logLikelihoodMap;

    /**
     * Public constructor.
     *
     * @param sampleName string identifer of the sample
     * @param sexGenotype string identifier for the sex genotype of the sample
     * @param sexGenotypesList list of all sex genotypes
     * @param logLikelihoods log likelihood of all sex genotypes in the same order as {@param sexGenotypesList}
     */
    public SexGenotypeData(@Nonnull final String sampleName,
                           @Nonnull final String sexGenotype,
                           @Nullable final List<String> sexGenotypesList,
                           @Nullable final double[] logLikelihoods) {
        this.sampleName = sampleName;
        this.sexGenotype = sexGenotype;

        if (sexGenotypesList == null || logLikelihoods == null) {
            logLikelihoodMap = null;
            hasExtendedGenotypingInfo = false;
        } else {
            ParamUtils.isPositive(sexGenotypesList.size(), "List of sex genotypes must be non-empty");
            if (logLikelihoods.length != sexGenotypesList.size()) {
                throw new UserException.BadInput("The size of genotype likelihood array must be equal to the number of" +
                        " genotypes");
            }
            hasExtendedGenotypingInfo = true;
            logLikelihoodMap = new HashMap<>();
            IntStream.range(0, sexGenotypesList.size()).forEach(i -> logLikelihoodMap.put(sexGenotypesList.get(i), logLikelihoods[i]));
        }
    }

    /**
     * Returns the set of sex genotypes if available; otherwise, throws an {@link IllegalStateException}
     * @return set of genotype string identifiers
     */
    public Set<String> getSexGenotypesSet() {
        if (logLikelihoodMap != null) {
            return logLikelihoodMap.keySet();
        } else {
            throw new IllegalStateException("Genotyping inference data is not available");
        }
    }

    /**
     * Returns the log likelihood for a given genotype if available; otherwise, throws an {@link IllegalStateException}
     * @param genotype genotype string identifier
     * @return log likelihood
     */
    public double getLogLikelihoodPerGenotype(@Nonnull final String genotype) {
        if (logLikelihoodMap != null) {
            return logLikelihoodMap.get(genotype);
        } else {
            throw new IllegalStateException("Genotyping inference data is not available");
        }
    }

    /**
     * Returns true if genotyping inference data is available; false otherwise
     * @return boolean
     */
    public boolean hasExtendedGenotypingInfo() {
        return hasExtendedGenotypingInfo;
    }

    public String getSampleName() {
        return sampleName;
    }

    public String getSexGenotype() {
        return sexGenotype;
    }

    /**
     * Equality test is based on sample name and sex genotype; genotyping inference data will
     * be ignored in equality comparison.
     *
     * @param obj object to compare with
     * @return boolean
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof SexGenotypeData)) {
            return false;
        }

        final SexGenotypeData dat = (SexGenotypeData) obj;
        return sampleName.equals(dat.getSampleName()) && sexGenotype.equals(dat.getSexGenotype());
    }

    @Override
    public int hashCode() {
        return 31 * sampleName.hashCode() + sexGenotype.hashCode();
    }
}
