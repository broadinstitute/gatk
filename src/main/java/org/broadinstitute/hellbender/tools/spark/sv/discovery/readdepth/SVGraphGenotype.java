package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import java.util.*;

public class SVGraphGenotype {

    protected final List<IndexedSVGraphPath> haplotypes;
    protected final int groupId;
    protected final int genotypeId;
    protected double depthLikelihood;
    protected double depthProbability;
    protected double evidenceProbability;
    protected double probability;
    private static final int DEFAULT_PLOIDY = 2;

    public SVGraphGenotype(final int groupId, final int genotypeId) {
        this.haplotypes = new ArrayList<>(DEFAULT_PLOIDY);
        this.groupId = groupId;
        this.genotypeId = genotypeId;
        this.depthLikelihood = Double.NaN;
        this.depthProbability = Double.NaN;
        this.evidenceProbability = Double.NaN;
        this.probability = Double.NaN;
    }

    public SVGraphGenotype(final int groupId, final int genotypeId, final Collection<IndexedSVGraphPath> pathList) {
        this(groupId, genotypeId);
        for (final IndexedSVGraphPath path : pathList) {
            addPath(path);
        }
    }

    public void addPath(final IndexedSVGraphPath path) {
        haplotypes.add(path);
        Collections.sort(haplotypes); //Maintain sorted haplotypes so genotypes can be easily deduplicated with a hash set
    }

    public double getDepthProbability() {
        return depthProbability;
    }

    public void setDepthProbability(double depthProbability) {
        this.depthProbability = depthProbability;
    }

    public double getEvidenceProbability() {
        return evidenceProbability;
    }

    public void setEvidenceProbability(double probability) {
        this.evidenceProbability = probability;
    }

    public double getProbability() {
        return probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

    public List<IndexedSVGraphPath> getHaplotypes() {
        return haplotypes;
    }

    public int getGroupId() {
        return groupId;
    }

    public int getGenotypeId() {
        return genotypeId;
    }

    public void setDepthLikelihood(double depthLikelihood) { this.depthLikelihood = depthLikelihood; }

    public double getDepthLikelihood() {
        return depthLikelihood;
    }

    /**
     * The genotype ID is not checked in equals() or hashCode()
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SVGraphGenotype)) return false;
        SVGraphGenotype that = (SVGraphGenotype) o;
        return groupId == that.groupId &&
                Double.compare(that.depthLikelihood, depthLikelihood) == 0 &&
                Double.compare(that.depthProbability, depthProbability) == 0 &&
                Double.compare(that.evidenceProbability, evidenceProbability) == 0 &&
                Double.compare(that.probability, probability) == 0 &&
                Objects.equals(haplotypes, that.haplotypes);
    }

    @Override
    public int hashCode() {

        return Objects.hash(haplotypes, groupId, depthLikelihood, depthProbability);
    }
}
