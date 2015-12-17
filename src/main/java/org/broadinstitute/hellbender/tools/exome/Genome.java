package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;

/**
 * Represents copy-number data from a single individual for exome analysis.  Contains coverage data from targets
 * and ref/alt allele counts at normal germline het SNP sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Genome {
    private final TargetCollection<TargetCoverage> targets;
    private final TargetCollection<AllelicCount> snps;
    private final String sampleName;

    /**
     * Constructs a genome from lists containing linear target-coverage and SNP-allele-count data.
     * @param targets       list of linear target coverages, cannot be {@code null}
     * @param snps          list of SNP allele counts, cannot be {@code null}
     * @param sampleName    name of the sample, cannot be {@code null}
     */
    public Genome(final List<TargetCoverage> targets, final List<AllelicCount> snps, final String sampleName) {
        Utils.nonNull(targets, "The list of targets cannot be null.");
        Utils.nonNull(snps, "The list of SNPs cannot be null.");
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        this.targets = new HashedListTargetCollection<>(targets);
        this.snps = new HashedListTargetCollection<>(snps);
        this.sampleName = sampleName;
    }

    /**
     * Constructs a genome from files containing linear target-coverage and SNP-allele-count data.
     * @param targetFile    linear target-coverage file
     * @param snpFile       SNP-allele-count file
     * @param sampleName    name of the sample
     */
    public Genome(final File targetFile, final File snpFile, final String sampleName) {
        this(TargetCoverageUtils.readTargetsWithCoverage(targetFile), new AllelicCountCollection(snpFile).getCounts(),
                sampleName);
    }

    public final TargetCollection<TargetCoverage> getTargets() {  return targets; }

    public final TargetCollection<AllelicCount> getSNPs() {  return snps; }

    public final String getSampleName() {   return sampleName;  }
}
