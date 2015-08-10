package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;

/**
 * Represents copy-number data from a single individual for exome analysis.  Contains interval and copy ratio
 * data from targets (units of coverage-data collection) and ref/alt allele counts at normal germline het SNP sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Genome {
    private final TargetCollection<TargetCoverage> targets;
    private final TargetCollection<AllelicCount> snps;
    private final String sampleName;

    /**
     * Constructs a genome from files containing target-coverage and SNP allele-count data.
     * @param targetFile    target-coverage file
     * @param snpFile       SNP allele-count file
     * @param sampleName    name of the sample
     */
    public Genome(final File targetFile, final File snpFile, final String sampleName) {
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        final List<TargetCoverage> targets = TargetCoverageUtils.readTargetsWithCoverage(targetFile);
        final List<AllelicCount> snps = new AllelicCountCollection(snpFile).getCounts();
        this.targets = new HashedListTargetCollection<>(targets);
        this.snps = new HashedListTargetCollection<>(snps);
        this.sampleName = sampleName;
    }

    public final TargetCollection<TargetCoverage> getTargets() {  return targets; }

    public final TargetCollection<AllelicCount> getSNPs() {  return snps; }
}
