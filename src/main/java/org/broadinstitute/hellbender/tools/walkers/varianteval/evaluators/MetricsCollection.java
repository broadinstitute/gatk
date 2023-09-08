package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;

/**
 * Created by knoblett on 9/15/15.
 */

@Analysis(description = "Metrics Collection")
public class MetricsCollection extends VariantEvaluator {
    public MetricsCollection(VariantEvalEngine engine) {
        super(engine);
    }

    @DataPoint(description = "The concordance rate from CompOverlap", format = "%.2f")
    public double concordantRate;
    @DataPoint(description = "Number of SNPs from IndelSummary", format = "%d")
    public int nSNPs;
    @DataPoint(description = "Number of SNP loci from CountVariants", format = "%d")
    public long nSNPloci;
    @DataPoint(description = "Number of indels from IndelSummary", format = "%d")
    public int nIndels;
    @DataPoint(description = "Number of indel loci from MultiallelicSummary", format = "%d")
    public int nIndelLoci;
    @DataPoint(description = "Insertion  to deletion ratio from IndelSummary")
    public String indelRatio;
    @DataPoint(description = "Insertion to deletion ratio from CountVariants", format = "%.2f")
    public double indelRatioLociBased;
    @DataPoint(description = "The transition to transversion ratio from TiTvVariantEvaluator", format = "%.2f")
    public double tiTvRatio;

    public int getComparisonOrder() {return 2;}

    public void setData(double concordantRate, int nSNPs, long nSNPloci, int nIndels, int nIndelLoci, String indelRatio, double indelRatioLociBased, double tiTvRatio){
        this.concordantRate = concordantRate;
        this.nSNPs = nSNPs;
        this.nSNPloci = nSNPloci;
        this.nIndels = nIndels;
        this.nIndelLoci = nIndelLoci;
        this.indelRatio = indelRatio;
        this.indelRatioLociBased = indelRatioLociBased;
        this.tiTvRatio = tiTvRatio;
    }
}
