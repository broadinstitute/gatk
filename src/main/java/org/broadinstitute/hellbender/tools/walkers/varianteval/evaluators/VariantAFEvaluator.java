package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

@Analysis(description = "Computes different estimates of theta based on variant sites and genotypes")

public class VariantAFEvaluator extends VariantEvaluator {
    @DataPoint(description = "Average variant allele fraction over all variant sites", format = "%.8f")
    public double avgVarAF = 0.0;
    @DataPoint(description = "Number of called sites over all variant sites;", format = "%d")
    public int totalCalledSites;
    @DataPoint(description = "Number of called heterozygous sites;", format = "%d")
    public int totalHetSites;
    @DataPoint(description = "Number of called homozygous variant sites;", format = "%d")
    public int totalHomVarSites;
    @DataPoint(description = "Number of called homozygous reference sites;", format = "%d")
    public int totalHomRefSites;

    private static final double PLOIDY = 2.0; // assume ploidy of 2
    private double sumVariantAFs; // this is the allele fraction we're summing over all sites, to be used to calculate the avgVarAlleles

    public VariantAFEvaluator(VariantEvalEngine engine) {
        super(engine);
        sumVariantAFs = 0.0;
        totalCalledSites = 0;
    }

    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void update1(final VariantContext vc, final VariantEvalContext context) {
        vc.getStart();

        if (vc == null || !vc.isSNP() || (getEngine().getVariantEvalArgs().ignoreAC0Sites() && vc.isMonomorphicInSamples())) {
            return;
        }

        for (final Genotype genotype : vc.getGenotypes()) {
             // eval array
            if (!genotype.isNoCall()) {
                if (genotype.getPloidy() != PLOIDY) {

                    throw new UserException.BadInput("This tool only works with ploidy 2");
                }
                // add AF at this site
                this.totalCalledSites += 1;
                int numReferenceAlleles= genotype.countAllele(vc.getReference());
                double varAFHere = (PLOIDY - numReferenceAlleles)/PLOIDY;
                this.sumVariantAFs += varAFHere;

                totalHetSites += numReferenceAlleles == 1 ? 1 : 0;
                totalHomVarSites += numReferenceAlleles == 0 ? 1 : 0;
                totalHomRefSites += numReferenceAlleles == 2 ? 1 : 0;

            }
        }

        if (!vc.hasGenotypes()) {
            // comp  ( sites only thousand genomes )
            this.totalCalledSites += 1;
            this.sumVariantAFs += vc.getAttributeAsDouble("AF", 0.0);
        }
    }

    @Override
    public void finalizeEvaluation() {
        this.avgVarAF = this.totalCalledSites == 0 ? 0 : this.sumVariantAFs / this.totalCalledSites;
    }
}

