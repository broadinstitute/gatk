package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;

@Analysis(description = "Computes different estimates of theta based on variant sites and genotypes")

public class PVariantEvaluator extends VariantEvaluator {
    @DataPoint(description = "Average number of variant alleles at variant sites; note that missing genotypes are ignored when computing this value", format = "%.8f")
    public double avgVarAlleles= 0.0;
    @DataPoint(description = "Sum of variant alleles over all variant sites; divide this by total target to get estimate of per base p", format = "%.8f")
    public double totalVarAlleles;

    @DataPoint(description = "Number of called alleles over all variant sites;", format = "%d")
    public int totalCalledAlleles;

    @DataPoint(description = "Number of heterozygous sites over all variant sites;", format = "%d")
    public int totalHetSites;

    @DataPoint(description = "Number of homozygous reference sites over all variant sites;", format = "%d")
    public int totalHomRefSites;

    @DataPoint(description = "Number of homozygous variant sites over all variant sites;", format = "%d")
    public int totalHomVarSites;

    private static final int RANDOM_SEED = 51; // used in theoretical sensitivity
    private static final RandomGenerator rg = new Well19937c(RANDOM_SEED);


    public PVariantEvaluator() {
        totalVarAlleles = 0.0;
    }

    public int getComparisonOrder() {
        return 1;
    }

    public void update1(VariantContext vc, final ReferenceContext referenceContext, final ReadsContext readsContext, final FeatureContext featureContext) {
        vc.getStart();

        if (vc == null || !vc.isSNP() || (getWalker().ignoreAC0Sites() && vc.isMonomorphicInSamples())) {
            return;
        }

        int numGenosHere = 0;
        int polymorphicChroms = 0;

        for (final Genotype genotype : vc.getGenotypes()) {
            if (!genotype.isNoCall()) {
                numGenosHere+=genotype.getPloidy();

                polymorphicChroms = genotype.getPloidy() - genotype.countAllele(vc.getReference());

                if (genotype.isHomRef()) {
                    this.totalHomRefSites += 1;
                } else if (genotype.isHet()) {
                    this.totalHetSites += 1;
                } else {
                    // het var genotypes will also count as hom var
                    this.totalHomVarSites+= 1;
                }
            }
        }

        if (!vc.hasGenotypes()){

            // try to get var allele counts using AF as a binomial

            double p = vc.getAttributeAsDouble("AF", 0.0);

            final BinomialDistribution bd;
            synchronized (rg) {
                bd = new BinomialDistribution(rg, 2, p);
            }

            int num_var_alleles = bd.sample();
            this.totalVarAlleles += num_var_alleles;
            this.totalCalledAlleles += vc.hasAttribute("AF") ? 2 : 0;
            if (num_var_alleles == 2) {
                this.totalHomVarSites += 1;
            } else {
                if (num_var_alleles == 1) {
                    this.totalHetSites += 1;
                } else {
                    this.totalHomRefSites += 1;
                }
            }

        }

        if (numGenosHere > 0) {
            //only if have one called genotype at least
            this.totalCalledAlleles += numGenosHere;
            this.totalVarAlleles += polymorphicChroms;

        }
    }

    @Override
    public void finalizeEvaluation() {

        if (this.totalCalledAlleles > 0) {
            this.avgVarAlleles = this.totalVarAlleles / totalCalledAlleles;

        }
    }
}