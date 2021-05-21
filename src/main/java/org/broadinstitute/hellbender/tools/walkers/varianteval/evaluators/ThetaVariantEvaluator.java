package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

@Analysis(description = "Computes different estimates of theta based on variant sites and genotypes")
public class ThetaVariantEvaluator extends VariantEvaluator {
    @DataPoint(description = "Average heterozygosity at variant sites; note that missing genotypes are ignored when computing this value", format = "%.8f")
    public double avgHet = 0.0;
    @DataPoint(description = "Average pairwise differences at aligned sequences; averaged over both number of sequeneces and number of variant sites; note that missing genotypes are ignored when computing this value", format = "%.8f")
    public double avgAvgDiffs = 0.0;
    @DataPoint(description = "Sum of heterozygosity over all variant sites; divide this by total target to get estimate of per base theta", format = "%.8f")
    public double totalHet = 0.0;
    @DataPoint(description = "Sum of pairwise diffs over all variant sites; divide this by total target to get estimate of per base theta", format = "%.8f")
    public double totalAvgDiffs = 0.0;
    @DataPoint(description = "Theta for entire region estimated based on number of segregating sites; divide ths by total target to get estimate of per base theta", format = "%.8f")
    public double thetaRegionNumSites = 0.0;

    //helper variables
    double numSites = 0;

    public ThetaVariantEvaluator(VariantEvalEngine engine) {
        super(engine);
    }

    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void update1(final VariantContext vc, final VariantEvalContext context) {
        if (vc == null || !vc.isSNP() || (getEngine().getVariantEvalArgs().ignoreAC0Sites() && vc.isMonomorphicInSamples())) {
            return;
        }

        //this maps allele to a count
        ConcurrentMap<String, Integer> alleleCounts = new ConcurrentHashMap<String, Integer>();

        int numHetsHere = 0;
        int numGenosHere = 0;
        int numIndsHere = 0;

        for (final Genotype genotype : vc.getGenotypes()) {
            numIndsHere++;
            if (!genotype.isNoCall()) {
                //increment stats for heterozygosity
                if (genotype.isHet()) {
                    numHetsHere++;
                }

                numGenosHere++;
                //increment stats for pairwise mismatches

                for (Allele allele : genotype.getAlleles()) {
                    if (allele.isCalled()) {
                        String alleleString = allele.toString();
                        alleleCounts.putIfAbsent(alleleString, 0);
                        alleleCounts.put(alleleString, alleleCounts.get(alleleString) + 1);
                    }
                }
            }
        }
        if (numGenosHere > 0) {
            //only if have one called genotype at least
            this.numSites++;

            this.totalHet += numHetsHere / (double)numGenosHere;

            //compute based on num sites
            float harmonicFactor = 0;
            for (int i = 1; i <= numIndsHere; i++) {
                harmonicFactor += 1.0 / i;
            }
            this.thetaRegionNumSites += 1.0 / harmonicFactor;

            //now compute pairwise mismatches
            float numPairwise = 0;
            int numDiffs = 0;
            for (String allele1 : alleleCounts.keySet()) {
                int allele1Count = alleleCounts.get(allele1);

                for (String allele2 : alleleCounts.keySet()) {
                    if (allele1.compareTo(allele2) < 0) {
                        continue;
                    }
                    if (allele1 .compareTo(allele2) == 0) {
                        numPairwise += allele1Count * (allele1Count - 1) * .5;

                    }
                    else {
                        int allele2Count = alleleCounts.get(allele2);
                        numPairwise += allele1Count * allele2Count;
                        numDiffs += allele1Count * allele2Count;
                    }
                }
            }

            if (numPairwise > 0) {
                this.totalAvgDiffs += numDiffs / numPairwise;
            }
        }
    }

    @Override
    public void finalizeEvaluation() {

        if (this.numSites > 0) {

            this.avgHet = this.totalHet / this.numSites;
            this.avgAvgDiffs = this.totalAvgDiffs / this.numSites;

        }
    }
}