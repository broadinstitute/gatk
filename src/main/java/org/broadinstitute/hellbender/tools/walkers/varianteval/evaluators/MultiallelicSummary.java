package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

@Analysis(description = "Evaluation summary for multi-allelic variants")
public class MultiallelicSummary extends VariantEvaluator implements StandardEval {
    final protected static Logger logger = LogManager.getLogger(MultiallelicSummary.class);

    public enum Type {
        SNP, INDEL
    }

    public MultiallelicSummary(VariantEvalEngine engine) {
        super(engine);
    }

    // basic counts on various rates found
    @DataPoint(description = "Number of processed loci", format = "%d")
    public long nProcessedLoci = 0;

    @DataPoint(description = "Number of SNPs", format = "%d")
    public int nSNPs = 0;
    @DataPoint(description = "Number of multi-allelic SNPs", format = "%d")
    public int nMultiSNPs = 0;
    @DataPoint(description = "% processed sites that are multi-allelic SNPs", format = "%.5f")
    public double processedMultiSnpRatio = 0;
    @DataPoint(description = "% SNP sites that are multi-allelic", format = "%.3f")
    public double variantMultiSnpRatio = 0;

    @DataPoint(description = "Number of Indels", format = "%d")
    public int nIndels = 0;
    @DataPoint(description = "Number of multi-allelic Indels", format = "%d")
    public int nMultiIndels = 0;
    @DataPoint(description = "% processed sites that are multi-allelic Indels", format = "%.5f")
    public double processedMultiIndelRatio = 0;
    @DataPoint(description = "% Indel sites that are multi-allelic", format = "%.3f")
    public double variantMultiIndelRatio = 0;

    @DataPoint(description = "Number of Transitions", format = "%d")
    public int nTi = 0;
    @DataPoint(description = "Number of Transversions", format = "%d")
    public int nTv = 0;
    @DataPoint(description = "Overall TiTv ratio", format = "%.2f")
    public double TiTvRatio = 0;

    @DataPoint(description = "Multi-allelic SNPs partially known", format = "%d")
    public int knownSNPsPartial = 0;
    @DataPoint(description = "Multi-allelic SNPs completely known", format = "%d")
    public int knownSNPsComplete = 0;
    @DataPoint(description = "Multi-allelic SNP Novelty Rate")
    public String SNPNoveltyRate = "NA";

    //TODO -- implement me
    //@DataPoint(description = "Multi-allelic Indels partially known", format = "%d")
    public int knownIndelsPartial = 0;
    //@DataPoint(description = "Multi-allelic Indels completely known", format = "%d")
    public int knownIndelsComplete = 0;
    //@DataPoint(description = "Multi-allelic Indel Novelty Rate")
    public String indelNoveltyRate = "NA";


    @Override public int getComparisonOrder() { return 2; }

    @Override
    public void update2(final VariantContext eval, final VariantContext comp, final VariantEvalContext context) {
        if ( eval == null || (getEngine().getVariantEvalArgs().ignoreAC0Sites() && eval.isMonomorphicInSamples()) )
            return;

        // update counts
        switch ( eval.getType() ) {
            case SNP:
                nSNPs++;
                if ( !eval.isBiallelic() ) {
                    nMultiSNPs++;
                    calculatePairwiseTiTv(eval);
                    calculateSNPPairwiseNovelty(eval, comp);
                }
                break;
            case INDEL:
                nIndels++;
                if ( !eval.isBiallelic() ) {
                    nMultiIndels++;
                    calculateIndelPairwiseNovelty(eval, comp);
                }
                break;
            default:
                //throw new UserException.BadInput("Unexpected variant context type: " + eval);
                break;
        }

        return;
    }

    private void calculatePairwiseTiTv(VariantContext vc) {
        if (!vc.isSNP())
            return;

        for ( Allele alt : vc.getAlternateAlleles() ) {
            if (BaseUtils.SNPSubstitutionType(vc.getReference().getBases()[0], alt.getBases()[0]) == BaseUtils.BaseSubstitutionType.TRANSITION)
                nTi++;
            else
                nTv++;
        }
    }

    private void calculateSNPPairwiseNovelty(VariantContext eval, VariantContext comp) {
        if ( comp == null )
            return;

        int knownAlleles = 0;
        for ( Allele alt : eval.getAlternateAlleles() ) {
            if ( comp.getAlternateAlleles().contains(alt) )
                knownAlleles++;
        }

        if ( knownAlleles == eval.getAlternateAlleles().size() )
            knownSNPsComplete++;
        else if ( knownAlleles > 0 )
            knownSNPsPartial++;
    }

    private void calculateIndelPairwiseNovelty(VariantContext eval, VariantContext comp) {
        // TODO -- implement me
    }

    @Override
    public void finalizeEvaluation() {
        nProcessedLoci = getEngine().getnProcessedLoci();
        processedMultiSnpRatio = (double)nMultiSNPs / (double)nProcessedLoci;
        variantMultiSnpRatio = (double)nMultiSNPs / (double)nSNPs;
        processedMultiIndelRatio = (double)nMultiIndels / (double)nProcessedLoci;
        variantMultiIndelRatio = (double)nMultiIndels / (double)nIndels;

        TiTvRatio = (double)nTi / (double)nTv;

        final int all = nMultiSNPs;
        final int known = knownSNPsPartial + knownSNPsComplete;
        SNPNoveltyRate = Utils.formattedPercent(all - known, all);
        indelNoveltyRate = Utils.formattedPercent(all - known, all);
    }

    @Override
    public boolean requiresTerritoryToBeSpecified() {
        return true;
    }
}
