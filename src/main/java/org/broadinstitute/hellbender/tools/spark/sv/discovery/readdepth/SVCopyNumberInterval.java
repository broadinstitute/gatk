package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

public final class SVCopyNumberInterval {

    private final SVInterval interval;
    private final double[] copyNumberLogPosteriors;

    public SVCopyNumberInterval(final VariantContext variantContext, final SAMSequenceDictionary dictionary) {
        this.interval = new SVInterval(dictionary.getSequenceIndex(variantContext.getContig()), variantContext.getStart(), variantContext.getEnd() + 1);
        final GenotypesContext genotypesContext = variantContext.getGenotypes();
        if (genotypesContext.isEmpty()) {
            throw new UserException.BadInput("No genotypes found in variant context " + variantContext.getID());
        }
        final Genotype genotype = genotypesContext.get(0);
        if (!genotype.hasExtendedAttribute("CNLP")) {
            throw new UserException.BadInput("Copy number genotype not found in variant context " + variantContext.getID());
        }
        final String[] phredScaledLikelihoodStrings = ((String) genotype.getExtendedAttribute("CNLP")).split(",");

        //Posteriors reported as integer phred-scaled likelihoods and need to be renormalized
        final double[] approximatePosteriors = new double[phredScaledLikelihoodStrings.length];
        double total = 0;
        for (int i = 0; i < phredScaledLikelihoodStrings.length; i++) {
            final double logLikelihood = -Double.valueOf(phredScaledLikelihoodStrings[i]) / (10.0 * Math.log(10.0));
            approximatePosteriors[i] = Math.max(Math.exp(logLikelihood), Double.MIN_NORMAL);
            total += approximatePosteriors[i];
        }

        copyNumberLogPosteriors = new double[phredScaledLikelihoodStrings.length];
        for (int i = 0; i < phredScaledLikelihoodStrings.length; i++) {
            copyNumberLogPosteriors[i] = Math.min(Math.log(approximatePosteriors[i] / total), Double.MIN_NORMAL);
        }
    }

    public SVInterval getInterval() {
        return interval;
    }

    public double[] getCopyNumberLogPosteriorsArray() {
        return copyNumberLogPosteriors;
    }
}
