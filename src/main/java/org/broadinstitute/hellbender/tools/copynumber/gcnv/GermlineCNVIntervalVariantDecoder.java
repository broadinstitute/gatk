package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Genotype;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledContigPloidyCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyNumberPosteriorDistribution;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntervalCopyNumberGenotypingData;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Decodes individual genotypes from a gCNV intervals VCF created with {@link GermlineCNVIntervalVariantComposer}.
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
public final class GermlineCNVIntervalVariantDecoder {

    final Map<String,Map<String,Integer>> sampleContigPloidyMap;

    public GermlineCNVIntervalVariantDecoder(final Collection<CalledContigPloidyCollection> contigPloidyCollections) {
        this.sampleContigPloidyMap = contigPloidyCollections.stream()
                .collect(Collectors.toMap(p -> p.getMetadata().getSampleName(), p -> getContigToPloidyCallMap(p)));
    }

    private Map<String,Integer> getContigToPloidyCallMap(final CalledContigPloidyCollection contigPloidyCollection) {
        return contigPloidyCollection.getRecords().stream().collect(Collectors.toMap(p -> p.getContig(), p -> p.getPloidy()));
    }

    public IntervalCopyNumberGenotypingData parseGenotype(final Genotype genotype, final SimpleInterval interval) {
        if (!genotype.hasExtendedAttribute(GermlineCNVIntervalVariantComposer.CNLP)) {
            throw new UserException.BadInput("Copy number genotype not found in genotype: " + genotype.toString());
        }
        final String[] phredScaledLikelihoodStrings = ((String)genotype.getExtendedAttribute(GermlineCNVIntervalVariantComposer.CNLP)).split(",");

        //Posteriors reported as integer phred-scaled likelihoods and need to be renormalized
        final double[] approximatePosteriors = new double[phredScaledLikelihoodStrings.length];
        double total = 0;
        for (int i = 0; i < phredScaledLikelihoodStrings.length; i++) {
            final double logLikelihood = -Double.valueOf(phredScaledLikelihoodStrings[i]) / (10.0 * Math.log(10.0));
            approximatePosteriors[i] = Math.max(Math.exp(logLikelihood), Double.MIN_VALUE);
            total += approximatePosteriors[i];
        }

        final Map<IntegerCopyNumberState,Double> copyNumberPosteriors = new HashMap<>(SVUtils.hashMapCapacity(phredScaledLikelihoodStrings.length));
        for (int i = 0; i < phredScaledLikelihoodStrings.length; i++) {
            copyNumberPosteriors.put(new IntegerCopyNumberState(i), FastMath.log(Math.max(approximatePosteriors[i] / total, Double.MIN_VALUE)));
        }

        final String sample = genotype.getSampleName();
        return new IntervalCopyNumberGenotypingData(
                interval,
                new CopyNumberPosteriorDistribution(copyNumberPosteriors),
                getBaselineState(sample, interval.getContig()));
    }

    private IntegerCopyNumberState getBaselineState(final String sample, final String contig) {
        if (!sampleContigPloidyMap.containsKey(sample)) {
            throw new IllegalArgumentException("Sample " + sample + " not found in ploidy contig calls.");
        }
        final Map<String,Integer> contigPloidyMap = sampleContigPloidyMap.get(sample);
        if (!contigPloidyMap.containsKey(contig)) {
            throw new IllegalArgumentException("Contig '" + contig + "' not found in ploidy contig calls for sample '" + sample + "'.");
        }
        return new IntegerCopyNumberState(contigPloidyMap.get(contig));
    }

    public Map<String,Map<String,Integer>> getSampleContigPloidyMap() {
        return sampleContigPloidyMap;
    }
}
