package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.spark_project.guava.collect.Range;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Potential polymorphic NuMT annotation compares the number of alts to the autosomal coverage
 *
 * <p>Polymorphic NuMTs (NuMT sequence that is not in the reference) will result in autosomal reads that map to the
 * mitochondria. This annotation notes when the number of alt reads falls within 80% of the autosomal coverage
 * distribution, assuming that autosomal coverage is a Poisson distribution given a central statistic of the depth
 * (median is recommended). This will also include true low allele fraction mitochondrial variants and should be used
 * as an annotation, rather than a filter.</p>
 *
 * <h3>Caveat</h3>
 * <p>This annotation can only be calculated in Mutect2 if median-autosomal-coverage argument is provided.</p>
 *
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of alts indicates it could be an autosomal false positive.")
public class PolymorphicNuMT extends GenotypeAnnotation implements Annotation {
    protected final OneShotLogger warning = new OneShotLogger(this.getClass());
    public static final String KEY = GATKVCFConstants.POTENTIAL_POLYMORPHIC_NUMT_KEY;
    private static final double LOWER_BOUND_PROB = .1;
    private Range<Long> autosomalHetRange;
    private Range<Long> autosomalHomAltRange;

    public PolymorphicNuMT(final double medianAutosomalCoverage){
        final PoissonDistribution autosomalCoverage = new PoissonDistribution(medianAutosomalCoverage);
        final long minAutosomalHomAlt = autosomalCoverage.inverseCumulativeProbability(LOWER_BOUND_PROB);
        final long maxAutosomalHomAlt = autosomalCoverage.inverseCumulativeProbability(1 - LOWER_BOUND_PROB);
        final long minAutosomalHet = minAutosomalHomAlt / 2;
        final long maxAutosomalHet = maxAutosomalHomAlt / 2;
        autosomalHetRange = Range.closed(minAutosomalHet, maxAutosomalHet);
        autosomalHomAltRange = Range.closed(minAutosomalHomAlt, maxAutosomalHomAlt);
    }

    // Barclay requires that each annotation define a constructor that takes no arguments
    public PolymorphicNuMT(){ }

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);
        final double[] lods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        if (lods==null) {
            warning.warn(String.format("One or more variant contexts is missing the 'TLOD' annotation, %s will not be computed for these VariantContexts", GATKVCFConstants.POTENTIAL_POLYMORPHIC_NUMT_KEY));
            return;
        }
        final int indexOfMaxLod = MathUtils.maxElementIndex(lods);
        final Allele altAlelle = vc.getAlternateAllele(indexOfMaxLod);
        Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        final long numAltReads = bestAlleles.stream().filter(ba -> ba.isInformative() && ba.allele.equals(altAlelle)).count();
        if ( autosomalHetRange.contains(numAltReads) || autosomalHomAltRange.contains(numAltReads) ) {
            gb.attribute(GATKVCFConstants.POTENTIAL_POLYMORPHIC_NUMT_KEY, "true");
        }
    }
    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(KEY));
    }
    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(KEY);
    }
}
