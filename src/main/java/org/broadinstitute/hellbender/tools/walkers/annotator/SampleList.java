package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * List samples that are non-reference at a given site
 *
 * <p>The output is a list of the samples that are genotyped as having one or more variant alleles. This allows you to easily determine which samples are non-reference (heterozygous or homozygous-variant) and compare them to samples that are homozygous-reference.</p>
 */
public final class SampleList extends InfoFieldAnnotation {

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( vc.isMonomorphicInSamples() || !vc.hasGenotypes() ) {
            return Collections.emptyMap();
        }

        final StringBuilder samples = new StringBuilder();
        for ( final Genotype genotype : vc.getGenotypesOrderedByName() ) {
            if ( genotype.isCalled() && !genotype.isHomRef() ){
                if ( samples.length() > 0 ) {
                    samples.append(",");
                }
                samples.append(genotype.getSampleName());
            }
        }

        if ( samples.length() == 0 ) {
            return Collections.emptyMap();
        }

        return Collections.singletonMap(getKeyNames().get(0), samples.toString());
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.SAMPLE_LIST_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() { return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}
