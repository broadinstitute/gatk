package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 *
 * Count of F1R2 alt reads at a variant for Mutect2 orientation artifact filter
 *
 * Created by tsato on 7/20/17.
 */
public class AltReadF1R2Count extends GenotypeAnnotation {
    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.ALT_READ_F1R2_COUNT_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(getKeyNames().get(0)));
    }

    @Override
    public void annotate( final ReferenceContext ref,
                          final VariantContext vc,
                          final Genotype g,
                          final GenotypeBuilder gb,
                          final ReadLikelihoods<Allele> likelihoods ) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        int altF1R2 = 0;
        // TODO: should this annotation be allele specific?
        Allele altAllele = vc.getAlternateAlleles().get(0);
        for (final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles(g.getSampleName())) {
            final GATKRead read = bestAllele.read;
            if (bestAllele.isInformative() && ReadUtils.readHasReasonableMQ(read) && read.isPaired()) {
                final Allele allele = bestAllele.allele;
                if (! allele.equals(altAllele, true)) {
                    continue;
                }

                if (readIsF1R2(read)){
                    altF1R2++;
                }
            }
        }

        gb.attribute(GATKVCFConstants.ALT_READ_F1R2_COUNT_KEY, altF1R2);
    }

    private boolean readIsF1R2(final GATKRead read){
        return (read.isFirstOfPair() && ! read.isReverseStrand()) || (read.isSecondOfPair() && read.isReverseStrand());
    }



}
